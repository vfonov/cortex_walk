/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2018,2006 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
//use for compatibility with older version of libminc
#define VIO_PREFIX_NAMES 1

#include "minc_wrappers.h"
#include "itkCommand.h"

#include "mincVectorMagnitudeFilter.h"
#include "mincVectorBSplineInterpolate.h"


#include <getopt.h>
#include  <bicpl.h>
#include <unistd.h>

#include <iostream>
#include <fstream>

using namespace  std;
using namespace  minc;

typedef minc::mincVectorBSplineInterpolate<minc::def3d,double> Vector_Interpolator;
typedef minc::VectorMagnitudeFilter<minc::def3d,minc::image3d> _MagFilterType;


class CommandProgressUpdate : public itk::Command
{
  public:
    typedef CommandProgressUpdate   Self;
    typedef itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer;
    itkNewMacro( Self );
  protected:
    CommandProgressUpdate() {};
  public:
    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      const itk::ProcessObject * filter =
          dynamic_cast< const itk::ProcessObject * >( object );
      if( ! itk::ProgressEvent().CheckEvent( &event ) )
      {
        return;
      }
      std::cout.precision(3);
      std::cout << filter->GetProgress()*100 << "%\t" <<std::flush;
    }
};


void show_usage (const char * prog)
{
  std::cerr<<"Usage:"<<prog<<" in_grid1.mnc in.obj out.obj "<<std::endl
      <<"[ "<<std::endl
      <<"  --distance <r> distance to travel, default 10.0"<<std::endl
      <<"  --distfile <f> distances file"<<std::endl
      <<"  --threshold <r> gradient magnitude threshold"<<std::endl
      <<"  --clobber clobber output file(s)" <<std::endl
      <<"  --step <f> step size, default 1.0" << std::endl
      <<"]"<<std::endl; 
}

int main (int argc, char **argv)
{
  int verbose=0,clobber=0;
  int glog=0,gexp=0,gsqrt=0,gsq=0,ginv=0;
  int iter=10;
  int gmag=0,gdet=0;
  int gcompat=1;
  int gcond=0,gscond=0; 
  int order=1;
  
  double mag_threshold=0.1;
  double distance=10.0;
  double step=1.0;
  std::string distfile;
  ifstream ifs;
  
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"distance", required_argument, 0, 'd'},
    {"distfile", required_argument, 0, 'D'},
    {"threshold",required_argument,0, 't'},
    {"help",     no_argument,0, 'h'},
    {0, 0, 0, 0}
  };
  
  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "hvd:t:s:D:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c)
    {
      case 0:
        break;
      case 'd': 
        distance=atof(optarg);
        break;
      case 'D': 
        distance=0;
        distfile=optarg;
        break;
      case 't':
        mag_threshold=atof(optarg);
        break;
      case 's':
        step=atof(optarg);
        break;
      case 'v':
        cout << "Version: 0.1" << endl;
        return 0;
      case '?':
        /* getopt_long already printed an error message. */
      case 'h':
         show_usage (argv[0]);
         return 0;
      default:
        show_usage (argv[0]);
        return 1;
    }
  }

  if((argc - optind) < 3) {
    show_usage (argv[0]);
    return 1;
  }
  std::string in_grid=argv[optind];
  std::string in_obj=argv[optind+1];
  std::string out_objf=argv[optind+2];
  
  if (!clobber && !access(out_objf.c_str(), F_OK))
  {
    cerr << out_objf.c_str() << " Exists!" << endl;
    return 1;
  }
  
  try
  {
    
    CommandProgressUpdate::Pointer observer = CommandProgressUpdate::New();
  
    //load gradient field
    std::cout<<"Loading:"<<in_grid.c_str()<<" ..."<<std::endl;
    
    def3d::Pointer  grid=minc::load_minc<def3d>(in_grid.c_str());
    
    std::cout<<"Building BSpline interpolator ..."<<std::endl;
    Vector_Interpolator::Pointer interpolator(Vector_Interpolator::New());
    interpolator->SetSplineOrder(order);
    interpolator->SetInputImage(grid);
    
    VIO_File_formats         format;
    object_struct        **object_list;
    int n_objects=0;
    
    //bicpl sucks!
    if( input_graphics_file( (char*)in_obj.c_str(), &format, &n_objects,
        &object_list ) != VIO_OK )
    {
      std::cerr << " Error reading "<<in_obj.c_str() << std::endl;
      return 1;
    }
    
    std::vector<minc::tag_point> lines;
    std::vector<int> line_index;
    
    std::cout<<"Processing "<<n_objects<<" objects"<<std::endl;
    if(!distfile.empty()) {
        std::cout<<"Tracking "<<distfile.c_str()<<" "<<std::endl;
        
        ifs.open(distfile.c_str());
        
    } else {
        std::cout<<"Tracking "<<distance<<" mm "<<std::endl;
    }
    
    std::cout<<"Gradient magnitude threshold "<< mag_threshold <<std::endl;
    std::cout<<"Step length "<< step <<std::endl;
    
    for(int i=0;i< n_objects;i++ )
    {
      if( get_object_type( object_list[i] ) == POLYGONS )
      {
        polygons_struct      *polygons;
        polygons = get_polygons_ptr(object_list[i]);
        
        double avg_distance=0.0;
        
        for(int pnt=0;pnt<polygons->n_points;pnt++ )
        {
          VIO_Point p_=polygons->points[pnt];//seeding point 
          minc::tag_point p_orig,p;
          p_orig[0]=p_.coords[0];p_orig[1]=p_.coords[1];p_orig[2]=p_.coords[2];
          //line_index.push_back(lines.size());
          //lines.push_back(p_orig);
          p=p_orig;

          if(!distfile.empty()) ifs>>distance;
          
          
          for(double i=0.0;i<distance;)
          {
            vnl_vector<double> _in(3);
            _in.set(interpolator->Evaluate(p).GetDataPointer());
            double _mag=sqrt(_in[0]*_in[0]+_in[1]*_in[1]+_in[2]*_in[2]);
            
            if(_mag<mag_threshold)
              break; 
            //_in.normalize();
            
            //walk against the gradient
            
            p[0]-=_in[0]*step/_mag;
            p[1]-=_in[1]*step/_mag;
            p[2]-=_in[2]*step/_mag;
            
            i+=step;
          }
          
          avg_distance+=p_orig.EuclideanDistanceTo(p);
          
          polygons->points[pnt].coords[0]=p[0];
          polygons->points[pnt].coords[1]=p[1];
          polygons->points[pnt].coords[2]=p[2];
          //int size = GET_OBJECT_SIZE( *polygons, poly );
          //if(size<3) continue; //?
          
/*          p1=polygons->points[POINT_INDEX( polygons->end_indices, poly, 0 )];
          p2=polygons->points[POINT_INDEX( polygons->end_indices, poly, 1 )];
          p3=polygons->points[POINT_INDEX( polygons->end_indices, poly, 2 )];*/
          /*
          std::cout<<poly<<"\t"<<POINT_INDEX( polygons->end_indices, poly, 0 )<<" "
              <<POINT_INDEX( polygons->end_indices, poly, 1 )<<" "
              <<POINT_INDEX( polygons->end_indices, poly, 2 )<<std::endl;*/
          //std::cout<<poly<<"\t"<<p1.coords[0]<<" "<<p1.coords[1]<<" "<<p1.coords[2]<<std::endl;
          //center
        }
        avg_distance/= polygons->n_points;
        std::cout<<"Average distance "<<avg_distance<<std::endl;
      }
    } 
    int status = output_graphics_file( (char*)out_objf.c_str(), format,n_objects, object_list );
    //free up memory
    delete_object_list( n_objects, object_list );
    return( status != VIO_OK );
    /*
    //outputting lines
    FILE *fp=fopen(out_objf.c_str(),"w");

  // IL want to write in binary

    fprintf(fp,"L 0.1 %d\n",lines.size());

  //WRITING cordinates

    for (int i=0;i<lines.size();i++) {
      fprintf(fp,"%f %f %f\n",lines[i][0],lines[i][1],lines[i][2]);
    }

    fprintf(fp,"%d\n1 ",line_index.size()); //1 color per line

    for (int i=0;i<line_index.size();i++) { //random colors
      fprintf(fp,"%f %f %f 1.0\n",rand()/(double)RAND_MAX,rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);
    }

    for (int i=0;i<line_index.size();i++) {
      
      if(i<(line_index.size()-1))
        fprintf(fp,"%d ",line_index[i+1]);
      else
        fprintf(fp,"%d ",lines.size());//last line lasts untill the end
    }

    fprintf(fp,"\n");

    for (int i=0;i<line_index.size();i++) 
    {
      int end=i<(line_index.size()-1)?line_index[i+1]:lines.size();
      for(int j=line_index[i];j<end;j++)
        fprintf(fp," %d",j);
      
      fprintf(fp,"\n");
    }
    fclose(fp);*/
    
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  }
  return 0;
};
