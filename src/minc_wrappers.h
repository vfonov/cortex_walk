#ifndef _MINC_WRAPPERS_H_
#define _MINC_WRAPPERS_H_

#include <complex>
#include <vector>
#include <algorithm>
#include <itkArray.h>
#include <iostream>


#include <stdlib.h>

#include "itkArray.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIOFactory.h"

#include <minc_io_exceptions.h>
#include <minc_io_fixed_vector.h>

#ifndef ITK_REPORT_ERROR
#define ITK_REPORT_ERROR(MSG) throw itk::ExceptionObject(__FILE__,__LINE__,MSG)
#endif //ITK_REPORT_ERROR


namespace minc
{
    //! default label voxel type
    typedef unsigned char minc_mask_voxel;

    //! default minc file voxel type
    typedef float voxel_type;

    //! default minc volume dimension
    const int volume_dimensions = 3;

    //! default minc grid volume voxel type
    typedef itk::Vector<float,volume_dimensions>    def_vector;
    //! minc tag point
    typedef itk::Point<double,volume_dimensions>    tag_point;
    typedef std::vector<tag_point>  tag_points;
    
    //! default minc complex voxel type
    typedef std::complex < voxel_type > complex;

    typedef itk::Image < complex, volume_dimensions >         image3d_complex;
    typedef itk::Image < voxel_type, volume_dimensions >      image3d;
    typedef itk::Image < minc_mask_voxel, volume_dimensions > mask3d;
    typedef itk::Image < def_vector, volume_dimensions >      def3d;

    
    typedef itk::ImageRegionIteratorWithIndex < image3d > image3d_iterator;
    typedef itk::ImageRegionConstIteratorWithIndex < image3d > image3d_const_iterator;
    
    typedef itk::ImageRegionIteratorWithIndex < mask3d > mask3d_iterator;
    typedef itk::ImageRegionConstIteratorWithIndex < mask3d > mask3d_const_iterator;
    
    typedef itk::ImageRegionIteratorWithIndex < def3d > def3d_iterator;
    typedef itk::ImageRegionConstIteratorWithIndex < def3d > def3d_const_iterator;
    
    typedef itk::ImageRegionIteratorWithIndex < image3d_complex > image3d_complex_iterator;
    typedef itk::ImageRegionConstIteratorWithIndex < image3d_complex > image3d_complex_const_iterator;

    //! find a maximum of elements
    template<class T> float v_max(const T& c)
    {
    float s=std::numeric_limits < float >::min ();;
    for(unsigned int i=0;i<3;i++)
        if(c[i]>s) s=c[i];
    return s;
    }

    //! find a minimum of elements
    template<class T> float v_min(const T& c) 
    {
    float s=std::numeric_limits < float >::max ();;
    for(unsigned int i=0;i<3;i++)
        if(c[i]<s) s=c[i];
    return s;
    }

    //! allocate volume of the same dimension,spacing and origin
    template<class T,class S> void allocate_same(typename T::Pointer &image,const typename S::Pointer &sample)
    {
                image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
                image->SetBufferedRegion(sample->GetLargestPossibleRegion());
                image->SetRequestedRegion(sample->GetLargestPossibleRegion());
                image->SetSpacing( sample->GetSpacing() );
                image->SetOrigin ( sample->GetOrigin() );
    image->SetDirection(sample->GetDirection());
                image->Allocate();
    }

    //! allocate volume of the same dimension,spacing and origin
    template<class T,class S> void allocate_same(typename T::Pointer &image,const typename S::ConstPointer &sample)
    {
    image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
    image->SetBufferedRegion(sample->GetLargestPossibleRegion());
    image->SetRequestedRegion(sample->GetLargestPossibleRegion());
    image->SetSpacing( sample->GetSpacing() );
    image->SetOrigin ( sample->GetOrigin() );
    image->SetDirection(sample->GetDirection());
    image->Allocate();
    }

    //! allocate volume of the same dimension,spacing and origin
    template<class T,class S> typename T::Pointer allocate_same(const typename S::Pointer &sample)
    {
    typename T::Pointer image=T::New();
    image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
    image->SetBufferedRegion(sample->GetLargestPossibleRegion());
    image->SetRequestedRegion(sample->GetLargestPossibleRegion());
    image->SetSpacing( sample->GetSpacing() );
    image->SetOrigin ( sample->GetOrigin() );
    image->SetDirection(sample->GetDirection());
    image->Allocate();
    return image;
    }


    //! allocate volume of the same dimension,spacing and origin and do nearest neighbour resampling
    template<class T,class S,class E,class D> void nearest_resample_like(typename T::Pointer &dst,const typename S::Pointer &sample,const typename E::Pointer &src, const D& def)
    {
    allocate_same<T,S>(dst,sample);
                itk::ImageRegionIteratorWithIndex<T> it(dst, dst->GetLargestPossibleRegion());
                for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        tag_point p;
        typename E::IndexType idx;
        dst->TransformIndexToPhysicalPoint(it.GetIndex(),p);
        if(src->TransformPhysicalPointToIndex(p,idx))
        {
        it.Value()=src->GetPixel(idx);
        }else{
        it.Value()=def;
        }
    }
    }

    //! check if volumes have the same dimensions, spacing and origin
    template<class T,class S> bool check_same(typename T::Pointer image,typename S::Pointer sample)
    {
    return
        (image->GetLargestPossibleRegion() == sample->GetLargestPossibleRegion()) &&
        (image->GetSpacing() == sample->GetSpacing()) &&
        (image->GetOrigin() == sample->GetOrigin()) &&
        (image->GetDirection().GetVnlMatrix()  == sample->GetDirection().GetVnlMatrix()); // carefull here! , maybe we should calculate some kind of difference here ?
        // this is warkaround a bug in itk
    }

    //! allocate volume
    //! \param[out] image - volume to allocate
    //! \param dims - dimensions (voxels)
    //! \param spacing - volume spacing (mm)
    //! \param origin  - volume origin (mm)
    template<class T> void allocate_image3d(typename T::Pointer &image, 
        const itk::Array<unsigned int> &dims, 
        const itk::Array<double>& spacing, 
        const itk::Array<double>& origin)
    {
    typename T::SizeType  imageSize3D = {{ dims[0], dims[1], dims[2]}};
    typename T::IndexType startIndex3D = { {0, 0, 0}};
    typename T::RegionType region;
    region.SetSize  (imageSize3D);
    region.SetIndex (startIndex3D);
    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);
    image->SetSpacing( spacing );
    image->SetOrigin( origin );
    image->Allocate ();
    }

    //! allocate volume
    //! \param[out] image - volume to allocate
    //! \param dims - dimensions (voxels)
    //! \param spacing - volume spacing (mm)
    //! \param origin  - volume origin (mm)
    template<class T> void allocate_image3d(typename T::Pointer &image, 
        const itk::Vector<unsigned int,3> &dims, 
        const itk::Vector<double,3>& spacing, 
        const itk::Vector<double,3>& origin)
    {
    typename T::SizeType  imageSize3D = {{ dims[0], dims[1], dims[2]}};
    typename T::IndexType startIndex3D = { {0, 0, 0}};
    typename T::RegionType region;
    region.SetSize  (imageSize3D);
    region.SetIndex (startIndex3D);
    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);
    image->SetSpacing( spacing );
    image->SetOrigin( origin.GetDataPointer () );
    image->Allocate ();
    }

    //! allocate volume
    //! \param[out] image - volume to allocate
    //! \param dims - dimensions (voxels)
    //! \param spacing - volume spacing (mm)
    //! \param origin  - volume origin (mm)
    template<class T> void allocate_image3d(typename T::Pointer &image, 
        const fixed_vec<3, unsigned int>&dims, 
        const fixed_vec<3, double>& spacing=fixed_vec<3, double>(1.0) , 
        const fixed_vec<3, double>& origin=fixed_vec<3, double>(0.0))
    {
    typename T::SizeType  imageSize3D = {{ dims[0], dims[1], dims[2]}};
    typename T::IndexType startIndex3D = { {0, 0, 0}};
    typename T::RegionType region;
    region.SetSize  (imageSize3D);
    region.SetIndex (startIndex3D);
    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);
    image->SetSpacing( spacing.c_buf() );
    image->SetOrigin( origin.c_buf() );
    image->Allocate ();
    }

    //! allocate volume
    //! \param[out] image - volume to allocate
    //! \param dims - dimensions (voxels)
    //! \param spacing - volume spacing (mm)
    //! \param origin  - volume origin (mm)
    template<class T> void allocate_image3d(typename T::Pointer &image, 
        const itk::Size<3> &dims)
    {
    //typename T::SizeType  imageSize3D = {{ dims[0], dims[1], dims[2]}};
    typename T::IndexType startIndex3D = { {0, 0, 0}};
    typename T::RegionType region;
    double spacing[3]= { 1.0,1.0,1.0};
    double origin[3]= { 0.0,0.0,0.0};
    region.SetSize  (dims);
    region.SetIndex (startIndex3D);
    image->SetLargestPossibleRegion (region);
    image->SetBufferedRegion (region);
    image->SetRequestedRegion (region);

    image->SetSpacing( spacing );
    image->SetOrigin( origin );
    image->Allocate ();
    }


    inline image3d::SizeType operator/= (image3d::SizeType & s, int d)
    {
    s[0] /= d;
    s[1] /= d;
    s[2] /= d;
    return s;
    }

    inline image3d::SizeType operator*= (image3d::SizeType & s, int d)
    {
    s[0] *= d;
    s[1] *= d;
    s[2] *= d;
    return s;
    }

    //! a helper function for minc reading
    template <class T> typename T::Pointer load_minc(const char *file)
    {
        
    typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();

    reader->SetFileName(file);
    reader->Update();

    return reader->GetOutput();
    }

    //! a helper function for minc writing
    template <class T> void save_minc(const char *file,typename T::Pointer img)
    {
    typename itk::ImageFileWriter< T >::Pointer writer = itk::ImageFileWriter<T>::New();
    writer->SetFileName(file);
    writer->SetInput( img );

    if( getenv("MINC_COMPRESS") != NULL)
        writer->SetUseCompression( true );

    writer->Update();
    } 

    //! minc4itk compatibility function
    template <class T> void  load_minc(const char *file,typename T::Pointer& img)
    {
    img=load_minc<T>(file);
    }

    template <class T> void  imitate_minc(const char *file,typename T::Pointer& img)
    {
    typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();

    reader->SetFileName(file);
    reader->Update();

    typename T::Pointer sample=reader->GetOutput();
    allocate_same<T,T>(img,sample);
    }


    
    
  //! minc4itk compatibility function
  template <class T> void  load_minc(const char *file,T &img)
  {
    img=load_minc<typename T::ObjectType>(file);
  }
};//minc
#endif //_MINC_WRAPPERS_H_
