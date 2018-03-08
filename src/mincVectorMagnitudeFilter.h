#ifndef __itkVectorMagnitudeFilter_h
#define __itkVectorMagnitudeFilter_h

#include "itkImageToImageFilter.h"

#if ( ITK_VERSION_MAJOR > 3 ) 
#include "itkImage.h"
#else
#include "itkOrientedImage.h"
#endif 

#include "itkVector.h"

namespace minc
{
  /** \class VectorMagnitudeFilter
  *
  * \brief Computes a scalar image from a vector image (e.g., deformation field)
  * input, where each output scalar at each pixel is the magnitude
  * of the vector field at that location.  
  *
  * \par Template Parameters (Input and Output)
  * This filter has one required template parameter which defines the input
  * image type.  The pixel type of the input image is assumed to be a vector
  * (e.g., itk::Vector, itk::RGBPixel, itk::FixedArray).  The scalar type of the
  * vector components must be castable to floating point.  Instantiating with an
  * image of RGBPixel<unsigned short>, for example, is allowed, but the filter
  * will convert it to an image of Vector<float,3> for processing.
  *
  * The second template parameter, TRealType, can be optionally specified to
  * define the scalar numerical type used in calculations.  This is the
  * component type of the output image, which will be of
  * itk::Vector<TRealType, N>, where N is the number of channels in the multiple
  * component input image.  The default type of TRealType is float.  For extra
  * precision, you may safely change this parameter to double.
  *
  * The third template parameter is the output image type.  The third parameter
  * will be automatically constructed from the first and second parameters, so
  * it is not necessary (or advisable) to set this parameter explicitly.  Given
  * an M-channel input image with dimensionality N, and a numerical type
  * specified as TRealType, the output image will be of type
  * itk::Image<TRealType, N>.
  *
  * \par Constraints
  *
  * The template parameter TRealType must be floating point (float or double) or
  * a user-defined "real" numerical type with arithmetic operations defined
  * sufficient to compute derivatives.
  *
  * \sa Image
  *
  * \author Vladimir S. Fonov
  */
  template < typename TInputImage,
             typename TOutputImage >
  class VectorMagnitudeFilter :
      public itk::ImageToImageFilter< TInputImage, TOutputImage >
  {
  public:
    /** Standard class typedefs. */
    typedef VectorMagnitudeFilter      Self;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self>                              Pointer;
    typedef itk::SmartPointer<const Self>                        ConstPointer;
    /** Method for creation through the object factory. */
    itkNewMacro(Self);
  
    /** Run-time type information (and related methods) */
    itkTypeMacro(VectorMagnitudeFilter, ImageToImageFilter);
  
    /** Extract some information from the image types.  Dimensionality
    * of the two images is assumed to be the same. */
    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef typename TInputImage::PixelType  InputPixelType;
  
    /** Image typedef support */
    typedef TInputImage                       InputImageType;
    typedef TOutputImage                      OutputImageType;
    typedef typename InputImageType::Pointer  InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;
  
    /** The dimensionality of the input and output images. */
    itkStaticConstMacro(ImageDimension, unsigned int,
                        TOutputImage::ImageDimension);
  
    /** Length of the vector pixel type of the input image. */
    itkStaticConstMacro(VectorDimension, unsigned int,
                        InputPixelType::Dimension);
  
    /** Define the data type and the vector of data type used in calculations. */
    typedef OutputPixelType RealType;
  
    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  
  
  protected:
    VectorMagnitudeFilter();
    virtual ~VectorMagnitudeFilter() {}
  
    /** Do any necessary casting/copying of the input data.  Input pixel types
      whose value types are not real number types must be cast to real number
      types. */
    void BeforeThreadedGenerateData ();
  
    /** VectorMagnitudeFilter can be implemented as a
    * multithreaded filter Therefore, this implementation provides a
    * ThreadedGenerateData() routine which is called for each
    * processing thread. The output image data is allocated
    * automatically by the superclass prior to calling
    * ThreadedGenerateData().  ThreadedGenerateData can only write to
    * the portion of the output image specified by the parameter
    * "outputRegionForThread"
    *
    * \sa ImageToImageFilter::ThreadedGenerateData(),
    *     ImageToImageFilter::GenerateData() */
#if (ITK_VERSION_MAJOR==3)
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                              int threadId );
#else
    void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                               itk::ThreadIdType threadId);    
#endif
    void PrintSelf(std::ostream& os, itk::Indent indent) const;
  
    typedef typename InputImageType::Superclass ImageBaseType;
  
    /** Get access to the input image casted as real pixel values */
    itkGetConstObjectMacro( RealValuedInputImage, ImageBaseType );
  
 
  private:
    int  m_RequestedNumberOfThreads;
  
    typename ImageBaseType::ConstPointer m_RealValuedInputImage;
  
    VectorMagnitudeFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "mincVectorMagnitudeFilter.txx"
#endif

#endif
