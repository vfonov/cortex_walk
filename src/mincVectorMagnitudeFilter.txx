#ifndef __itkVectorMagnitudeFilter_txx
#define __itkVectorMagnitudeFilter_txx

#include "mincVectorMagnitudeFilter.h"

#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace minc
{

template <typename TInputImage, typename TOutputImage>
VectorMagnitudeFilter<TInputImage,  TOutputImage>
::VectorMagnitudeFilter()
{
  m_RequestedNumberOfThreads = this->GetNumberOfThreads();
}


template< typename TInputImage,  typename TOutputImage >
void
VectorMagnitudeFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  Superclass::BeforeThreadedGenerateData();
}
#if (ITK_VERSION_MAJOR==3)
template< typename TInputImage,  typename TOutputImage >
void
VectorMagnitudeFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
#else
template< typename TInputImage,  typename TOutputImage >
void
VectorMagnitudeFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       itk::ThreadIdType threadId)
#endif
{

  itk::ImageRegionIterator<TOutputImage> it(this->GetOutput(), outputRegionForThread);
  itk::ImageRegionConstIterator<TInputImage> iit(this->GetInput(), outputRegionForThread);
  
  // Support progress methods/callbacks
  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  while ( ! it.IsAtEnd() )
  {
    double sum=0;
    for(unsigned int i=0;i<VectorDimension;i++)
      sum+=(double)iit.Get()[i]*(double)iit.Get()[i];
    
    it.Set( sqrt(sum) );
    
    ++iit;
    ++it;
    progress.CompletedPixel();
  }
}


template <typename TInputImage, typename TOutputImage>
void
VectorMagnitudeFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "m_RequestedNumberOfThreads = " << m_RequestedNumberOfThreads
     << std::endl;
}

} // end namespace itk

#endif
