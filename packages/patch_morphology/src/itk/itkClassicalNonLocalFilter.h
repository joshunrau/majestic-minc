/* ----------------------------- MNI Header -----------------------------------
@COPYRIGHT  :
              Copyright 2014 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

              This program is free software: you can redistribute it and/or modify
              it under the terms of the GNU General Public License as published by
              the Free Software Foundation, either version 3 of the License, or
              (at your option) any later version.
---------------------------------------------------------------------------- */

#ifndef __itkClassicalNonLocalFilter_h
#define __itkClassicalNonLocalFilter_h

#include "itkNonLocalPatchesFilter.h"
#include "itkPatchDistance.h"
#include "itkPatchCostFunction.h"

namespace itk
{

  template<class TInputImage, class TOutputImage, class TSearch, class TPatch,
          class TDistance=L2PatchDistance<TInputImage,TPatch>,class TWeight=InvExpWeight<double>,
          class TPreselectionFilter=NOOPPreselection<3> >
  class  ClassicalNonLocalFilter : 
      public NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch,TPreselectionFilter>
  {
  public:
    /** Standard class typedefs. */
    typedef ClassicalNonLocalFilter Self;
    typedef NonLocalPatchesFilter<TInputImage, TOutputImage, TSearch, TPatch,TPreselectionFilter>
                                      Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;
    
    /** Standard New method. */
    itkNewMacro(Self);  

    /** Runtime information support. */
    itkTypeMacro(ClassicalNonLocalFilter, 
                NonLocalPatchesFilter);
    
    /** Declaration of pixel type. */
    typedef typename Superclass::PixelType PixelType;

    /** Kernel (structuring element) iterator. */
    typedef typename Superclass::SearchKernelIteratorType SearchKernelIteratorType;
    
    /** Kernel (structuring element) iterator. */
    typedef typename Superclass::PatchKernelIteratorType  PatchKernelIteratorType;

    /** Neighborhood iterator type. */
    typedef typename Superclass::NeighborhoodIteratorType NeighborhoodIteratorType;

    /** Kernel typedef. */
    typedef typename Superclass::SearchKernelType SearchKernelType;
    typedef typename Superclass::PatchKernelType  PatchKernelType;
    
    typedef typename Superclass::PreselectionFilterPointerType PreselectionFilterPointerType;
    typedef typename TDistance::Pointer DistancePointerType;

    /** Default boundary condition type */
    typedef typename Superclass::DefaultBoundaryConditionType DefaultBoundaryConditionType;

    /** ImageDimension constants */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);
    itkStaticConstMacro(KernelDimension, unsigned int,
                        TSearch::NeighborhoodDimension);

    /** Type of the pixels in the Kernel. */
    typedef typename TSearch::PixelType            KernelPixelType;

    /** Set squared sigma. */
    itkSetMacro(sigma2, double);
    /** Get squared sigma. */
    itkGetMacro(sigma2, double);

    /** Set weight functor. */
    itkSetMacro(Weight, TWeight);
    /** Get weight functor. */
    itkGetMacro(Weight, TWeight);

    /** Set distance functor. */
    itkSetMacro(Distance, DistancePointerType);
    /** Get distance functor. */
    itkGetMacro(Distance, DistancePointerType);

    /** set output weight flag */
    itkSetMacro(OutputMeanWeight,bool);
    /** get output weight flag */
    itkGetMacro(OutputMeanWeight,bool);
    
    

  // #ifdef ITK_USE_CONCEPT_CHECKING
  //   /** Begin concept checking */
  //   itkConceptMacro(InputConvertibleToOutputCheck,
  //     (Concept::Convertible<PixelType, typename TOutputImage::PixelType>));
  //   itkConceptMacro(SameDimensionCheck1,
  //      (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
  //   itkConceptMacro(SameDimensionCheck2,
  //     (Concept::SameDimension<InputImageDimension, KernelDimension>));
  //   itkConceptMacro(InputGreaterThanComparableCheck,
  //     (Concept::GreaterThanComparable<PixelType>));
  //   itkConceptMacro(KernelGreaterThanComparableCheck,
  //     (Concept::GreaterThanComparable<KernelPixelType>));
  //   /** End concept checking */
  // #endif

  protected:
    ClassicalNonLocalFilter();
    ~ClassicalNonLocalFilter() {};

    /** Evaluate image neighborhood with kernel to find the new value 
    * for the center pixel value
    *
    */
    
    PixelType Evaluate(const NeighborhoodIteratorType &searchIt,
                      const  NeighborhoodIteratorType &patchIt1,
                      NeighborhoodIteratorType &patchIt2,
                      const SearchKernelIteratorType searchKernelBegin,
                      const SearchKernelIteratorType searchKernelEnd,
                      const PatchKernelIteratorType patchKernelBegin,
                      const PatchKernelIteratorType patchKernelEnd,
                      PreselectionFilterPointerType flt
                      );

    
    void PrintSelf(std::ostream& os, itk::Indent indent) const;

  private:
    ClassicalNonLocalFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
    
    double              m_sigma2;
    DistancePointerType m_Distance;
    
    TWeight             m_Weight;
    bool                m_OutputMeanWeight;

  }; // end of class

} // end namespace itk

#include "itkClassicalNonLocalFilter.txx"

#endif //__itkClassicalNonLocalFilter_h

// kate: space-indent on; indent-width 2; indent-mode C++;replace-tabs on;word-wrap-column 80;show-tabs on;tab-width 2