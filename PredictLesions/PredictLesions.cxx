/*=========================================================================

  Program:   Lesion Segmentation CLIs for Slicer4
  Module:    $HeadURL$
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Biomedical Mining, LLC. All rights reserved.
  See https://github.com/msscully/LesionSegmentation for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>
#include <vector>
#include <utility>
#include "itkPluginUtilities.h" 
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"
#include "itkPluginFilterWatcher.h"
#include "itkImageToVectorImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "nanoflann.hpp"
#include "LesionSegmentationModel.h"
#include "PredictLesionsCLP.h"

/* Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
*/ 
namespace
{

typedef short      PixelType;
const unsigned int Dimension = 3;
typedef itk::Image< PixelType,  Dimension >  ImageType;
typedef LesionSegmentationModel::TrainingArrayType MeasurementArrayType;
typedef LesionSegmentationModel::FLANNMatrixType FLANNMatrixType;

const PixelType imageExclusion = itk::NumericTraits<PixelType>::min( PixelType() );

template<class DetectedPixelType> ImageType::Pointer castImage(std::string& imageFileName,DetectedPixelType)
{
  typedef itk::Image<DetectedPixelType,Dimension> DetectedImageType;
  typedef itk::ImageFileReader< DetectedImageType  >  ImageReaderType;
  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(imageFileName);

  typedef itk::CastImageFilter<DetectedImageType,ImageType> CastImageFilterType;
  typename CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
  castImageFilter->SetInput(imageReader->GetOutput());
  castImageFilter->Update();

  return castImageFilter->GetOutput();
}

ImageType::Pointer ReadAndConvertImage(std::string& imageFileName)
{
  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  itk::GetImageType(imageFileName, pixelType, componentType);

  try
    {
    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        return castImage<unsigned char>(imageFileName,static_cast<unsigned char>(0));
        break;
      case itk::ImageIOBase::CHAR:
        return castImage<char>(imageFileName,static_cast<char>(0));
        break;
      case itk::ImageIOBase::USHORT:
        return castImage<unsigned short>(imageFileName,static_cast<unsigned short>(0));
        break;
      case itk::ImageIOBase::SHORT:
        return castImage<short>(imageFileName,static_cast<short>(0));
        break;
      case itk::ImageIOBase::UINT:
        return castImage<unsigned int>(imageFileName,static_cast<unsigned int>(0));
        break;
      case itk::ImageIOBase::INT:
        return castImage<int>(imageFileName,static_cast<int>(0));
        break;
      case itk::ImageIOBase::ULONG:
        return castImage<unsigned long>(imageFileName,static_cast<unsigned long>(0));
        break;
      case itk::ImageIOBase::LONG:
        return castImage<long>(imageFileName,static_cast<long>(0));
        break;
      case itk::ImageIOBase::FLOAT:
        return castImage<float>(imageFileName,static_cast<float>(0));
        break;
      case itk::ImageIOBase::DOUBLE:
        return castImage<double>(imageFileName,static_cast<double>(0));
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "  ERROR: Unknown pixel type!" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  catch( itk::ExceptionObject &excep)
    {
    std::cerr << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    exit(EXIT_FAILURE);
    }
}

ImageType::Pointer HistogramMatch(ImageType::Pointer referenceImage, ImageType::Pointer movingImage)
{
  size_t numMatchPoints = 1000;

  typedef itk::HistogramMatchingImageFilter< ImageType, ImageType> HistogramMatchingFilterType;
  HistogramMatchingFilterType::Pointer histogramMatchingFilter = HistogramMatchingFilterType::New();

  typedef HistogramMatchingFilterType::OutputImageType HistMatchType;

  typedef itk::StatisticsImageFilter< ImageType > StatisticsFilterType;
  StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
  typedef itk::CastImageFilter< HistMatchType,ImageType > CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();

  size_t numOfBins = 0;

  statisticsFilter->SetInput(referenceImage);
  statisticsFilter->Update();
  numOfBins = statisticsFilter->GetMaximum() - statisticsFilter->GetMinimum() + 1;

  histogramMatchingFilter->SetSourceImage(movingImage);
  histogramMatchingFilter->SetReferenceImage(referenceImage);
  histogramMatchingFilter->SetNumberOfHistogramLevels(numOfBins);
  histogramMatchingFilter->SetNumberOfMatchPoints(numMatchPoints);

  castFilter->SetInput(histogramMatchingFilter->GetOutput());
  castFilter->Update();
  return castFilter->GetOutput();
}

/* Adapted from nanoflann/examples/vector_of_vectors_example.cpp */
template <class FLANNMatrixType, typename num_t = double, int DIM = -1, class Distance = nanoflann::metric_L2, typename IndexType = size_t>
struct KDTreeVectorOfVectorsAdaptor
{
	typedef KDTreeVectorOfVectorsAdaptor<FLANNMatrixType,num_t,DIM,Distance> self_t;
	typedef typename Distance::template traits<num_t,self_t>::distance_t metric_t;
	typedef nanoflann::KDTreeSingleIndexAdaptor< metric_t,self_t,DIM,IndexType>  index_t;

	index_t* index; //! The kd-tree index for the user to call its methods as usual with any other FLANN index.

	/// Constructor: takes a const ref to the vector of vectors object with the data points
	KDTreeVectorOfVectorsAdaptor(const int dimensionality, const FLANNMatrixType &mat, const int leaf_max_size = 10) : m_data(mat)
	{
		assert(mat.size()!=0 && mat[0].Size()!=0);
		const size_t dims = mat[0].Size();
		if (DIM>0 && static_cast<int>(dims)!=DIM)
			throw std::runtime_error("Data set dimensionality does not match the 'DIM' template argument");
		index = new index_t( dims, *this /* adaptor */, nanoflann::KDTreeSingleIndexAdaptorParams(leaf_max_size, dims ) );
		index->buildIndex();
	}

	~KDTreeVectorOfVectorsAdaptor() {
		delete index;
	}

	const FLANNMatrixType &m_data;

	/** Query for the \a num_closest closest points to a given point (entered as query_point[0:dim-1]).
	  *  Note that this is a short-cut method for index->findNeighbors().
	  *  The user can also call index->... methods as desired.
	  * \note nChecks_IGNORED is ignored but kept for compatibility with the original FLANN interface.
	  */
	inline void QueryTree(const num_t *query_point, const size_t num_closest, IndexType *out_indices, num_t *out_distances_sq, const int nChecks_IGNORED = 10) const
	{
		//nanoflann::KNNResultSet<typename FLANNMatrixType::Scalar,IndexType> resultSet(num_closest);
		nanoflann::KNNResultSet<float ,IndexType> resultSet(num_closest);
		resultSet.init(out_indices, out_distances_sq);
		index->findNeighbors(resultSet, query_point, nanoflann::SearchParams());
	}

	/** @name Interface expected by KDTreeSingleIndexAdaptor
	  * @{ */

	const self_t & derived() const {
		return *this;
	}
	self_t & derived()       {
		return *this;
	}

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const {
		return m_data.size();
	}

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline num_t kdtree_distance(const num_t *p1, const size_t idx_p2,size_t size) const
	{
		num_t s=0;
		for (size_t i=0; i<size; i++) {
			const num_t d= p1[i]-m_data[idx_p2][i];
			s+=d*d;
		}
		return s;
	}

	// Returns the dim'th component of the idx'th point in the class:
	inline num_t kdtree_get_pt(const size_t idx, int dim) const {
		return m_data[idx][dim];
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const {
		return false;
	}

	/** @} */

}; // end of KDTreeVectorOfVectorsAdaptor

int DoIt(std::string inputT1Volume, std::string inputT2Volume, std::string inputFLAIRVolume, std::string inputMaskVolume, std::string inputT1RefVolume, std::string inputT2RefVolume, std::string inputFLAIRRefVolume, std::string inputMaskRefVolume, std::string inputModel, int inputLesionThreshold, std::string outputLesionVolume, std::string outputLesionProbVolume)
{
  typedef itk::ImageFileReader< ImageType  >  ReaderType;

  typedef itk::MaskImageFilter< ImageType, ImageType, ImageType > MaskFilterType;
  MaskFilterType::Pointer flairMaskFilter = MaskFilterType::New();
  MaskFilterType::Pointer t1MaskFilter = MaskFilterType::New();
  MaskFilterType::Pointer t2MaskFilter = MaskFilterType::New();
  MaskFilterType::Pointer clippedT1Filter = MaskFilterType::New();
  MaskFilterType::Pointer refT1MaskFilter = MaskFilterType::New();
  MaskFilterType::Pointer refT2MaskFilter = MaskFilterType::New();
  MaskFilterType::Pointer refFLAIRMaskFilter = MaskFilterType::New();

  LesionSegmentationModel lesionSegmentationModel = LesionSegmentationModel();
  lesionSegmentationModel.ReadModel(inputModel);
  const unsigned char numFeatures = lesionSegmentationModel.GetNumFeatures();
  typedef LesionSegmentationModel::LabelVectorType LabelVectorType;

  MeasurementArrayType testMeans; testMeans.Fill(0);  
  MeasurementArrayType testSigmas; testSigmas.Fill(0);  
  MeasurementArrayType trainMins = lesionSegmentationModel.GetTrainingMins();
  MeasurementArrayType trainMaxes = lesionSegmentationModel.GetTrainingMaxes();
  MeasurementArrayType trainSignedRangeInverse = lesionSegmentationModel.GetTrainingSignedRangeInverse();
  LabelVectorType trainLabels = lesionSegmentationModel.GetTrainingLabels();

  typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageRegionIteratorType; 

  typedef itk::BinaryThresholdImageFilter<ImageType,ImageType> BinaryThresholdFilterType;
  BinaryThresholdFilterType::Pointer maskBinaryFilter = BinaryThresholdFilterType::New();

  // Indexes are -1 because the input is 1 indexed, not zero indexed.
  ImageType::Pointer refT1 = ReadAndConvertImage(inputT1RefVolume);
  ImageType::Pointer refT2 = ReadAndConvertImage(inputT2RefVolume);
  ImageType::Pointer refFLAIR = ReadAndConvertImage(inputFLAIRRefVolume);
  ImageType::Pointer refMask = ReadAndConvertImage(inputMaskRefVolume);

  maskBinaryFilter->SetInput( refMask );
  maskBinaryFilter->SetLowerThreshold( 1 );
  maskBinaryFilter->SetOutsideValue( 0 );
  maskBinaryFilter->SetInsideValue( 1 );

  refT1MaskFilter->SetInput1(refT1);
  refT1MaskFilter->SetInput2(maskBinaryFilter->GetOutput());
  refT1MaskFilter->Update();

  refT2MaskFilter->SetInput1(refT2);
  refT2MaskFilter->SetInput2(maskBinaryFilter->GetOutput());
  refT2MaskFilter->Update();

  refFLAIRMaskFilter->SetInput1(refFLAIR);
  refFLAIRMaskFilter->SetInput2(maskBinaryFilter->GetOutput());
  refFLAIRMaskFilter->Update();

  typedef itk::BinaryBallStructuringElement<PixelType,Dimension > StructuringElementType;
  StructuringElementType structuringElement2;  
  structuringElement2.SetRadius( 2 ); 
  structuringElement2.CreateStructuringElement();  

  StructuringElementType structuringElement3;  
  structuringElement3.SetRadius( 3 ); 
  structuringElement3.CreateStructuringElement();  

  typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructuringElementType > DilateFilterType; 
  DilateFilterType::Pointer flairGrayscaleDilate2 = DilateFilterType::New();
  flairGrayscaleDilate2->SetKernel( structuringElement2 );

  typedef itk::FlipImageFilter <ImageType> FlipImageFilterType;
  FlipImageFilterType::Pointer flipT1Filter = FlipImageFilterType::New();
  flipT1Filter->FlipAboutOriginOff();
  FlipImageFilterType::FlipAxesArrayType flipAxes;
  flipAxes[0] = true;
  flipAxes[1] = false;
  flipAxes[2] = false;
  flipT1Filter->SetFlipAxes(flipAxes);

  typedef itk::SubtractImageFilter <ImageType,ImageType> SubtractFilterType;
  SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();

  typedef itk::ChangeInformationImageFilter < ImageType > ChangeImageFilterType;
  ChangeImageFilterType::Pointer changeImageFilter = ChangeImageFilterType::New();

  typedef itk::LabelStatisticsImageFilter< ImageType, ImageType > LabelStatisticsFilterType;
  typedef LabelStatisticsFilterType::RealType StatisticRealType;
  LabelStatisticsFilterType::Pointer statisticsFilter = LabelStatisticsFilterType::New();

  typedef itk::ScalarImageKmeansImageFilter< ImageType > KMeansFilterType;
  typedef KMeansFilterType::RealPixelType RealPixelType;
  KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();
  typedef KMeansFilterType::OutputImageType LabelImageType;

  typedef itk::CastImageFilter<LabelImageType,ImageType> CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();

  typedef itk::BinaryThresholdImageFilter<LabelImageType,LabelImageType> BinaryThresholdLabelFilterType;
  BinaryThresholdLabelFilterType::Pointer grayBinaryFilter = BinaryThresholdLabelFilterType::New();
  BinaryThresholdLabelFilterType::Pointer whiteBinaryFilter = BinaryThresholdLabelFilterType::New();
  BinaryThresholdLabelFilterType::Pointer csfBinaryFilter = BinaryThresholdLabelFilterType::New();

  typedef itk::SignedMaurerDistanceMapImageFilter< LabelImageType, LabelImageType> DistanceMapFilterType;
  DistanceMapFilterType::Pointer grayDistanceMapFilter = DistanceMapFilterType::New();
  DistanceMapFilterType::Pointer whiteDistanceMapFilter = DistanceMapFilterType::New();
  typedef DistanceMapFilterType::OutputImageType DistanceMapImageType;

  typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementType > ErodeFilterType; 
  ErodeFilterType::Pointer flairGrayscaleErode3 = ErodeFilterType::New();
  flairGrayscaleErode3->SetKernel( structuringElement3 );

  typedef itk::MedianImageFilter< ImageType, ImageType > MedianFilterType;
  MedianFilterType::Pointer flairMedian3Filter = MedianFilterType::New();

  MedianFilterType::InputSizeType radius3;
  radius3[0] = 3;
  radius3[1] = 3;
  radius3[2] = 3;

  try
    {

    ImageType::Pointer t1Image = ReadAndConvertImage(inputT1Volume);
    ImageType::Pointer t2Image = ReadAndConvertImage(inputT2Volume);
    ImageType::Pointer flairImage = ReadAndConvertImage(inputFLAIRVolume);
    ImageType::Pointer tempMaskImage = ReadAndConvertImage(inputMaskVolume);

    maskBinaryFilter->SetInput( tempMaskImage );
    maskBinaryFilter->SetLowerThreshold( 1 );
    maskBinaryFilter->SetOutsideValue( 0 );
    maskBinaryFilter->SetInsideValue( 1 );
    maskBinaryFilter->Modified();
    ImageType::Pointer maskImage = maskBinaryFilter->GetOutput();

    t1MaskFilter->SetInput1( t1Image );
    t1MaskFilter->SetInput2( maskImage );
    t1MaskFilter->SetOutsideValue(0);
    t1MaskFilter->Modified();
    t1MaskFilter->Update();

    t2MaskFilter->SetInput1( t2Image );
    t2MaskFilter->SetInput2( maskImage );
    t2MaskFilter->SetOutsideValue(0);
    t2MaskFilter->Modified();
    t2MaskFilter->Update();

    flairMaskFilter->SetInput1( flairImage );
    flairMaskFilter->SetInput2( maskImage );
    flairMaskFilter->SetOutsideValue(0);
    flairMaskFilter->Modified();
    flairMaskFilter->Update();

    ImageType::Pointer t1ImageHistMatched = HistogramMatch(refT1MaskFilter->GetOutput(),t1MaskFilter->GetOutput());
    ImageType::Pointer t2ImageHistMatched = HistogramMatch(refT2MaskFilter->GetOutput(),t2MaskFilter->GetOutput());
    ImageType::Pointer flairImageHistMatched = HistogramMatch(refFLAIRMaskFilter->GetOutput(),flairMaskFilter->GetOutput());

    /* Iterate over the flair, get it's intensity and it's location.  Use those
     ** values to calculate the mean and standard deviation.  Iterate over the 
     ** flair again, zero-mean and sigma correct, and push the measurement vector
     ** onto a listsample.
     */

    flairGrayscaleDilate2->SetInput( flairImageHistMatched );
    flairGrayscaleDilate2->Modified();
    flairGrayscaleDilate2->Update();

    flipT1Filter->SetInput( t1ImageHistMatched );
    flipT1Filter->Modified();
    flipT1Filter->Update();

    changeImageFilter->SetInput( flipT1Filter->GetOutput());
    changeImageFilter->SetOutputOrigin( t1ImageHistMatched->GetOrigin() );
    changeImageFilter->SetOutputDirection( t1ImageHistMatched->GetDirection() );
    changeImageFilter->SetChangeOrigin( true );
    changeImageFilter->SetChangeDirection( true );
    changeImageFilter->Modified();
    changeImageFilter->Update();

    subtractFilter->SetInput1( t1ImageHistMatched );
    subtractFilter->SetInput2( changeImageFilter->GetOutput() );
    subtractFilter->Modified();
    subtractFilter->Update();

    clippedT1Filter->SetInput1( t1ImageHistMatched );
    clippedT1Filter->SetInput2( maskImage );
    clippedT1Filter->SetOutsideValue( imageExclusion );
    clippedT1Filter->Modified();
    clippedT1Filter->Update();

    statisticsFilter->SetInput( t1ImageHistMatched );
    statisticsFilter->SetLabelInput( maskImage );
    statisticsFilter->Modified();
    statisticsFilter->Update();

    const StatisticRealType imageMean = statisticsFilter->GetMean( 1 );
    const StatisticRealType imageSigma = statisticsFilter->GetSigma( 1 );

    const RealPixelType csfInitialMean = imageMean - 2*imageSigma;
    const RealPixelType whiteInitialMean = imageMean + imageSigma;
    const RealPixelType grayInitialMean = imageMean - imageSigma/5;

    RealPixelType backgroundInitialMean = imageExclusion;
    kmeansFilter->AddClassWithInitialMean( backgroundInitialMean );
    kmeansFilter->AddClassWithInitialMean( csfInitialMean );
    kmeansFilter->AddClassWithInitialMean( grayInitialMean );
    kmeansFilter->AddClassWithInitialMean( whiteInitialMean );
    kmeansFilter->SetInput( clippedT1Filter->GetOutput() );
    kmeansFilter->Modified();
    castFilter->SetInput(kmeansFilter->GetOutput());
    castFilter->Modified();
    castFilter->Update();

    grayBinaryFilter->SetInput(kmeansFilter->GetOutput());
    grayBinaryFilter->SetLowerThreshold(2);
    grayBinaryFilter->SetUpperThreshold(2);
    grayBinaryFilter->Modified();
    grayDistanceMapFilter->SetInput( grayBinaryFilter->GetOutput() );
    grayDistanceMapFilter->UseImageSpacingOff();
    grayDistanceMapFilter->SquaredDistanceOff();
    grayDistanceMapFilter->Modified();
    grayDistanceMapFilter->Update();

    whiteBinaryFilter->SetInput(kmeansFilter->GetOutput());
    whiteBinaryFilter->SetLowerThreshold(1);
    whiteBinaryFilter->SetUpperThreshold(1);
    whiteBinaryFilter->Modified();
    whiteDistanceMapFilter->SetInput( whiteBinaryFilter->GetOutput() );
    whiteDistanceMapFilter->UseImageSpacingOff();
    whiteDistanceMapFilter->SquaredDistanceOff();
    whiteDistanceMapFilter->Modified();
    whiteDistanceMapFilter->Update();

    flairGrayscaleErode3->SetInput( flairImageHistMatched );
    flairGrayscaleErode3->Modified();
    flairGrayscaleErode3->Update();

    flairMedian3Filter->SetInput( flairImageHistMatched );
    flairMedian3Filter->SetRadius( radius3 );
    flairMedian3Filter->Modified();
    flairMedian3Filter->Update();

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    /* The means and stds should be calculated over the lesion and non-lesion voxels
     ** at the same time, and then used seperately to get the lesion and non-lesion
     ** samples
     */
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    /* Iterator over the lesion masked, brain masked, flair */
    ImageRegionIteratorType flairItr( flairImageHistMatched,flairImageHistMatched->GetRequestedRegion() ); 

    size_t count = 0;

    /* First pass is to calculate the mean and standard deviation of the x, y, z
     ** voxel locations and the mean and std of the flair intensity.
     */
    MeasurementArrayType tempSamples; tempSamples.Fill(0);

    for ( flairItr.GoToBegin(); !flairItr.IsAtEnd(); ++flairItr) 
      { 
      ImageType::IndexType idx = flairItr.GetIndex();

      if(maskImage->GetPixel(idx) != 0)
        {
        tempSamples[0] = idx[0];
        tempSamples[1] = idx[1];
        tempSamples[2] = idx[2];
        tempSamples[3] = flairGrayscaleDilate2->GetOutput()->GetPixel(idx);
        tempSamples[4] = subtractFilter->GetOutput()->GetPixel(idx);
        tempSamples[5] = grayDistanceMapFilter->GetOutput()->GetPixel(idx);
        tempSamples[6] = whiteDistanceMapFilter->GetOutput()->GetPixel(idx);
        tempSamples[7] = flairGrayscaleErode3->GetOutput()->GetPixel(idx);
        tempSamples[8] = t2ImageHistMatched->GetPixel(idx);
        tempSamples[9] = flairMedian3Filter->GetOutput()->GetPixel(idx);

        for( size_t i=0; i<tempSamples.Size(); i++ )
          {
          testMeans[i] += tempSamples[i];
          testSigmas[i] += tempSamples[i] * tempSamples[i];
          }

        ++count;
        }
      }


    const size_t numNeighbors = 30;
    typedef KDTreeVectorOfVectorsAdaptor< FLANNMatrixType, float >  my_kd_tree_t;
    FLANNMatrixType trainingDataset = lesionSegmentationModel.GetFLANNDataset();
    std::cout << trainingDataset.size() << ":" << trainingDataset[0].Size() << "\n";

    my_kd_tree_t   trainingIndex(numFeatures, trainingDataset, 10 /* max leaf */ );

    std::vector<size_t> indices(numNeighbors,0);
    std::vector<float> distsSquared(numNeighbors,0);
    std::vector<float> query(numFeatures,0);

    nanoflann::KNNResultSet<float> resultSet(numNeighbors);

    for( size_t w=0; w<testMeans.Size(); w++ )
      {
      testMeans[w] = testMeans[w] / count;
      testSigmas[w] = testSigmas[w] / count;
      testSigmas[w] = sqrt(testSigmas[w] - (testMeans[w]*testMeans[w]));
      }

    typedef itk::Image<short,Dimension> OutputImageType;
    OutputImageType::Pointer lesionMask = OutputImageType::New();
    OutputImageType::RegionType region;
    region.SetSize(t1Image->GetLargestPossibleRegion().GetSize());
    region.SetIndex(t1Image->GetLargestPossibleRegion().GetIndex());
    lesionMask->SetRegions(region);
    lesionMask->SetSpacing(t1Image->GetSpacing());
    lesionMask->SetOrigin(t1Image->GetOrigin());
    lesionMask->SetDirection(t1Image->GetDirection());
    lesionMask->Allocate();

    OutputImageType::Pointer percentLesionImage = OutputImageType::New();
    percentLesionImage->SetRegions(region);
    percentLesionImage->SetSpacing(t1Image->GetSpacing());
    percentLesionImage->SetOrigin(t1Image->GetOrigin());
    percentLesionImage->SetDirection(t1Image->GetDirection());
    percentLesionImage->Allocate();


    /* Grab the index and value, zero-mean and sigma correct, create a measurement vector,
     ** and push it onto the list of samples while pushing the label onto the list of labels.
     */
    size_t lesionCount = 0;
    for ( flairItr.GoToBegin(); !flairItr.IsAtEnd(); ++flairItr)
      {

      ImageType::IndexType idx = flairItr.GetIndex();
      if(maskImage->GetPixel(idx) != 0)
        {
        query[0] = (idx[0]-testMeans[0])/testSigmas[0];
        query[1] = (idx[1]-testMeans[1])/testSigmas[1];
        query[2] = (idx[2]-testMeans[2])/testSigmas[2];
        query[3] = (flairGrayscaleDilate2->GetOutput()->GetPixel(idx)-testMeans[3]-1)/testSigmas[3];
        query[4] = (subtractFilter->GetOutput()->GetPixel(idx)-testMeans[4]-1)/testSigmas[4];
        query[5] = (grayDistanceMapFilter->GetOutput()->GetPixel(idx)-testMeans[5]-1)/testSigmas[5];
        query[6] = (whiteDistanceMapFilter->GetOutput()->GetPixel(idx)-testMeans[6]-1)/testSigmas[6];
        query[7] = (flairGrayscaleErode3->GetOutput()->GetPixel(idx)-testMeans[7]-1)/testSigmas[7];
        query[8] = (t2ImageHistMatched->GetPixel(idx) -testMeans[8]-1)/testSigmas[8];
        query[9] = (flairMedian3Filter->GetOutput()->GetPixel(idx)-testMeans[9]-1)/testSigmas[9];

        for(size_t j=0;j<numFeatures;j++)
          {
          query[j] = trainSignedRangeInverse[j] * (query[j] - trainMins[j]) - 1;
          }

        trainingIndex.QueryTree(&query[0], numNeighbors, &indices[0], &distsSquared[0]);

        size_t numLesion = 0;
        for(size_t j=0;j<numNeighbors;j++)
          {
          if(trainLabels[indices[j]] > 0)
            {//lesion
            numLesion++;
            }
          }
        OutputImageType::PixelType chanceLesion = (OutputImageType::PixelType)(100*((float)(numLesion)/(float)(numNeighbors)));

        percentLesionImage->SetPixel(idx,chanceLesion);

        if(chanceLesion > inputLesionThreshold)
          {
          lesionMask->SetPixel(idx,1);
          lesionCount++;
          }
        else
          {
          lesionMask->SetPixel(idx,0);
          }
        }
      else
        {
        lesionMask->SetPixel(idx,0);
        percentLesionImage->SetPixel(idx,0);
        }
      }

    std::cout << "lesionCount=" << lesionCount << "\n";

    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    WriterType::Pointer lesionMaskWriter = WriterType::New();
    lesionMaskWriter->SetInput(lesionMask);
    lesionMaskWriter->SetUseCompression(true);
    lesionMaskWriter->SetFileName(outputLesionVolume);
    lesionMaskWriter->Update();

    WriterType::Pointer lesionPercentWriter = WriterType::New();
    lesionPercentWriter->SetInput(percentLesionImage);
    lesionPercentWriter->SetUseCompression(true);
    lesionPercentWriter->SetFileName(outputLesionProbVolume);
    lesionPercentWriter->Update();


    }
  catch (itk::ExceptionObject &excep)
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
} // End anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  bool violated=false;
  if (inputT1Volume.size() == 0) { violated = true; std::cout << "  --inputT1Volume Required! "  << std::endl; }
  if (inputT2Volume.size() == 0) { violated = true; std::cout << "  --inputT2Volume Required! "  << std::endl; }
  if (inputFLAIRVolume.size() == 0) { violated = true; std::cout << "  --inputFLAIRVolume Required! "  << std::endl; }
  if (inputMaskVolume.size() == 0) { violated = true; std::cout << "  --inputMaskVolume Required! "  << std::endl; }
  if (inputT1RefVolume.size() == 0) { violated = true; std::cout << "  --inputT1RefVolume Required! "  << std::endl; }
  if (inputT2RefVolume.size() == 0) { violated = true; std::cout << "  --inputT2RefVolume Required! "  << std::endl; }
  if (inputFLAIRRefVolume.size() == 0) { violated = true; std::cout << "  --inputFLAIRRefVolume Required! "  << std::endl; }
  if (inputMaskRefVolume.size() == 0) { violated = true; std::cout << "  --inputMaskRefVolume Required! "  << std::endl; }
  if (inputModel.size() == 0) { violated = true; std::cout << "  --inputModel Required! "  << std::endl; }
  if (outputLesionVolume.size() == 0) { violated = true; std::cout << "  --outputLesionVolume Required! "  << std::endl; }
  if (outputLesionProbVolume.size() == 0) { violated = true; std::cout << "  --outputLesionProbVolume Required! "  << std::endl; }

  if (violated) exit(EXIT_FAILURE);

  try
    {
    return DoIt(inputT1Volume, inputT2Volume, inputFLAIRVolume, inputMaskVolume, inputT1RefVolume, inputT2RefVolume, inputFLAIRRefVolume, inputMaskRefVolume, inputModel, inputLesionThreshold, outputLesionVolume, outputLesionProbVolume);
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
}
