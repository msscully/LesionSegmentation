#include <iostream>
#include <vector>
#include <utility>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkListSample.h"
#include "itkMembershipSample.h"
#include "itkVector.h"
#include "itkPluginFilterWatcher.h"
#include "itkImageToVectorImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkKdTreeGenerator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkAddImageAdaptor.h"
#include "TrainModelCLP.h"

/* Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
*/ 
namespace
{

template <class T>
int DoIt(int argc, char * argv [], T)
{
  PARSE_ARGS;

  bool violated=false;
  if (inputFLAIRVolumes.size() == 0) { violated = true; std::cout << "  --inputFLAIRVolumes Required! "  << std::endl; }
  if (inputLesionVolumes.size() == 0) { violated = true; std::cout << "  --inputLesionVolumes Required! "  << std::endl; }
  if (inputMaskVolumes.size() == 0) { violated = true; std::cout << "  --inputMaskVolumes Required! "  << std::endl; }
  if (inputT1Volumes.size() == 0) { violated = true; std::cout << "  --inputT1Volumes Required! "  << std::endl; }
  if (inputT2Volumes.size() == 0) { violated = true; std::cout << "  --inputT2Volumes Required! "  << std::endl; }
  if ((inputFLAIRVolumes.size() != inputLesionVolumes.size()) && 
    (inputMaskVolumes.size() != inputFLAIRVolumes.size()) &&
    (inputT1Volumes.size() != inputFLAIRVolumes.size()) &&
    (inputT2Volumes.size() != inputFLAIRVolumes.size())) 
    { violated = true; std::cout << "  the number of files after --inputT1Volumes, --inputT2Volumes, --inputFLAIRVolumes, --inputLesionVolumes, and --inputMaskVolumes must all be equal! "  << std::endl;
    }
  if (violated) exit(1);

  typedef T PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image< PixelType,  Dimension >   InputImageType;
  typedef itk::Image< float, Dimension > OutputImageType;
  typedef itk::ImageFileReader< InputImageType  >  ReaderType;

  typename ReaderType::Pointer lesionArrayReader = ReaderType::New();
  typename ReaderType::Pointer maskArrayReader = ReaderType::New();
  typename ReaderType::Pointer flairArrayReader = ReaderType::New();
  typename ReaderType::Pointer t1ArrayReader = ReaderType::New();
  typename ReaderType::Pointer t2ArrayReader = ReaderType::New();

  typename InputImageType::Pointer LoadImage( std::string );
   
  /* TODO Need to read in the images indicated by, inputIndexOfBestImages, build a joint histogram, 
   * then stash the histogram so it can be used to intensity standardize the images as they are loaded.
   * When loading these images again the intensity standardization step can be skipped.
   * */

  typedef itk::MaskImageFilter< InputImageType, InputImageType, InputImageType > MaskFilterType;
  typename MaskFilterType::Pointer brainMaskFilter = MaskFilterType::New();
  typename MaskFilterType::Pointer lesionMaskFilter = MaskFilterType::New();
  typename MaskFilterType::Pointer nonLesionMaskFilter = MaskFilterType::New();

  typedef itk::InvertIntensityImageFilter< InputImageType, InputImageType > InvertFilterType;
  typename InvertFilterType::Pointer invertFilter = InvertFilterType::New();

  typedef itk::Vector< float, 6 > MeasurementVectorType ;
  typedef itk::Vector< unsigned short, 1 > LabelVectorType;

  MeasurementVectorType trainMeans;
  MeasurementVectorType trainSigmas;

  typedef itk::Statistics::ListSample< MeasurementVectorType > ListSampleType;
  ListSampleType::Pointer lesionSamples = ListSampleType::New();
  ListSampleType::Pointer nonLesionSamples = ListSampleType::New();

  typedef itk::ImageRegionIteratorWithIndex< InputImageType > ImageRegionIteratorType; 

  srand((unsigned)time(0)); 

  if (inputPercentNonLesion == 0)
    {
    inputPercentNonLesion = floor(100/(inputFLAIRVolumes.size()*2));
    }

  /* Need to loop over all the input images, brain mask them, combine the flair and lesion
   ** values into a vector image, calculate stats, normalize the data, and construct a 
   ** sample list.
   */

  try
    {

    for (unsigned int i=0; i<inputFLAIRVolumes.size(); i++)
      {
      /* Load the ith training images */
      t1ArrayReader->SetFileName( inputT1Volumes[i].c_str() );
      t1ArrayReader->Modified();

      t2ArrayReader->SetFileName( inputT2Volumes[i].c_str() );
      t2ArrayReader->Modified();

      flairArrayReader->SetFileName( inputFLAIRVolumes[i].c_str() );
      flairArrayReader->Modified();

      maskArrayReader->SetFileName( inputMaskVolumes[i].c_str() );
      maskArrayReader->Modified();

      lesionArrayReader->SetFileName( inputLesionVolumes[i].c_str() );
      lesionArrayReader->Modified();

      /* Brain mask the ith trianing flair */
      brainMaskFilter->SetInput1( flairArrayReader->GetOutput() );
      brainMaskFilter->SetInput2( maskArrayReader->GetOutput() );
      brainMaskFilter->SetOutsideValue( 0 );
      brainMaskFilter->Update();

      /* Iterate over the flair, get it's intensity and it's location.  Use those
       ** values to calculate the mean and standard deviation.  Iterate over the 
       ** flair again, zero-mean and sigma correct, and push the measurement vector
       ** onto a listsample.
       */

      /* Use the ith training lesion map to mask the brain masked flair */
      lesionMaskFilter->SetInput1( brainMaskFilter->GetOutput() );
      //lesionMaskFilter->SetInput2( LoadImage( inputLesionVolumes[i].c_str() ) );
      lesionMaskFilter->SetInput2( lesionArrayReader->GetOutput() );
      lesionMaskFilter->SetOutsideValue( 0 );

      /* Invert the ith training lesion mask to get all non lesions */
      invertFilter->SetInput( lesionArrayReader->GetOutput() );

      /* Top 10 most informative features
         Dilate FLAIR radius=2
         T1w flipped difference
         Normalized Z location
         Normalized X location
         Distance to white matter
         Normalized Y location
         Distance to gray matter
         Erode Flair radius=3
         Normalized T2w intensity
         Median FLAIR radius=3
        */

      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      /* The means and stds should be calculated over the lesion and non-lesion voxels
       ** at the same time, and then used seperately to get the lesion and non-lesion
       ** samples
       */
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      /* Iterator over the lesion masked, brain masked, flair */
      ImageRegionIteratorType flairItr( brainMaskFilter->GetOutput(),brainMaskFilter->GetOutput()->GetRequestedRegion() ); 

      unsigned int count = 0;

      /* First pass is to calculate the mean and standard deviation of the x, y, z
       ** voxel locations and the mean and std of the flair intensity.
       */
      for ( flairItr.GoToBegin(); !flairItr.IsAtEnd(); ++flairItr) 
        { 
        typename InputImageType::IndexType idx = flairItr.GetIndex();

        trainMeans[0] += idx[0];
        trainMeans[1] += idx[1];
        trainMeans[2] += idx[2];
        trainMeans[3] += flairItr.Get();

        trainSigmas[0] = idx[0] * idx[0];
        trainSigmas[1] = idx[1] * idx[1];
        trainSigmas[2] = idx[2] * idx[2];
        trainSigmas[3] = flairItr.Get() * flairItr.Get();
        ++count;
        }

      trainMeans[0] /= count;
      trainMeans[1] /= count;
      trainMeans[2] /= count;
      trainMeans[3] /= count;

      trainSigmas[0] = sqrt((trainSigmas[0] / count) - (trainMeans[0]*trainMeans[0]));
      trainSigmas[1] = sqrt((trainSigmas[1] / count) - (trainMeans[1]*trainMeans[1]));
      trainSigmas[2] = sqrt((trainSigmas[2] / count) - (trainMeans[2]*trainMeans[2]));
      trainSigmas[3] = sqrt((trainSigmas[3] / count) - (trainMeans[3]*trainMeans[3]));


      /* Grab the index and value, zero-mean and sigma correct, create a measurement vector,
       ** and push it onto the list of samples while pushing the label onto the list of labels.
       */
      for ( flairItr.GoToBegin(); !flairItr.IsAtEnd(); ++flairItr)
        {
        typename InputImageType::IndexType idx = flairItr.GetIndex();

        MeasurementVectorType tempMeasurement;
        tempMeasurement.SetElement( 0, (idx[0]-trainMeans[0])/trainSigmas[0] ); 
        tempMeasurement.SetElement( 1, (idx[1]-trainMeans[1])/trainSigmas[1] ); 
        tempMeasurement.SetElement( 2, (idx[2]-trainMeans[2])/trainSigmas[2] ); 
        tempMeasurement.SetElement( 3, (flairItr.Get()-trainMeans[3]-1)/trainSigmas[3] ); 
        
        if(maskArrayReader->GetOutput()->GetPixel(idx) !=0)          
          {
          lesionSamples->PushBack(tempMeasurement);
          }
        else if ((rand()%100+1) <= inputPercentNonLesion) /* We only want a fraction of non-lesion voxels */
          {
          nonLesionSamples->PushBack(tempMeasurement);
          }
        } 
      }



    }
  catch (itk::ExceptionObject &excep)
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
} // End anonymous namespace

int main( int argc, char * argv[] )
{
    PARSE_ARGS;

    try
      {
      return DoIt( argc, argv, static_cast<unsigned char>(0) );
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << argv[0] << ": exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
      }
    return EXIT_SUCCESS;
}
