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
#include "TrainModelCLP.h"

/* Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
*/ 
namespace
{

typedef unsigned short      PixelType;
const unsigned int          Dimension = 3;
typedef itk::Image< PixelType,  Dimension >  ImageType;

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
  typename ImageType::Pointer castImage = castImageFilter->GetOutput();

  return castImage;
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

int DoIt(int argc, char * argv [])
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
  if (violated) exit(EXIT_FAILURE);

  typedef itk::Image< PixelType,  Dimension >   InputImageType;
  typedef itk::Image< float, Dimension > OutputImageType;
  typedef itk::ImageFileReader< InputImageType  >  ReaderType;

  ReaderType::Pointer lesionArrayReader = ReaderType::New();
  ReaderType::Pointer maskArrayReader = ReaderType::New();

  /* TODO Need to read in the images indicated by, inputIndexOfBestImages, build a joint histogram, 
   * then stash the histogram so it can be used to intensity standardize the images as they are loaded.
   * When loading these images again the intensity standardization step can be skipped.
   * */

  typedef itk::MaskImageFilter< InputImageType, InputImageType, InputImageType > MaskFilterType;
  MaskFilterType::Pointer brainMaskFilter = MaskFilterType::New();
  MaskFilterType::Pointer lesionMaskFilter = MaskFilterType::New();
  MaskFilterType::Pointer nonLesionMaskFilter = MaskFilterType::New();

  typedef itk::InvertIntensityImageFilter< InputImageType, InputImageType > InvertFilterType;
  InvertFilterType::Pointer invertFilter = InvertFilterType::New();

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
      InputImageType::Pointer t1Image = ReadAndConvertImage(inputT1Volumes[i]);

      InputImageType::Pointer t2Image = ReadAndConvertImage(inputT2Volumes[i]);

      InputImageType::Pointer flairImage = ReadAndConvertImage(inputFLAIRVolumes[i]);

      maskArrayReader->SetFileName( inputMaskVolumes[i].c_str() );
      maskArrayReader->Modified();

      lesionArrayReader->SetFileName( inputLesionVolumes[i].c_str() );
      lesionArrayReader->Modified();

      /* Brain mask the ith trianing flair */
      brainMaskFilter->SetInput1( flairImage );
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
        InputImageType::IndexType idx = flairItr.GetIndex();

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
        InputImageType::IndexType idx = flairItr.GetIndex();

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
  if (violated) exit(EXIT_FAILURE);

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;
  itk::ImageIOBase::IOPixelType     firstPixelType;
  itk::ImageIOBase::IOComponentType firstComponentType;

  try
    {

    itk::GetImageType(inputFLAIRVolumes[0].c_str(), firstPixelType, firstComponentType);

    for (unsigned int i=0; i<inputFLAIRVolumes.size(); i++)
      {
      itk::GetImageType(inputFLAIRVolumes[i].c_str(),  pixelType, componentType);
      if(firstComponentType != componentType)
        {
        std::cout << "  The pixel types of the input FLAIRs are not the same." << std::endl;
        exit(EXIT_FAILURE);
        }
      itk::GetImageType(inputT1Volumes[i].c_str(),  pixelType, componentType);
      if(firstComponentType != componentType)
        {
        std::cout << "  The pixel types of the input T1s are not the same, or do not match the pixel type of the FLAIRS." << std::endl;
        exit(EXIT_FAILURE);
        }
      itk::GetImageType(inputT1Volumes[i].c_str(),  pixelType, componentType);
      if(firstComponentType != componentType)
        {
        std::cout << "  The pixel types of the input T2s are not the same, or do not match the pixel type of the FLAIRS." << std::endl;
        exit(EXIT_FAILURE);
        }
      }

    return DoIt( argc, argv );
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
