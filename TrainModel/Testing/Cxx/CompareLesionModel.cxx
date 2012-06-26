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

#include "algorithm"
#include "LesionSegmentationModel.h"
#include "CompareLesionModelCLP.h"

namespace{

bool DoIt(std::string inputModel1, std::string inputModel2, int inputPercentDiff)
  {
  LesionSegmentationModel lesionModel1 = LesionSegmentationModel();
  LesionSegmentationModel lesionModel2 = LesionSegmentationModel();
  lesionModel1.ReadModel(inputModel1);
  lesionModel2.ReadModel(inputModel2);
  if(lesionModel1 == lesionModel2)
    {
    return EXIT_SUCCESS;
    }
  else // Models aren't exactly equal.  Are they equal within a threshold?
    {

    if(lesionModel1.GetNumFeatures() != lesionModel2.GetNumFeatures())
      {
      std::cerr << "Number of features do not match!" << std::endl;
      return EXIT_FAILURE;
      }
    if(lesionModel1.GetTrainingMins() != lesionModel2.GetTrainingMins())
      {
      std::cerr << "TrainingMins do not match!" << std::endl;
      return EXIT_FAILURE;
      }
    if(lesionModel1.GetTrainingMaxes() != lesionModel2.GetTrainingMaxes())
      {
      std::cerr << "TrainingMaxes do not match!" << std::endl;
      return EXIT_FAILURE;
      }

    
    int labelFraction = 100 * (int)abs(lesionModel1.GetTrainingLabels().size() - lesionModel2.GetTrainingLabels().size())/((lesionModel1.GetTrainingLabels().size() + lesionModel2.GetTrainingLabels().size())/2);
    if(labelFraction > inputPercentDiff)
      {
      std::cerr << "Label count not within " << inputPercentDiff << "% of each other!" << std::endl;
      std::cerr << "model1LabelCount: " << lesionModel1.GetTrainingLabels().size();
      std::cerr << " model2LesionCount: " << lesionModel2.GetTrainingLabels().size();
      return EXIT_FAILURE;
      }
    size_t model1LesionCount = (size_t)std::count(lesionModel1.GetTrainingLabels().begin(),lesionModel1.GetTrainingLabels().end(),1);
    size_t model2LesionCount = (size_t)std::count(lesionModel2.GetTrainingLabels().begin(),lesionModel2.GetTrainingLabels().end(),1);
    if(model1LesionCount != model2LesionCount)
      {
      std::cerr << "Lesion count not within " << inputPercentDiff << "% of each other!" << std::endl;
      std::cerr << "model1LesionCount: " << model1LesionCount << " model2LesionCount: " << model2LesionCount << std::endl;
      return EXIT_FAILURE;
      }

    size_t model1NonLesionCount = (size_t)std::count(lesionModel1.GetTrainingLabels().begin(),lesionModel1.GetTrainingLabels().end(),0);
    size_t model2NonLesionCount = (size_t)std::count(lesionModel2.GetTrainingLabels().begin(),lesionModel2.GetTrainingLabels().end(),0);
    int nonLesionFraction = 100* (int)abs(model1NonLesionCount - model2NonLesionCount) / ((model1NonLesionCount+model2NonLesionCount)/2);
    if(nonLesionFraction > inputPercentDiff*2)
      {
      std::cerr << "Nonlesion count not within " << inputPercentDiff*2 << "% of each other!" << std::endl;
      std::cerr << "model1NonLesionCount: " << model1NonLesionCount << " model2NonLesionCount: " << model2NonLesionCount << std::endl;
      return EXIT_FAILURE;
      }

    return EXIT_SUCCESS;
    }
  }
}
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  bool violated=false;
  if (inputModel1.size() == 0) { violated = true; std::cout << "  --inputModel1 Required! "  << std::endl; }
  if (inputModel2.size() == 0) { violated = true; std::cout << "  --inputModel2 Required! "  << std::endl; }
  if (violated) exit(EXIT_FAILURE);

  return DoIt(inputModel1, inputModel2, inputPercentDiff);
}
 
