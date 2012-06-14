#ifndef LesionSegmentationModel_h
#define LesionSegmentationModel_h

#include <vector>
//#include <vnl_vector>
#include <fstream>
#include <iostream>
#include "itkArray.h"
#include "itkVariableSizeMatrix.h"
#include "itkByteSwapper.h"
#include "itkVariableLengthVector.h"
#include <math.h>

class LesionSegmentationModel
{
private:
  typedef LesionSegmentationModel Self;
  typedef enum { readFail,writeFail } ioErr;
  typedef enum { file_signature=0x12345678,
    swapped_file_signature=0x78563412 } fileSig;

public:
  static const unsigned char m_NumFeatures = 10;
  typedef itk::FixedArray< float, m_NumFeatures > TrainingArrayType;
  typedef std::vector< char > LabelVectorType;
  //typedef flann::Index< flann::L2<float> > FLANNIndexType;
  //typedef flann::Matrix< float > FlannMatrixType;

  LesionSegmentationModel() : m_Swapped(false)
    {
    InitializeModel();
    }

  /*TODO: Get FLANN saving to work.
  FLANNIndexType & GetFLANNIndex() {return m_flannIndex;}
  void SetFLANNIndex(FLANNIndexType &inputFLANNIndex)
    {
    m_flannIndex = inputFLANNIndex;
    }

  FLANNIndexType & GetFLANNDataset() {return m_flannDataset;}
  void SetFLANNDataset(FLANNMatrixType &inputFLANNDataset)
    {
    m_flannDataset = inputFLANNDataset;
    }
  */
  TrainingArrayType & GetTrainingMins() {return m_trainingMins;}
  void SetTrainingMins(TrainingArrayType &inputTrainingMins)
    {
    m_trainingMins = inputTrainingMins;
    }

  TrainingArrayType & GetTrainingMaxes() {return m_trainingMaxes;}
  void SetTrainingMaxes(TrainingArrayType &inputTrainingMaxes)
    {
    m_trainingMaxes = inputTrainingMaxes;
    }

  TrainingArrayType & GetTrainingSignedRangeInverse() {return m_trainingSignedRangeInverse;}
  void SetTrainingSignedRangeInverse(TrainingArrayType &inputTrainingSignedRangeInverse)
    {
    m_trainingSignedRangeInverse = inputTrainingSignedRangeInverse;
    }

  LabelVectorType & GetTrainingLabels() {return m_trainingLabels;}
  void SetTrainingLabels(LabelVectorType &inputTrainingLabels)
    {
    m_trainingLabels = inputTrainingLabels;
    }

  size_t & GetLabelsSize() {return m_trainingLabelsSize;}
  void SetLabelsSize(size_t inputLabelsSize)
    {
    m_trainingLabelsSize = inputLabelsSize;
    }

  int GetNumFeatures(){return m_NumFeatures;}

  void SaveModel(const std::string &filename)
    {
    ////////////////////////////////////////////////////////////////////////////

    std::ofstream output(filename.c_str());
    if (!output.is_open())
      {
      std::cerr << "Can't write " << filename << std::endl;
      std::cerr.flush();
      exit(-1);
      }
    try
      {
      this->Write<unsigned int>(output,file_signature); //Write out the signature first
 
      this->Write(output,this->GetTrainingMins());
      this->Write(output,this->GetTrainingMaxes());
      this->Write(output,this->GetTrainingSignedRangeInverse());
      this->Write(output,this->GetLabelsSize());
      this->Write(output,this->GetTrainingLabels());
      //TODO: Write out Histograms or images for intensity standardization? 
      //this->Write(input,this->GetT1RefImage());
      //this->Write(input,this->GetT2RefImage());
      //this->Write(input,this->GetFLAIRRefImage());
      //this->Write(output,this->GetFLANNDataset());
      //this->Write(output,this->GetFLANNIndex());
     
      }
    catch (ioErr e)
      {
      std::cerr << "Write failed for " << filename << std::endl;
      }
    output.close();
    }

  void ReadModel(const std::string &filename)
    {
    ////////////////////////////////////////////////////////////////////////////

    std::ifstream input(filename.c_str()); 
    if (!input.is_open())
      {
      std::cerr << "Can't read " << filename << std::endl;
      std::cerr.flush();
      exit(-1);
      }
    try
      {
      unsigned int sig;
      this->Read<unsigned int>(input,sig);
      if(sig != file_signature && sig != swapped_file_signature)
        {
        this->m_Swapped = false;
        }
      else if(sig == swapped_file_signature)
        {
        this->m_Swapped = true;
        }
      this->Read(input,this->GetTrainingMins());
      this->Read(input,this->GetTrainingMaxes());
      this->Read(input,this->GetTrainingSignedRangeInverse());
      this->Read(input,this->GetLabelsSize());
      // This can be variable length, so we have to get the size from the 
      // file and use it to allocate memory for the read step.
      this->m_trainingLabels = LabelVectorType(this->GetLabelsSize());
      this->Read(input,this->GetTrainingLabels());
      //TODO: Read Histograms or images for intensity standardization? 
      //this->Read(input,this->GetT1RefImage());
      //this->Read(input,this->GetT2RefImage());
      //this->Read(input,this->GetFLAIRRefImage());
      //this->Read(input,this->GetFLANNDataset());
      //this->Read(input,this->GetFLANNIndex());
     
      }
    catch (ioErr e)
      {
      std::cerr << "Read failed for " << filename << std::endl;
      std::cerr << e << std::endl;
      exit(-1);
      }
    input.close();
    }
  
  void InitializeModel()
    {
    }

private:
  bool m_Swapped;

  template <class T>
    void
    Write(std::ofstream &f,T var)
      {
      if(f.bad() || f.eof())
        {
        throw writeFail;
        }
      f.write(reinterpret_cast<char *>(&var),sizeof(T));
      }
  template <class T>
    void
    Read(std::ifstream &f, T &var)
      {
      if(f.bad() || f.eof())
        {
        throw readFail;
        }
      f.read(reinterpret_cast<char *>(&var),sizeof(T));
      if(this->m_Swapped)
        {
        if(itk::ByteSwapper<T>::SystemIsBigEndian())
          {
          itk::ByteSwapper<T>::SwapFromSystemToLittleEndian(&var);
          }
        else
          {
          itk::ByteSwapper<T>::SwapFromSystemToBigEndian(&var);
          }
        }
      }

  void Write(std::ofstream &f,const TrainingArrayType &vec)
    {
    for(unsigned int y=0;y<vec.Size();y++)
      {
      this->Write<TrainingArrayType::ValueType>(f,vec[y]);
      }
    }
  void Read(std::ifstream &f,TrainingArrayType &vec)
    {
    for(unsigned int y=0;y<vec.Size();y++)
      {
      this->Read<TrainingArrayType::ValueType>(f,vec[y]);
      }
    }

  void Write(std::ofstream &f,const LabelVectorType &vec)
    {
    for(unsigned int y=0;y<vec.size();y++)
      {
      this->Write<LabelVectorType::value_type>(f,vec[y]);
      }
    }
  void Read(std::ifstream &f,LabelVectorType &vec)
    {
    for(unsigned int y=0;y<vec.size();y++)
      {
      this->Read<LabelVectorType::value_type>(f,vec[y]);
      }
    }

  /*
  void Write(std::ofstream &f,const FLANNIndexType &flannIndex)
    {
    flannIndex.saveIndex(f);
    }
  void Read(std::ifstream &f,FLANNIndexType &flannIndex)
    {
    flannIndex.loadIndex(f);
    }
  void Write(std::ofstream &f,const FLANNIndexType &flannIndex)
    {
    flannIndex.saveIndex(f);
    }
  void Read(std::ifstream &f,FLANNIndexType &flannIndex)
    {
    flannIndex.loadIndex(f);
    }
  */

  TrainingArrayType m_trainingMins;
  TrainingArrayType m_trainingMaxes;
  TrainingArrayType m_trainingSignedRangeInverse;
  LabelVectorType  m_trainingLabels;
  unsigned long m_trainingLabelsSize;
  //FLANNIndexType m_flannIndex;
  //FLANNDatasetType m_flannDataset
};

#endif // LesionSegmentationModel_h
