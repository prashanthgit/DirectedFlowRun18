// Accessing individual TLeaf objects doesn't (easily) work with TChain objects
// the address changes every time a new file is loaded (which makes sense if you think of it)
// should stick with TTrees.
// https://root-forum.cern.ch/t/reading-exactly-1-leaf-from-tchain/14129/3

#include <string>
#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "tools/PicoAnalyzer.h"


/*
void AnalyzePico(int nRootFilesToProcess=1, std::string FileListName = "DataFile.list"){

  TString rootFileName;
  TTree* picoDst;
  PicoAnalyzer* anal = new PicoAnalyzer();
  anal->Init();

  std::ifstream ifs(FileListName,std::ifstream::in);
  if (!ifs){std::cout << "File list " << FileListName << " does not exist!  I quit.\n ";  return;}
  for (int iFile=0; iFile<nRootFilesToProcess; iFile++){
    ifs >> rootFileName;
    if (!(ifs.good())) break;
    std::cout << "Reading data from " << rootFileName << " (File " << iFile << " of " << nRootFilesToProcess << ")" << std::endl;
    TFile* inFile = new TFile(rootFileName);
    if(inFile->IsZombie()){
       std::cout << rootFileName << " is a zombie!!!" << std::endl;
      delete inFile;
      continue;
    }
    inFile->GetObject("PicoDst",picoDst);
    anal->SetPicoDst(picoDst);              // this must be done for each new file.  that's why we don't use TChain
    for (int ievent=0; ievent<picoDst->GetEntries(); ievent++){
      anal->Make(ievent);
    }
    inFile->Close();
    delete inFile;
  }
  anal->Finish();
}
*/

void AnalyzePico(int nRootFilesToProcess=1, std::string FileListName = "DataFile.list"){


  std::cout<<" file name of list of files "<<FileListName<<std::endl;
  std::ifstream ifs(FileListName,std::ifstream::in);
  if (!ifs){std::cout << "File list " << FileListName << " does not exist!  I quit.\n ";  return;}


  TChain* picoDst = new TChain("PicoDst");
  PicoAnalyzer* anal = new PicoAnalyzer();
  anal->Init();
  TString rootFileName;

  for (int iFile=0; iFile<nRootFilesToProcess; iFile++){
    ifs >> rootFileName;
    std::cout << iFile<<"\t"<<rootFileName<<std::endl;
    if (!(ifs.good())) break;
    picoDst->Add(rootFileName.Data());
  }
  anal->SetPicoDst(picoDst);              // this must be done for each new file.  that's why we don't use TChain
  int NeventsToAnalyze = picoDst->GetEntries();
  std::cout << "Preparing to analyze " << NeventsToAnalyze << " events...\n";
  for (int ievent=0; ievent<picoDst->GetEntries(); ievent++){
    if (ievent%1000==0) std::cout << "Processed " << ievent << " events - " << 100.0*(double)ievent/(double)NeventsToAnalyze << "% \n";
    anal->Make(ievent);
  }
  anal->Finish();
  
}
  


  

int main(int argc, char** argv)
{
  int nFilesToProcess=1;                 // default
  if (argc>1) nFilesToProcess = atoi(argv[1]);

  TString FileList = "DataFile.list";    // default
  if (argc>2) FileList = argv[2];
  //  AnalyzePico(1, "DataFile.list");
  AnalyzePico(nFilesToProcess, FileList.Data());

  return EXIT_SUCCESS;
}
