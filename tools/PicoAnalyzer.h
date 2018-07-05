#ifndef PicoAnalyzer__
#define PicoAnalyzer__

#include "TObject.h"

class TChain;
class TClonesArray;
class StPicoEvent;
class TH1D;
class TH2D;
class TProfile;
class TProfile2D;
class TNtuple;
class TFile;
class StEpdGeom;
class StBbcGeom;
class TRandom3;
class StPicoEpdHit; 

class PicoAnalyzer : public TObject {
 public:
  PicoAnalyzer();
  ~PicoAnalyzer();

  void SetPicoDst(TChain*);
  short Init();
  short Make(int iEvent);
  short Finish();

 private:

  // internal methods
  int FindCent(int RefMult);   // utility class just giving centrality bin.  Copied directly from Isaac 1 May 2018
  double GetBbcPmtPhi(short PmtId);
  void ReadInSystems();    // reads a text file that idenfies the collision system for every run
  short WhichSystem(int runId);
  void NewRun(int runId);    // invoked when Make() detects that a new run has been loaded
  void FillPhiWeightHistos(StPicoEpdHit* epdHit, double weight);     // fills the histograms used for (a later job's) phi weighting

  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double v1Weight(int CentId, double eta);

  // useful objects kept by PicoAnalyzer
  StEpdGeom* mEpdGeom;
  StBbcGeom* mBbcGeom;
  TRandom3* mRan;
  char mCollidingSystem[365][500];    // index1=day of year;  index2=run of day
  int mRunId;                         // when this changes, refresh some information.
  short mRunCollisionSystem;

  // the data objects
  TChain*   mPicoDst;
  TClonesArray* mEpdHits;
  TClonesArray* mBbcHits;
  TClonesArray* mTracks;
  TClonesArray* mEventClonesArray;  // kind of hilarious that the StPicoEvent is stored as a one-element TClonesArray :-)

  // parameters relevant to my analysis
  double mNmipQtB;  // ADC value of MIP peak on rings 6-16 (read out thru QT32Bs)
  double mNmipQtC;  // ADC value of MIP peak on rings 1-5  (read out thru QT32Cs)
  double mnMipThreshold;  // low-signal threshold, to cut out noise basically.



  //  static const int mNumberOfEpdSubEvents = 6;
  //  double mEpdEtaBbounds[mNumberOfEpdSubEvents+1] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  static const int mNumberOfEpdSubEvents = 3;
  double mEpdEtaBbounds[mNumberOfEpdSubEvents+1] = {2.0, 3.0, 4.0, 5.0};
  static const int mNumberOfTpcSubEvents = 16;

  TProfile* mEastLongDecorProfile[mNumberOfEpdSubEvents][4];    // second index is order of event plane
  TProfile* mWestLongDecorProfile[mNumberOfEpdSubEvents][4];


  // ---------- Now, my histograms, ntuples, TFiles, etc.  All stuff particular to my analysis
  // 1D histograms
  TH1D* mNmipDists[2][12][31];   // nMIP distributions for all tiles
  TH1D* mAdcDists[2][12][31];   // ADC distributions for all tiles

  // ntuples
  TNtuple* mQ1vectorNtuple;      // Q1 vectors ring-by-ring. For offline weight optimization
  TNtuple* mQ2vectorNtuple;      // Q2 vectors ring-by-ring. For offline weight optimization

  // TFiles (to store histograms and data)
  TFile* mHistoFile;
  TFile* mQ1NtupleFile;
  TFile* mQ2NtupleFile;
  TFile* mCorrectionInputFile;
  TFile* mCorrectionOutputFile;

  TFile* mTpcWeightFile;      // just holds eta-phi 2D histogram
  TFile* mTpcWeightInputFile;
};


#endif
