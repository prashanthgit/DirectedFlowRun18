
/*
  This version of PicoAnalyzer is only for the flow-decorrelation stuff
  I got rid of HUGE AMOUNT of other very useful stuff (keeping it in a
  copied version of this file), just to make things simpler (and faster).
  - MAL 16june2018
*/


#include "tools/PicoAnalyzer.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBbcHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StBbcGeom.h"

#include "TChain.h"
//#include "TTree.h"
//#include "TLeaf.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
using namespace std;


double PicoAnalyzer::GetBbcPmtPhi(short PmtId){
  unsigned short nBbcTiles;        // how many (1 or 2) tiles are associated with a given BBC phototube
  unsigned short tileNumbers[2];   // what are the tile ids of tiles associated with a given BBC phototube
  int chooseOne=0;
  // here we go with the bloody nightmare of the BBC's "shared phototubes"....
  mBbcGeom->GetTilesOfPmt(abs(PmtId),&nBbcTiles,tileNumbers);
  if (nBbcTiles>1) chooseOne = (mRan->Rndm()<0.5)?0:1;
  double phi = mBbcGeom->TileCenter(tileNumbers[chooseOne]*PmtId/abs(PmtId)).Phi();
  return phi;
}




//=================================================
PicoAnalyzer::PicoAnalyzer() : mNmipQtB(115), mNmipQtC(160), mnMipThreshold(0.3), mPicoDst(0), mEpdHits(0), mBbcHits(0), mEventClonesArray(0), mQ1NtupleFile(0), mQ2NtupleFile(0), mRunId(0), mRunCollisionSystem(0)
{
  mEpdGeom = new StEpdGeom;
  mBbcGeom = new StBbcGeom;
  mRan = new TRandom3;
  mRan->GetSeed();
}

//=================================================
PicoAnalyzer::~PicoAnalyzer(){
  /* no-op */
}


//=================================================
//void PicoAnalyzer::SetPicoDst(TTree* PicoDst){
void PicoAnalyzer::SetPicoDst(TChain* PicoDst){
  mPicoDst        = PicoDst;
  //  mEpdBranch      = mPicoDst->GetBranch("EpdHit");
  //  mBbcBranch      = mPicoDst->GetBranch("BbcHit");
  //  mEventBranch    = mPicoDst->GetBranch("Event");

  mEpdHits = new TClonesArray("StPicoEpdHit");
  mBbcHits = new TClonesArray("StPicoBbcHit");
  mTracks  = new TClonesArray("StPicoTrack");
  mEventClonesArray = new TClonesArray("StPicoEvent");


  mPicoDst->SetBranchStatus("*",0);         // turns OFF all branches  (speeds it up :-)
  unsigned int found;
  mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  cout << "EpdHit Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);

  mPicoDst->SetBranchStatus("Event*",1,&found);
  cout << "Event Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Event",&mEventClonesArray);

  /*
     mPicoDst->SetBranchStatus("BbcHit*",1,&found);
     cout << "BbcHit Branch returned found= " << found << endl;
     mPicoDst->SetBranchAddress("BbcHit",&mBbcHits);
     */

  mPicoDst->SetBranchStatus("Track*",1,&found);
  cout << "Track Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Track",&mTracks);

  /*
     IMPORTANT!  When you use TTrees, you can grab the TBranches and set the address
     But when you are using a TChain, then you have to use TChain::SetBranchAddress as I do above
     mEpdBranch->SetAddress(&mEpdHits);
     mBbcBranch->SetAddress(&mBbcHits);
     mEventBranch->SetAddress(&mEventClonesArray);
     */
  cout << "I have set all branches, leaves, etc \n";
}

//=================================================
short PicoAnalyzer::Init(){
  mHistoFile = new TFile("EpdHists.root","RECREATE");

  // 1D histograms for nMIP distributions
  for (int ew=0; ew<2; ew++){
    for (int pp=1; pp<13; pp++){
      for (int tt=1; tt<32; tt++){
	mNmipDists[ew][pp-1][tt-1] = new TH1D(Form("NmipEW%dPP%dTT%d",ew,pp,tt),Form("NmipEW%dPP%dTT%d",ew,pp,tt),500,0,20);
	mAdcDists[ew][pp-1][tt-1]  = new TH1D(Form("AdcEW%dPP%dTT%d",ew,pp,tt),Form("AdcEW%dPP%dTT%d",ew,pp,tt),512,0,4096);
      }
    }
  }

  //ReadInSystems();  // read in a text file of what system is being collided in what run.

  cout << "Init done\n";
  return 0;
}

//==============================================
// some things might need resetting when there is a new run
void PicoAnalyzer::NewRun(int runId){
  mRunCollisionSystem = WhichSystem(runId);
  mRunId = runId;
}



//=================================================
short PicoAnalyzer::Make(int iEvent){

  //----------------- get data --------------
  mPicoDst->GetEntry(iEvent);
  StPicoEvent* event = (StPicoEvent*)((*mEventClonesArray)[0]);
  if (event->runId()!=mRunId){
    NewRun(event->runId());        // some things should be reset when there is a new run loaded
    cout << "New run detected: " << mRunId << " and it is collision system #" << mRunCollisionSystem << endl;
  }
  //----- done getting data; have fun! ------

  StThreeVectorF primaryVertex = event->primaryVertex();
  TVector3 PV(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());

  double pi = TMath::Pi();

  int mRefMult = event->refMult();
  int mRunId = event->runId();
  double mVzVpd = event->vzVpd() + 36.4046 - 3.13415;  // numbers from Rosi's email of 3may2018
  int CentId = FindCent(mRefMult);   // returns an integer between 0 (70-80%) and 8 (0-5%)
  int RunYear = 19;
  int RunDay = floor( (mRunId - RunYear*pow(10,6))/pow(10,3) );
  int DayBinId = RunDay-89;

  if (CentId<0) return 0;            // 80-100% - very peripheral
  if (CentId<5) return 0;  // I only want to look at 0-30%
  if (fabs(PV.Z())>30.0) return 0;
  if (sqrt(pow(PV.X(),2)+pow(PV.Y(),2))>1.0) return 0;

  StPicoEpdHit* epdHit;

  double Q[2][3][2] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // indices are east/west, order (1, 2 or 3), component (x or y)
  double nMipSum(0.0);
  double maxWeight = 2;


  double QxRing[2][16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double QyRing[2][16] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

  // these weights optimized to midcentral Au+Au collisions at 27 GeV
  // ---> THEY SHOULD BE USED FOR FIRST-ORDER PLANE ONLY!!
  // https://drupal.star.bnl.gov/STAR/blog/lisa/phi-weighting-and-optimizing-ring-weights-auau-27-gev
  //  double RingWeight[16] = {+1.00000e+00, +6.76445e-01, +4.43371e-01, +3.13981e-01, +2.85431e-02, -6.11067e-02, -1.74905e-01, -1.62882e-01,
  //			   -2.17714e-02, -1.74471e-01, -1.52217e-01, -1.80971e-01, -1.42348e-01, -2.08882e-01, -1.20629e-01, -1.15527e-01};

  double QsubEventTpc[4][2][mNumberOfTpcSubEvents];
  for (int order=1; order<=4; order++){
    for (int xy=0; xy<2; xy++){
      for (int sub=0; sub<mNumberOfTpcSubEvents; sub++){
	QsubEventTpc[order-1][xy][sub] = 0.0;
      }
    }
  }

  double QsubEventEpd[2][4][3][mNumberOfEpdSubEvents];      // the third index is xy.  xy=0 means x-component (weight*cos); xy=1 means y-component (weight*sin); xy=2 means normalization (weight)
  for (int ew=0; ew<2; ew++){
    for (int order=1; order<=4; order++){
      for (int xy=0; xy<3; xy++){
	for (int sub=0; sub<mNumberOfEpdSubEvents; sub++){
	  QsubEventEpd[ew][order-1][xy][sub] = 0.0;
	}
      }
    }
  }


  //------------------------------------------------------
  //--------------- Begin loop over EPD hits --------------
  for (int hit=0; hit<mEpdHits->GetEntries(); hit++){
    epdHit = (StPicoEpdHit*)((*mEpdHits)[hit]);
    int ew = (epdHit->id()<0)?0:1;
    int pp = epdHit->position();
    int tile = epdHit->tile();
    double adc = epdHit->adc();
    //std::cout<<ew<<"\t"<<pp<<"\t"<<tile<<"\t"<<adc<<std::endl;
    double nMip = (epdHit->tile()<10)?epdHit->adc()/mNmipQtC:epdHit->adc()/mNmipQtB;
    if (nMip<mnMipThreshold) continue;
    double weight = (nMip>maxWeight)?maxWeight:nMip;
    FillPhiWeightHistos(epdHit,weight);  // calculates phi-weights for NEXT use.  (Obviously do this BEFORE phi-weighting here!)
    //    weight /= mPhiWeightInput[ew]->GetBinContent(epdHit->position(),epdHit->tile());     // note!  Phi weighting!
    //    TVector3 pos = mEpdGeom->RandomPointOnTile(epdHit->id());
    //    TVector3 StraightLine = pos - PV;
    TVector3 CenterOfTile = mEpdGeom->TileCenter(epdHit->id());
    //    double phi = CenterOfTile.Phi();
    //    double eta = StraightLine.Eta();

    TVector3 StraightLineToTileCenter = CenterOfTile - PV;
    double phi             = StraightLineToTileCenter.Phi();
    double etaOfTileCenter = StraightLineToTileCenter.Eta();

    int row = epdHit->row();

    if (fabs(etaOfTileCenter)<mEpdEtaBbounds[0]) continue;
    int EpdSubEventBin = -99;
    for (int iEpdBin=1; iEpdBin<=mNumberOfEpdSubEvents; iEpdBin++){
      if (fabs(etaOfTileCenter)<mEpdEtaBbounds[iEpdBin]){
	EpdSubEventBin = iEpdBin-1;
	break;
      }
    }
    if (EpdSubEventBin<0){
      continue;
    }

  }

  return 0;
}

//=================================================
short PicoAnalyzer::Finish(){
  cout << "Finish!!\n\n";
  return 0;
}

//-----------------------------------
// I copied this directly from Isaac
// it's for the 2018 isobars
int PicoAnalyzer::FindCent(int Multiplicity){
  int CentId = -1;
  if(Multiplicity <= 19) CentId = -1; // > 80%                                 
  else if (Multiplicity > 19 && Multiplicity <= 31) CentId = 0; // 70-80%      
  else if (Multiplicity > 31 && Multiplicity <= 46) CentId = 1; // 60-70%      
  else if (Multiplicity > 46 && Multiplicity <= 67) CentId = 2; // 50-60%      
  else if (Multiplicity > 67 && Multiplicity <= 93) CentId = 3; // 40-50%      
  else if (Multiplicity > 93 && Multiplicity <= 128) CentId = 4; // 30-40%     
  else if (Multiplicity > 128 && Multiplicity <= 172) CentId = 5; // 20-30%    
  else if (Multiplicity > 172 && Multiplicity <= 230) CentId = 6; // 10-20%    
  else if (Multiplicity > 230 && Multiplicity <= 267) CentId = 7; // 5-10%     
  else if (Multiplicity > 267) CentId = 8; // 0-5%                             
  return CentId;
}



/*  
    To make the file that is read in by the method below, there are two ways:

    First, from Dmitry:
    How-To: invoke RunLog, select appropriate filters ("physics" - check, "tpx" 
    - check, "filter bad runs": check three times etc), click "select". Then click 
    on the number of runs displayed below the select button: "Runs: 3513" (click on 
    number). In the middle pane you will get an ascii list like this:

    ...
    19125030 1525539034 Ru 100.00 Ru 100.00
    19125029 1525537184 Ru 100.00 Ru 100.00
    19125028 1525535327 Ru 100.00 Ru 100.00
    19125027 1525534083 Ru 100.00 Ru 100.00
    19125022 1525526621 Zr 100.00 Zr 100.00
    19125021 1525524762 Zr 100.00 Zr 100.00
    19125020 1525522902 Zr 100.00 Zr 100.00
    19125019 1525520656 Zr 100.00 Zr 100.00
    19125018 1525518764 Zr 100.00 Zr 100.00
    ...
columns: run number, run start time (unixtime), blue beam species, blue energy, yellow beam species, yellow energy


Next, from Gene, though I didn't try it:
On RCF:
mysql -h heston.star.bnl.gov -P 3501 -C RunLog -e "select runNumber,blueSpecies from beamInfo"
Output that to a file and modify as you please.
*/



void PicoAnalyzer::ReadInSystems(){
  // 2018 run ("run 19") colliding system:  1=Au27Au / 2=Ru200Ru / 3=Zr200Zr / 4-Au27fixedTarget
  for (int day=0; day<365; day++){
    for (int runInDay=0; runInDay<500; runInDay++){
      mCollidingSystem[day][runInDay] = 0;
    }
  }
  int runId,junk,iSys;
  int iReadIn(0);
  ifstream RunSystems("Run2018systems.txt");
  if (!RunSystems){
    cout << "Error opening systems file!\n";
    return;
  }
  do{
    RunSystems >> runId >> junk >> iSys;
    if (RunSystems.eof()) break;
    if (runId!=0){
      iReadIn++;
      mCollidingSystem[(runId%1000000)/1000][runId%1000] = iSys;
    }
    else{
      RunSystems.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }while(kTRUE);

  cout << "I have learned the colliding system for " << iReadIn << " runs\n";
  RunSystems.close();
}

//------------------------------------------------------
short PicoAnalyzer::WhichSystem(int runId){
  if (runId/1000000 != 19){
    cout << "Sorry, I only know about 2018 run.  Do some coding!\n";
    return 0;
  }
  short shColSys = mCollidingSystem[(runId%1000000)/1000][runId%1000];
  return shColSys;
}

void PicoAnalyzer::FillPhiWeightHistos(StPicoEpdHit* epdHit, double weight){
  int ew = (epdHit->id()<0)?0:1;

  /* 
     mPhiWeightOutput are the histograms that will be used for Phi-Weighting in the next run
     (mPhiAveraged is just used for normalization of mPhiWeightOutput at the end, then it is deleted)
     These histograms are binned in PP and TT, so GetBinContent(10,4) gives the weight for PP10TT04.
     */
  //  mPhiWeightOutput[ew]->Fill(epdHit->position(),epdHit->tile(),weight);
  if (epdHit->tile()==1){
    //  for (int pp=1; pp<13; pp++) mPhiAveraged[ew]->Fill(pp,1,weight/12.0);
  }
  else{
    for (int pp=1; pp<13; pp++){
      // for (int tt=2*(epdHit->row()-1); tt<2*epdHit->row(); tt++) mPhiAveraged[ew]->Fill(pp,tt,weight/24.0);
    }
  }

  /*
     The histograms below are the same as the ones above, except that they are binned
     in row and phi.  This makes it really nice to make plots in polar mode.

     Unfortunately, unless you want rotatable lego/surface plots (which probably you don't),
     the way you have to draw those "polar" histograms is as follows:
     gPad->DrawFrame(-17,-17,17,17);
     EpdWeightWestPolar->Draw("colz,pol,same");    
     */
  double phi = mEpdGeom->TileCenter(epdHit->id()).Phi();
  double phiDraw = (phi>0.0)?phi:phi+2.0*TMath::Pi();
  double dPhi = 2.0*TMath::Pi()/24.0;
  int row = epdHit->row();
  if (row==1){      // special case, tile 1.  Display it as two "tiles"
    //mHisto2D[12+ew]->Fill(phiDraw+dPhi/2.0,row,weight);
    //mHisto2D[12+ew]->Fill(phiDraw-dPhi/2.0,row,weight);
    // for (double phiTemp=dPhi/2.0; phiTemp<2.0*TMath::Pi(); phiTemp+=dPhi) mHisto2D[14+ew]->Fill(phiTemp,row,weight/12.0);
  }
  else{
    //mHisto2D[12+ew]->Fill(phiDraw,row,weight);    
    // for (double phiTemp=dPhi/2.0; phiTemp<2.0*TMath::Pi(); phiTemp+=dPhi) mHisto2D[14+ew]->Fill(phiTemp,row,weight/24.0);
  }

}

double PicoAnalyzer::v1Weight(int CentId, double eta){
  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
  double cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};

  double etaMagnitude = fabs(eta);   // we have the same weights for both wheels.  Any sign flips (East versus West etc) need to be done in the main code.

  return lin[CentId]*etaMagnitude + cub[CentId]*pow(etaMagnitude,3);
}
