#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/VectorUtil.h"
 
#define USE_VAR_BASELINE

#define USE_LEP1 
#define USE_LEP2
#define USE_JETS 
#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_GLOBAL_VAR
#define USE_EVT_VAR 
#define USE_GEN_LEP
#define ROOT_Math_GenVector_VectorUtil

#include "../common/TFFactory.h"
// #include "goodrun.cc"
#include "../Selection/moriond.h"

using namespace std;
//using namespace theDoctor;



// ----------------------------------------------
// Should be called only here because many
// struct and fuctions have to be declare first
// ----------------------------------------------
#include "../sonicScrewdriver/interface/BabyScrewdriver.h"

uint32_t counter = 0;
string empty = "";
string storedDataset = "";
TH2D *h2 = NULL;
TH3D *h3 = NULL;
TAxis *xaxis = NULL;
TAxis *yaxis = NULL;
bool checkNegativeYields = false;
uint32_t nthentry = 0;
string outputName = "";
double scale1fbS2 =1;
 

	
//TH2 *moyenne2 = new TH2D("moyenne2", "moyenne2", 40, 200, 1200, 26, 0, 650);



 
//float getWeight(string currentProcessType, float lumi, float s1fb2=0);
double getWeight(string currentProcessType, float lumi, double s1fb2=1,bool mediumbtag=true);
map< pair<uint32_t,uint32_t>, string > scanMap;

TFile *fileX = new TFile("../common/xsec_stop_13TeV.root");
TH1D* stopXSEC = (TH1D*)fileX->Get("stop")->Clone();

bool lepChannel() 
{ 
    return true; 
} 
    
//Add this as a global variable 

void BabyScrewdriver::Init()
{

  PrintBoxedMessage("Initializing babyScrewdriver");
  
  
  //path where all the files are
  babyTuplePath = "/opt/sbg/cms/ui3_data1/mjansova/Stop1lSharedBabies/v22/skim/";
  
  //for paralelization
  //    totalNumberOfWorkers = 1;
  totalNumberOfWorkers = 5;
  
  
  //creation of the histogram
  TFile *ftmp = NULL;
  TH2D *htmp = NULL;
  TString fNameTmp =  babyTuplePath+"signal_T2tt_forBinning.root"; //stop
  ftmp = new TFile(fNameTmp);
  htmp = (TH2D*)ftmp->Get("histNEvts")->Clone();
  TH1D* pX = htmp->ProjectionX(); 
  TH1D* pY = htmp->ProjectionY();

    //Add the variables needed to plot and to make the cut for the region the definition is on the ../common/Reader_CommonFormat_CommonBabies.h document, it makes the reation between the variable in the root plot and the variable in the code
  AddVariable("MET", "MET",  "MET", 4 ,250,650,  &(myEvent.pfmet), "noUnderflowInFirstBin");
  AddVariable("MT2W", "MT2W",  "MT2W", 10 ,0,500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
  AddVariable("MT", "MT",  "MT", 20 ,0,1000,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
  AddVariable("nJets","nJets","nJets",10,1,10,&(myEvent.ngoodjets),"noUnderflowInFirstBin");
   AddVariable("nBJets","nBJets","nBJets",5,1,5,&(myEvent.ngoodbtags),"noUnderflowInFirstBin");
  AddVariable("topnessMod","topnessMod","topnessMod",20,-20,20,&(myEvent.topnessMod),"noUnderflowInFirstBin");
  AddVariable("dphi","dphi","dphi", 7,0,3.5,&(myEvent.dphi_ak4pfjets_met),"noUnderflowInFirstBin");//angle between the leading jet and the MET
  AddVariable("Mlb","Mlb","Mlb", 10,0, 500,&(myEvent.Mlb),"noUnderflowInFirstBin");
  AddVariable("lep1_pt","lep1_pt","lep1_pt",10,0,500,&(myEvent.lep1_pt),"noUnderflowInFirstBin");
  AddVariable("StopMass","stop mass", "GeV", pX->GetNbinsX(), pX->GetBinLowEdge(1) ,  (pX->GetBinLowEdge( pX->GetNbinsX())) + (pX->GetBinWidth(pX->GetNbinsX())) ,&myEvent.mass_stop);
  AddVariable("NeutralinoMass","lsp mass", "GeV", pY->GetNbinsX(), pY->GetBinLowEdge(1), (pY->GetBinLowEdge( pY->GetNbinsX())) + (pY->GetBinWidth(pY->GetNbinsX())) , &myEvent.mass_lsp);
  // double *jetB = &(myEvent.ak4pfjets_leadbtag_p4.Pt());
  //cout<<&(myEvent.ak4pfjets_leadbtag_p4.Pt())<<end;
  AddVariable("jetB_pt","jetB_pt","jetB_pt",10,0,500,&(jetB) ,"noUnderflowInFirstBin"); //Pt of the jet that is identified as a B jets
  AddVariable("pt_lead_jet","pt_lead_jet","pt_lead_jet",16,0,800,&(myEvent.ak4pfjets_pt),"noUnderflowInFirstBin"); //Pt of the leading jet
  AddVariable("jet_b_indetif","jet_b_indetif","jet_b_indetif",10,0,1,&(myEvent.ak4pfjets_leadJet_identification),"noUnderflowInFirstBin");//CSV variable for the leading jet, usefull for te study o the compressed regions
  AddVariable("dphiLep","dphiLep","dphiLep", 7,0,3.5,&(myEvent.lep1_dphiMET),"noUnderflowInFirstBin");// angle between the lepton and the MET
 

  // -----------------
  // Datasets
  // ------------------
  //data
  /*   AddProcessClass("data1", "data1", "data", 1); ///@MJ@ TODO discared this!!
       AddDataset("data_met_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
       AddDataset("data_met_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
       AddDataset("data_met_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
       AddDataset("data_met_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
       AddDataset("data_met_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
       AddDataset("data_met_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
       AddDataset("data_met_Run2016H_MINIAOD_03Feb2017_ver2-v1","data1",0,0);
       AddDataset("data_met_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);

        AddDataset("data_single_muon_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_single_muon_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016H_MINIAOD_03Feb2017_ver2-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);

	AddDataset("data_single_electron_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
	AddDataset("data_single_electron_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
	AddDataset("data_single_electron_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
	AddDataset("data_single_electron_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
	AddDataset("data_single_electron_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
	AddDataset("data_single_electron_Run2016H_MINIAOD_03Feb2017_ver2-v1","data1",0,0);
	AddDataset("data_single_electron_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);
  */
  
    
  
    //---------------------------
    ///background simulation with background separated
	    //---------------------------
   AddProcessClass("bkgOneLepFromW", "bkgOneLepFromW", "background", kBlue); 
  AddDataset("W1JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0); 
  AddDataset("W2JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
  AddDataset("W3JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
  AddDataset("W4JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
        
  AddDataset("W1JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
  AddDataset("W2JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
  AddDataset("W3JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
  AddDataset("W4JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
  AddDataset("WWToLNuQQ_powheg_25ns","bkgOneLepFromW",0,0);
  
  AddProcessClass("bkgOneLepFromTop", "bkgOneLepFromTop", "background", kRed);
  AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
  AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); 
  AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
  AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); 
  
  AddProcessClass("bkgLostLepton", "bkgLostLepton", "background", kGreen+2); 
  AddDataset("ttbar_diLept_madgraph_pythia8_25ns","bkgLostLepton",0,0);
  AddDataset("ttbar_diLept_madgraph_pythia8_ext1_25ns","bkgLostLepton",0,0);
  AddDataset("tbar_tch_4f_powheg_pythia8_inclDecays_25ns","bkgLostLepton",0,0);
  AddDataset("t_sch_4f_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
  AddDataset("t_tW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
  AddDataset("t_tbarW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
  AddDataset("WWTo2l2Nu_powheg_25ns","bkgLostLepton",0,0);
  AddDataset("ttWJets_13TeV_madgraphMLM","bkgLostLepton",0,0);
  
  AddProcessClass("bkgZnunu", "bkgZnunu", "background", kViolet);
  AddDataset("ZZTo2L2Nu_powheg_pythia8_25ns","bkgZnunu",0,0);
  AddDataset("ttZJets_13TeV_madgraphMLM","bkgZnunu",0,0);
  AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","bkgZnunu",0,0);
 

 

  
  //----------------------
  //Background all in one in order to plot with only one background
  /*  AddProcessClass("bkgTotal", "bkgTotal", "background", kRed);
  AddDataset("W1JetsToLNu_madgraph_pythia8_25ns","bkgTotal",0,0); 
  AddDataset("W2JetsToLNu_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("W3JetsToLNu_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("W4JetsToLNu_madgraph_pythia8_25ns","bkgTotal",0,0);
  
  AddDataset("W1JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("W2JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("W3JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("W4JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("WWToLNuQQ_powheg_25ns","bkgTotal",0,0);
  
  AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns","bkgTotal",0,0); 
  AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns","bkgTotal",0,0); 
  
  AddDataset("ttbar_diLept_madgraph_pythia8_25ns","bkgTotal",0,0);
  AddDataset("ttbar_diLept_madgraph_pythia8_ext1_25ns","bkgTotal",0,0);
  AddDataset("tbar_tch_4f_powheg_pythia8_inclDecays_25ns","bkgTotal",0,0);
  AddDataset("t_sch_4f_amcnlo_pythia8_25ns","bkgTotal",0,0);
  AddDataset("t_tW_5f_powheg_pythia8_noHadDecays_25ns","bkgTotal",0,0);
  AddDataset("t_tbarW_5f_powheg_pythia8_noHadDecays_25ns","bkgTotal",0,0);
  AddDataset("WWTo2l2Nu_powheg_25ns","bkgTotal",0,0);
  AddDataset("ttWJets_13TeV_madgraphMLM","bkgTotal",0,0);
  
  AddDataset("ZZTo2L2Nu_powheg_pythia8_25ns","bkgTotal",0,0);
  AddDataset("ttZJets_13TeV_madgraphMLM","bkgTotal",0,0);
  AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","bkgTotal",0,0);
  */


  

  //------------------  
  //stop simulation of the signal
  //------------------
  //test for multiple mass for the plots, not to use if you want only one point
  /*   AddProcessClass("2Chargino_800_600", "2Chargino_800_600", "signal", kOrange+4);
  AddDataset("Signal_T2bW_1","2Chargino_800_600",0,0);
  AddDataset("Signal_T2bW_2","2Chargino_800_600",0,0);
  AddDataset("Signal_T2bW","2Chargino_800_600",0,0);*/
  
  /*   AddProcessClass("mixt_800_600", "mixt_800_600", "signal", kGreen+8);
  AddDataset("Signal_T2tb","mixt_800_600",0,0);
  AddDataset("Signal_T2tb_1","mixt_800_600",0,0); */
  
 
  /*  AddProcessClass("NoChargino_800_600", "NoChargino_800_600", "signal", kViolet+8);
  AddDataset("Signal_T2tt_mStop_150to250", "NoChargino_800_600", 0, 0 );
  AddDataset("Signal_T2tt_mStop_250to350", "NoChargino_800_600", 0, 0 );
  AddDataset("Signal_T2tt_mStop_350to400", "NoChargino_800_600", 0, 0 );
  AddDataset("Signal_T2tt_mStop_350to400_1", "NoChargino_800_600", 0, 0 );
  AddDataset("Signal_T2tt_mStop_400to1200_1", "NoChargino_800_600", 0, 0 );
  AddDataset("Signal_T2tt_mStop_400to1200_2", "NoChargino_800_600", 0, 0 );
  AddDataset("Signal_T2tt_mStop_400to1200", "NoChargino_800_600", 0, 0 );*/



    AddProcessClass("2Chargino", "2Chargino", "signal", kOrange+4);
  AddDataset("Signal_T2bW_1","2Chargino",0,0);
  AddDataset("Signal_T2bW_2","2Chargino",0,0);
  AddDataset("Signal_T2bW","2Chargino",0,0);
  
    AddProcessClass("mixt", "mixt", "signal", kGreen+8);
  AddDataset("Signal_T2tb","mixt",0,0);
  AddDataset("Signal_T2tb_1","mixt",0,0); 
  
  
   AddProcessClass("NoChargino", "NoChargino", "signal", kViolet+8);
  AddDataset("Signal_T2tt_mStop_150to250", "NoChargino", 0, 0 );
  AddDataset("Signal_T2tt_mStop_250to350", "NoChargino", 0, 0 );
  AddDataset("Signal_T2tt_mStop_350to400", "NoChargino", 0, 0 );
  AddDataset("Signal_T2tt_mStop_350to400_1", "NoChargino", 0, 0 );
  AddDataset("Signal_T2tt_mStop_400to1200_1", "NoChargino", 0, 0 );
  AddDataset("Signal_T2tt_mStop_400to1200_2", "NoChargino", 0, 0 );
  AddDataset("Signal_T2tt_mStop_400to1200", "NoChargino", 0, 0 );

  //in order to separate multiple point on the final plot   
  // AddProcessClass("2Chargino_600_400", "2Chargino_600_400", "signal", kOrange+6);
  // AddProcessClass("2Chargino_300_100", "2Chargino_300_100", "signal", kOrange+4); 
  // AddProcessClass("2Chargino_800_600", "2Chargino_800_600", "signal", kOrange+2);


 
 // AddProcessClass("mixt_600_400", "mixt_600_400", "signal", kGreen+6); 
 //  AddProcessClass("mixt_300_100", "mixt_300_100", "signal", kGreen+4); 
  // AddProcessClass("mixt_800_600", "mixt_800_600", "signal", kGreen+2);
 

 
    //  AddProcessClass("NoChargino_600_400", "NoChargino_600_400", "signal", kViolet+6); 
    //  AddProcessClass("NoChargino_800_600", "NoChargino_800_600", "signal", kViolet+4); 
      //   AddProcessClass("NoChargino_300_100", "NoChargino_300_100", "signal", kViolet+2);
  
  
    
  
  /*     AddProcessClass("signal_500_325",  "T2tt (500/325)",             "signal",kViolet);
	 AddProcessClass("signal_800_400",  "T2tt (800/400)",             "signal",kBlue+3);
	 AddProcessClass("signal_1000_50",  "T2tt (1000/50)",             "signal",kBlue+6);
   */
  
  
    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;
    
    //-------------------------
    //region cut if you want to see the definition of all this region go in the ../Selection/moriond.h documents
    //-------------------------
    //region without cut
     AddRegion("baseline","baseline",&baseline);
     
     //addition of all regions
     AddRegion("Baseline","Baseline",&Baseline);
 
     
     //first basics cut
     AddRegion("SR1l","SR1l",&SR1l);
     
     //
     AddRegion("test","test", &test);
     
     
     //region in the 13TeV article from CMS
     AddRegion("A1","A1",&SR1l_A_250lessMETless350);
     AddRegion("A2","A2",&SR1l_A_350lessMETless450);
     AddRegion("A3","A3",&SR1l_A_450lessMETless600);
     AddRegion("A4","A4",&SR1l_A_600lessMETlessInf);
     AddRegion("B1","B1",&SR1l_B_250lessMETless450);
     AddRegion("B2","B2",&SR1l_B_450lessMETless600);
     AddRegion("B3","B3",&SR1l_B_600lessMETlessInf);
     AddRegion("C1","C1",&SR1l_C_250lessMETless350);
     AddRegion("C2","C2",&SR1l_C_350lessMETless450);
     AddRegion("C3","C3",&SR1l_C_450lessMETless550);
     AddRegion("C4","C4",&SR1l_C_550lessMETless650);
     AddRegion("C5","C5",&SR1l_C_650lessMETlessInf);
     AddRegion("D1","D1",&SR1l_D_250lessMETless350);
     AddRegion("D2","D2",&SR1l_D_350lessMETless450);
     AddRegion("D3","D3",&SR1l_D_450lessMETless550);
     AddRegion("D4","D4",&SR1l_D_550lessMETlessInf);
     AddRegion("E1","E1",&SR1l_E_250lessMETless350);
     AddRegion("E2","E2",&SR1l_E_350lessMETless550);
     AddRegion("E3","E3",&SR1l_E_550lessMETlessInf);
     AddRegion("F1","F1",&SR1l_F_250lessMETless450);
     AddRegion("F2","F2",&SR1l_F_450lessMETlessInf);
     AddRegion("G1","G1",&SR1l_G_250lessMETless350);
     AddRegion("G2","G2",&SR1l_G_350lessMETless450);
     AddRegion("G3","G3",&SR1l_G_450lessMETless600);
     AddRegion("G4","G4",&SR1l_G_600lessMETlessInf);
     AddRegion("H1","H1",&SR1l_H_250lessMETless450);
     AddRegion("H2","H2",&SR1l_H_450lessMETlessInf);
     AddRegion("I1","I1",&SR1l_I_250lessMETless350);
     AddRegion("I2","I2",&SR1l_I_350lessMETless450);
     AddRegion("I3","I3",&SR1l_I_450lessMETless550);
     AddRegion("I4","I4",&SR1l_I_550lessMETlessInf);
     
     /*
    
     //---------------------------
     /// test of new region for H for the 900,50
     //---------------------------
     // region with normal odified topness >10
     AddRegion("H3","H3",&SR1l_H_250lessMETless500);
     AddRegion("H4","H4",&SR1l_H_500lessMETlessInf);
     AddRegion("H5","H5",&SR1l_H_250lessMETless550);
     AddRegion("H6","H6",&SR1l_H_550lessMETlessInf);
     AddRegion("H7","H7",&SR1l_H_250lessMETless600);
     AddRegion("H8","H8",&SR1l_H_600lessMETlessInf);
     AddRegion("H9","H9",&SR1l_H_250lessMETless650);
     AddRegion("H10","H10",&SR1l_H_650lessMETlessInf);
     ////new region for H with 12 modified topness
     AddRegion("H1(12)","H1(12)",&SR1l_H_250lessMETless450_12);
     AddRegion("H2(12)","H2(12)",&SR1l_H_450lessMETlessInf_12);
     AddRegion("H3(12)","H3(12)",&SR1l_H_250lessMETless500_12);
     AddRegion("H4(12)","H4(12)",&SR1l_H_500lessMETlessInf_12);
     AddRegion("H5(12)","H5(12)",&SR1l_H_250lessMETless550_12);
     AddRegion("H6(12)","H6(12)",&SR1l_H_550lessMETlessInf_12);
     AddRegion("H7(12)","H7(12)",&SR1l_H_250lessMETless600_12);
     AddRegion("H8(12)","H8(12)",&SR1l_H_600lessMETlessInf_12);
     AddRegion("H9(12)","H9(12)",&SR1l_H_250lessMETless650_12);
     AddRegion("H10(12)","H10(12)",&SR1l_H_650lessMETlessInf_12);
    
    
    


    //-------------
    // test New regions for 600-325 2Chargino
    //------------
    AddRegion("E1B","E1B",&SR1l_E_250lessMETless350B);
    AddRegion("E2B","E2B",&SR1l_E_350lessMETless550B);
    AddRegion("G1B","G2B",&SR1l_G_250lessMETless350B);
    AddRegion("G2B","G2B",&SR1l_G_250lessMETless350B);
    AddRegion("I2B","I2B",&SR1l_I_350lessMETless450B);
    AddRegion("I3B","I3B",&SR1l_I_450lessMETless550B);


    //-------------
    //New regions for 300-50 2Chargino
    //------------
    AddRegion("CLOW","CLOW",&SR1l_C_150lessMETless250);
    AddRegion("DLOW","DLOW",&SR1l_D_150lessMETless250);
    AddRegion("ELOW","ELOW",&SR1l_E_150lessMETless250);
    AddRegion("ILOW","ILOW",&SR1l_I_150lessMETless250);
     */

    
    //-------------------------
    //New regions in order to include the compressed regions
    //-------------------------
     //Region for the R0 set of region to add to I1 -> I4
    AddRegion("C1CMPR","C1CMPR",&SR1l_C_250lessMETless350CMPR);
    AddRegion("C2CMPR","C2CMPR",&SR1l_C_350lessMETless450CMPR);
    AddRegion("C3CMPR","C3CMPR",&SR1l_C_450lessMETless550CMPR);
    AddRegion("C4CMPR","C4CMPR",&SR1l_C_550lessMETless650CMPR);
    AddRegion("C5CMPR","C5CMPR",&SR1l_C_650lessMETlessInfCMPR);
    AddRegion("D1CMPR","D1CMPR",&SR1l_D_250lessMETless350CMPR);
    AddRegion("D2CMPR","D2CMPR",&SR1l_D_350lessMETless450CMPR);
    AddRegion("D3CMPR","D3CMPR",&SR1l_D_450lessMETless550CMPR);
    AddRegion("D4CMPR","D4CMPR",&SR1l_D_550lessMETlessInfCMPR);
    AddRegion("E1CMPR","E1CMPR",&SR1l_E_250lessMETless350CMPR);
    AddRegion("E2CMPR","E2CMPR",&SR1l_E_350lessMETless550CMPR);
    AddRegion("E3CMPR","E3CMPR",&SR1l_E_550lessMETlessInfCMPR);
    AddRegion("F1CMPR","F1CMPR",&SR1l_F_250lessMETless450CMPR);
    AddRegion("F2CMPR","F2CMPR",&SR1l_F_450lessMETlessInfCMPR);
    AddRegion("G1CMPR","G1CMPR",&SR1l_G_250lessMETless350CMPR);
    AddRegion("G2CMPR","G2CMPR",&SR1l_G_350lessMETless450CMPR);
    AddRegion("G3CMPR","G3CMPR",&SR1l_G_450lessMETless600CMPR);
    AddRegion("G4CMPR","G4CMPR",&SR1l_G_600lessMETlessInfCMPR);
    AddRegion("H1CMPR","H1CMPR",&SR1l_H_250lessMETless450CMPR);
    AddRegion("H2CMPR","H2CMPR",&SR1l_H_450lessMETlessInfCMPR);
    

    //Region for the R1 set of region to add to the A to H from the paper 
    AddRegion("I108","I108",&SR1l_I_250lessMETless350_0_8);
    AddRegion("I208","I208",&SR1l_I_350lessMETless450_0_8);
    AddRegion("I308","I308",&SR1l_I_450lessMETless550_0_8);
    AddRegion("I408","I408",&SR1l_I_550lessMETlessInf_0_8);

    


    /*
    //for the mother W from Chargino, test for the origin of the lepton in the mixt scenario (T2tb)
    AddRegion("A1A","A1A",&SR1l_A_250lessMETless350A);
    AddRegion("A2A","A2A",&SR1l_A_350lessMETless450A);
    AddRegion("A3A","A3A",&SR1l_A_450lessMETless600A);
    AddRegion("A4A","A4A",&SR1l_A_600lessMETlessInfA);
    AddRegion("B1A","B1A",&SR1l_B_250lessMETless450A);
    AddRegion("B2A","B2A",&SR1l_B_450lessMETless600A);
    AddRegion("B3A","B3A",&SR1l_B_600lessMETlessInfA);
    AddRegion("C1A","C1A",&SR1l_C_250lessMETless350A);
    AddRegion("C2A","C2A",&SR1l_C_350lessMETless450A);
    AddRegion("C3A","C3A",&SR1l_C_450lessMETless550A);
    AddRegion("C4A","C4A",&SR1l_C_550lessMETless650A);
    AddRegion("C5A","C5A",&SR1l_C_650lessMETlessInfA);
    AddRegion("D1A","D1A",&SR1l_D_250lessMETless350A);
    AddRegion("D2A","D2A",&SR1l_D_350lessMETless450A);
    AddRegion("D3A","D3A",&SR1l_D_450lessMETless550A);
    AddRegion("D4A","D4A",&SR1l_D_550lessMETlessInfA);
    AddRegion("E1A","E1A",&SR1l_E_250lessMETless350A);
    AddRegion("E2A","E2A",&SR1l_E_350lessMETless550A);
    AddRegion("E3A","E3A",&SR1l_E_550lessMETlessInfA);
    AddRegion("F1A","F1A",&SR1l_F_250lessMETless450A);
    AddRegion("F2A","F2A",&SR1l_F_450lessMETlessInfA);
    AddRegion("G1A","G1A",&SR1l_G_250lessMETless350A);
    AddRegion("G2A","G2A",&SR1l_G_350lessMETless450A);
    AddRegion("G3A","G3A",&SR1l_G_450lessMETless600A);
    AddRegion("G4A","G4A",&SR1l_G_600lessMETlessInfA);
    AddRegion("H1A","H1A",&SR1l_H_250lessMETless450A);
    AddRegion("H2A","H2A",&SR1l_H_450lessMETlessInfA);
    AddRegion("I1A","I1A",&SR1l_I_250lessMETless350A);
    AddRegion("I2A","I2A",&SR1l_I_350lessMETless450A);
    AddRegion("I3A","I3A",&SR1l_I_450lessMETless550A);
    AddRegion("I4A","I4A",&SR1l_I_550lessMETlessInfA);
    




    //mother W from stop, test for the origin of the leptin in the mixt scenario  (T2tb)
    AddRegion("A1B","A1B",&SR1l_A_250lessMETless350B);
    AddRegion("A2B","A2B",&SR1l_A_350lessMETless450B);
    AddRegion("A3B","A3B",&SR1l_A_450lessMETless600B);
    AddRegion("A4B","A4B",&SR1l_A_600lessMETlessInfB);
    AddRegion("B1B","B1B",&SR1l_B_250lessMETless450B);
    AddRegion("B2B","B2B",&SR1l_B_450lessMETless600B);
    AddRegion("B3B","B3B",&SR1l_B_600lessMETlessInfB);
    AddRegion("C1B","C1B",&SR1l_C_250lessMETless350B);
    AddRegion("C2B","C2B",&SR1l_C_350lessMETless450B);
    AddRegion("C3B","C3B",&SR1l_C_450lessMETless550B);
    AddRegion("C4B","C4B",&SR1l_C_550lessMETless650B);
    AddRegion("C5B","C5B",&SR1l_C_650lessMETlessInfB);
    AddRegion("D1B","D1B",&SR1l_D_250lessMETless350B);
    AddRegion("D2B","D2B",&SR1l_D_350lessMETless450B);
    AddRegion("D3B","D3B",&SR1l_D_450lessMETless550B);
    AddRegion("D4B","D4B",&SR1l_D_550lessMETlessInfB);
    AddRegion("E1B","E1B",&SR1l_E_250lessMETless350B);
    AddRegion("E2B","E2B",&SR1l_E_350lessMETless550B);
    AddRegion("E3B","E3B",&SR1l_E_550lessMETlessInfB);
    AddRegion("F1B","F1B",&SR1l_F_250lessMETless450B);
    AddRegion("F2B","F2B",&SR1l_F_450lessMETlessInfB);
    AddRegion("G1B","G1B",&SR1l_G_250lessMETless350B);
    AddRegion("G2B","G2B",&SR1l_G_350lessMETless450B);
    AddRegion("G3B","G3B",&SR1l_G_450lessMETless600B);
    AddRegion("G4B","G4B",&SR1l_G_600lessMETlessInfB);
    AddRegion("H1B","H1B",&SR1l_H_250lessMETless450B);
    AddRegion("H2B","H2B",&SR1l_H_450lessMETlessInfB);
    AddRegion("I1B","I1B",&SR1l_I_250lessMETless350B);
    AddRegion("I2B","I2B",&SR1l_I_350lessMETless450B);
    AddRegion("I3B","I3B",&SR1l_I_450lessMETless550B);
    AddRegion("I4B","I4B",&SR1l_I_550lessMETlessInfB);
    

  
    //other mothers, test for the origin of the lepton in the mixt scenario (T2tb)
    AddRegion("A1C","A1C",&SR1l_A_250lessMETless350C);
    AddRegion("A2C","A2C",&SR1l_A_350lessMETless450C);
    AddRegion("A3C","A3C",&SR1l_A_450lessMETless600C);
    AddRegion("A4C","A4C",&SR1l_A_600lessMETlessInfC);
    AddRegion("B1C","B1C",&SR1l_B_250lessMETless450C);
    AddRegion("B2C","B2C",&SR1l_B_450lessMETless600C);
    AddRegion("B3C","B3C",&SR1l_B_600lessMETlessInfC);
    AddRegion("C1C","C1C",&SR1l_C_250lessMETless350C);
    AddRegion("C2C","C2C",&SR1l_C_350lessMETless450C);
    AddRegion("C3C","C3C",&SR1l_C_450lessMETless550C);
    AddRegion("C4C","C4C",&SR1l_C_550lessMETless650C);
    AddRegion("C5C","C5C",&SR1l_C_650lessMETlessInfC);
    AddRegion("D1C","D1C",&SR1l_D_250lessMETless350C);
    AddRegion("D2C","D2C",&SR1l_D_350lessMETless450C);
    AddRegion("D3C","D3C",&SR1l_D_450lessMETless550C);
    AddRegion("D4C","D4C",&SR1l_D_550lessMETlessInfC);
    AddRegion("E1C","E1C",&SR1l_E_250lessMETless350C);
    AddRegion("E2C","E2C",&SR1l_E_350lessMETless550C);
    AddRegion("E3C","E3C",&SR1l_E_550lessMETlessInfC);
    AddRegion("F1C","F1C",&SR1l_F_250lessMETless450C);
    AddRegion("F2C","F2C",&SR1l_F_450lessMETlessInfC);
    AddRegion("G1C","G1C",&SR1l_G_250lessMETless350C);
    AddRegion("G2C","G2C",&SR1l_G_350lessMETless450C);
    AddRegion("G3C","G3C",&SR1l_G_450lessMETless600C);
    AddRegion("G4C","G4C",&SR1l_G_600lessMETlessInfC);
    AddRegion("H1C","H1C",&SR1l_H_250lessMETless450C);
    AddRegion("H2C","H2C",&SR1l_H_450lessMETlessInfC);
    AddRegion("I1C","I1C",&SR1l_I_250lessMETless350C);
    AddRegion("I2C","I2C",&SR1l_I_350lessMETless450C);
    AddRegion("I3C","I3C",&SR1l_I_450lessMETless550C);
    AddRegion("I4C","I4C",&SR1l_I_550lessMETlessInfC);*/

    /*
    
    //lower cut for Delta Phi ->0.5 test for new cut with delta phi
    AddRegion("A105","A105",&SR1l_A_250lessMETless350_0_5);
    AddRegion("A205","A205",&SR1l_A_350lessMETless450_0_5);
    AddRegion("A305","A305",&SR1l_A_450lessMETless600_0_5);
    AddRegion("A405","A405",&SR1l_A_600lessMETlessInf_0_5);
    AddRegion("B105","B105",&SR1l_B_250lessMETless450_0_5);
    AddRegion("B205","B205",&SR1l_B_450lessMETless600_0_5);
    AddRegion("B305","B305",&SR1l_B_600lessMETlessInf_0_5);
    AddRegion("C105","C105",&SR1l_C_250lessMETless350_0_5);
    AddRegion("C205","C205",&SR1l_C_350lessMETless450_0_5);
    AddRegion("C305","C305",&SR1l_C_450lessMETless550_0_5);
    AddRegion("C405","C405",&SR1l_C_550lessMETless650_0_5);
    AddRegion("C505","C505",&SR1l_C_650lessMETlessInf_0_5);
    AddRegion("D105","D105",&SR1l_D_250lessMETless350_0_5);
    AddRegion("D205","D205",&SR1l_D_350lessMETless450_0_5);
    AddRegion("D305","D305",&SR1l_D_450lessMETless550_0_5);
    AddRegion("D405","D405",&SR1l_D_550lessMETlessInf_0_5);
    AddRegion("E105","E105",&SR1l_E_250lessMETless350_0_5);
    AddRegion("E205","E205",&SR1l_E_350lessMETless550_0_5);
    AddRegion("E305","E305",&SR1l_E_550lessMETlessInf_0_5);
    AddRegion("F105","F105",&SR1l_F_250lessMETless450_0_5);
    AddRegion("F205","F205",&SR1l_F_450lessMETlessInf_0_5);
    AddRegion("G105","G105",&SR1l_G_250lessMETless350_0_5);
    AddRegion("G205","G205",&SR1l_G_350lessMETless450_0_5);
    AddRegion("G305","G305",&SR1l_G_450lessMETless600_0_5);
    AddRegion("G405","G405",&SR1l_G_600lessMETlessInf_0_5);
    AddRegion("H105","H105",&SR1l_H_250lessMETless450_0_5);
    AddRegion("H205","H205",&SR1l_H_450lessMETlessInf_0_5);
   



 
    //new regions with no cut on ntight 
    AddRegion("B1NoNTight","B1NoNTight",&SR1l_B_250lessMETless450_NoNTight);
    AddRegion("B2NoNTight","B2NoNTight",&SR1l_B_450lessMETless600_NoNTight);
    AddRegion("B3NoNTight","B3NoNTight",&SR1l_B_600lessMETlessInf_NoNTight);
    AddRegion("D1NoNTight","D1NoNTight",&SR1l_D_250lessMETless350_NoNTight);
    AddRegion("D2NoNTight","D2NoNTight",&SR1l_D_350lessMETless450_NoNTight);
    AddRegion("D3NoNTight","D3NoNTight",&SR1l_D_450lessMETless550_NoNTight);
    AddRegion("D4NoNTight","D4NoNTight",&SR1l_D_550lessMETlessInf_NoNTight);
    AddRegion("F1NoNTight","F1NoNTight",&SR1l_F_250lessMETless450_NoNTight);
    AddRegion("F2NoNTight","F2NoNTight",&SR1l_F_450lessMETlessInf_NoNTight);
    AddRegion("H1NoNTight","H1NoNTight",&SR1l_H_250lessMETless450_NoNTight);
    AddRegion("H2NoNTight","H2NoNTight",&SR1l_H_450lessMETlessInf_NoNTight);

    
     






    //lower cut in Delta Phi ->0.5 and no cut in Ntight
    AddRegion("B1NoNTight05","B1NoNTight05",&SR1l_B_250lessMETless450_0_5_NoNTight);
    AddRegion("B2NoNTight05","B2NoNTight05",&SR1l_B_450lessMETless600_0_5_NoNTight);
    AddRegion("B3NoNTight05","B3NoNTight05",&SR1l_B_600lessMETlessInf_0_5_NoNTight);
    AddRegion("D1NoNTight05","D1NoNTight05",&SR1l_D_250lessMETless350_0_5_NoNTight);
    AddRegion("D2NoNTight05","D2NoNTight05",&SR1l_D_350lessMETless450_0_5_NoNTight);
    AddRegion("D3NoNTight05","D3NoNTight05",&SR1l_D_450lessMETless550_0_5_NoNTight);
    AddRegion("D4NoNTight05","D4NoNTight05",&SR1l_D_550lessMETlessInf_0_5_NoNTight);
    AddRegion("F1NoNTight05","F1NoNTight05",&SR1l_F_250lessMETless450_0_5_NoNTight);
    AddRegion("F2NoNTight05","F2NoNTight05",&SR1l_F_450lessMETlessInf_0_5_NoNTight);
    AddRegion("H1NoNTight05","H1NoNTight05",&SR1l_H_250lessMETless450_0_5_NoNTight);
    AddRegion("H2NoNTight05","H2NoNTight05",&SR1l_H_450lessMETlessInf_0_5_NoNTight);

    */

   
    
     //big region here the region are just added to form a big one, i.e. no cut on the MET that make a separation
    AddRegion("A","A",&SR1l_A);
    AddRegion("B","B",&SR1l_B);
    AddRegion("C","C",&SR1l_C);
    AddRegion("D","D",&SR1l_D);
    AddRegion("E","E",&SR1l_E);
    AddRegion("F","F",&SR1l_F);
    AddRegion("G","G",&SR1l_G);
    AddRegion("H","H",&SR1l_H);
    AddRegion("I","I",&SR1l_I);


    /*  AddRegion("C_mod","C_mod",&SR1l_C_mod);
    AddRegion("D_mod","D_mod",&SR1l_D_mod);
    AddRegion("E_mod","E_mod",&SR1l_E_mod);
    AddRegion("F_mod","F_mod",&SR1l_F_mod);
    AddRegion("G_mod","G_mod",&SR1l_G_mod);
    AddRegion("H_mod","H_mod",&SR1l_H_mod);
    AddRegion("I_mod","I_mod",&SR1l_I_mod);

      AddRegion("I_equal_4_no_ptCut","I_equal_4_no_ptCut",&SR1l_I_equal_4_no_ptCut);
      AddRegion("I_equal_4","I_equal_4",&SR1l_I_equal_4);*/

    
    /*

    //------------------------
    //regions for hight MET, that was for a test in order to change the binning in MET for the higher luminosity
    //-----------
    AddRegion("A3HIGHT","A3HIGHT",&SR1l_A_450lessMETless700);
    AddRegion("A4HIGHT","A4HIGHT",&SR1l_A_700lessMETlessInf);
    AddRegion("B2HIGHT","B2HIGHT",&SR1l_B_450lessMETless650);
    AddRegion("B3HIGHT","B3HIGHT",&SR1l_B_650lessMETlessInf);
    AddRegion("C4HIGHT","C4HIGHT",&SR1l_C_550lessMETless750);
    AddRegion("C5HIGHT","C5HIGHT",&SR1l_C_750lessMETlessInf);
    AddRegion("D3HIGHT","D3HIGHT",&SR1l_D_450lessMETless700);
    AddRegion("D4HIGHT","D4HIGHT",&SR1l_D_700lessMETlessInf);
    AddRegion("E2HIGHT","E2HIGHT",&SR1l_E_350lessMETless650);
    AddRegion("E3HIGHT","E3HIGHT",&SR1l_E_650lessMETlessInf);
    AddRegion("F1HIGHT","F1HIGHT",&SR1l_F_250lessMETless550);
    AddRegion("F2HIGHT","F2HIGHT",&SR1l_F_550lessMETlessInf);
    AddRegion("G3HIGHT","G3HIGHT",&SR1l_G_450lessMETless650);
    AddRegion("G4HIGHT","G4HIGHT",&SR1l_G_650lessMETlessInf);
    AddRegion("H1HIGHT","H1HIGHT",&SR1l_H_250lessMETless650);
    AddRegion("H2HIGHT","H2HIGHT",&SR1l_H_650lessMETlessInf);
    AddRegion("I3HIGHT","I3HIGHT",&SR1l_I_450lessMETless700);
    AddRegion("I4HIGHT","I4HIGHT",&SR1l_I_700lessMETlessInf);
*/


    //regions to test the compressed regions, test for new cut in the I regions
    /*  AddRegion("I_no_pt","I_no_pt", &I_no_pt);
    AddRegion("I_no_dphiMET","I_no_dphiMET", &I_no_dphiMET);
    AddRegion("I_no_jets","I_no_pjets", &I_no_jets);
    AddRegion("less_175","less_175", &SR1l_less175Mlb);
    AddRegion("more_175","more_175", &SR1l_more175Mlb);*/

    /*//add modified topness 0-10 for the regions B
    AddRegion("B1_topMod","B1_topMod",&SR1l_B_topMod_250lessMETless450);
    AddRegion("B2_topMod","B2_topMod",&SR1l_B_topMod_450lessMETless600);
    AddRegion("B3_topMod","B3_topMod",&SR1l_B_topMod_600lessMETlessInf);
      */

    //regions 4 jets I new definition
    /*  AddRegion("I14","I14",&SR1l_I4_250lessMETless350);
    AddRegion("I24","I24",&SR1l_I4_350lessMETless450);
    AddRegion("I34","I34",&SR1l_I4_450lessMETless550);
    AddRegion("I44","I44",&SR1l_I4_550lessMETlessInf);

    //regions 3 jets I 
    AddRegion("I13","I13",&SR1l_I3_250lessMETless350);
    AddRegion("I23","I23",&SR1l_I3_350lessMETless450);
    AddRegion("I33","I33",&SR1l_I3_450lessMETless550);
    AddRegion("I43","I43",&SR1l_I3_550lessMETlessInf);

    AddRegion("I_mod_top_inf0","I_mod_top_inf0",&SR1l_I_top_mod_inf0);
    AddRegion("I_mod_top_sup0","I_mod_top_sup0",&SR1l_I_top_mod_sup0);
    AddRegion("I_Mlb_inf175","I_mod_Mlb_inf175",&SR1l_I_Mlb_inf175);
    AddRegion("I_Mlb_sup175","I_mod_Mlb_sup175",&SR1l_I_Mlb_sup175);

    AddRegion("I1_test_mother","I1_test_mother",&SR1l_I_250lessMETless350_test_mother);
    AddRegion("I2_test_mother","I2_test_mother",&SR1l_I_350lessMETless450_test_mother);
    AddRegion("I3_test_mother","I3_test_mother",&SR1l_I_450lessMETless550_test_mother);
    AddRegion("I4_test_mother","I4_test_mother",&SR1l_I_550lessMETlessInf_test_mother);

    AddRegion("I1_Pt","I1_Pt",&SR1l_I_250lessMETless350_Pt);
    AddRegion("I2_Pt","I2_Pt",&SR1l_I_350lessMETless450_Pt);
    AddRegion("I3_Pt","I3_Pt",&SR1l_I_450lessMETless550_Pt);
    AddRegion("I4_Pt","I4_Pt",&SR1l_I_550lessMETlessInf_Pt);



 AddRegion("I1_3jets_Pt","I1_3jets_Pt",&SR1l_I3_250lessMETless350_Pt);
    AddRegion("I2_3jets_Pt","I2_3jets_Pt",&SR1l_I3_350lessMETless450_Pt);
    AddRegion("I3_3jets_Pt","I3_3jets_Pt",&SR1l_I3_450lessMETless550_Pt);
    AddRegion("I4_3jets_Pt","I4_3jets_Pt",&SR1l_I3_550lessMETlessInf_Pt);*/
    
    //all event are added in one plot, usefull when you want to compare what type of event can be in
    AddRegion("all_Event","all_Event",&SR1l_all_event);
    
    //test for new definition in modified topness in order to either improve or to mix it with the I region
    AddRegion("NJets3_Modtopness","NJets3_Modtopness",&NJets3_Modtopness);
    AddRegion("NJets3_ModtopnessA","NJets3_ModtopnessA",&NJets3_ModtopnessA);
    AddRegion("NJets3_ModtopnessB","NJets3_ModtopnessB",&NJets3_ModtopnessB);

    //test for new regions
    AddRegion("NJets3_Mlb","NJets3_Mlb", &NJets3_Mlb);
    AddRegion("NJets3_MET","NJets3_MET", &NJets3_MET);
    AddRegion("NJets4_MET","NJets4_MET", &NJets4_MET);


    //The N-1 plots for the I regions the variable written is the one on wich the cut is take off here I take alle the I1->I4 into only one region
    AddRegion("I_no_cut_MET","I_no_cut_MET", &I_no_cut_MET);
    AddRegion("I_no_cut_dphimet","I_no_cut_dphimet", &I_no_cut_dphimet);
    AddRegion("I_no_cut_dphilep","I_no_cut_dphilep", &I_no_cut_dphilep);
    AddRegion("I_no_cut_Ptlep","I_no_cut_Ptlep", &I_no_cut_Ptlep);
    AddRegion("I_no_cut_NJets","I_no_cut_NJets", &I_no_cut_NJets);
    AddRegion("I_no_cut_ident","I_no_cut_ident", &I_no_cut_ident);
    

    //try for a new definition for the I region with less jets
    AddRegion("I_more_4","I_more_4", &SR1l_I_more_4);
    AddRegion("I_more_3","I_more_3", &SR1l_I_more_3);

    AddRegion("I_equal_4","I_equal_4", &SR1l_I_equal_4);
    AddRegion("I_equal_3","I_equal_3", &SR1l_I_equal_3);
    
    AddRegion("I_equal_34","I_equal_34", &I_equal_34);

    AddRegion("I_equal_34_region1","I_equal_34_region1", &SR1l_I_equal_34_250lessMETless350);
    AddRegion("I_equal_34_region2","I_equal_34_region2", &SR1l_I_equal_34_350lessMETless450);
    AddRegion("I_equal_34_region3","I_equal_34_region3", &SR1l_I_equal_34_450lessMETless550);
    AddRegion("I_equal_34_region4","I_equal_34_region4", &SR1l_I_equal_34_550lessMETlessInf);



    AddRegion("I_equal_3_region1","I_equal_3_region1", &SR1l_I_equal3_250lessMETless350);
    AddRegion("I_equal_3_region2","I_equal_3_region2", &SR1l_I_equal3_350lessMETless450);
    AddRegion("I_equal_3_region3","I_equal_3_region3", &SR1l_I_equal3_450lessMETless550);
    AddRegion("I_equal_3_region4","I_equal_3_region4", &SR1l_I_equal3_550lessMETlessInf);



    AddRegion("I_equal_4_region1","I_equal_4_region1", &SR1l_I_equal4_250lessMETless350);
    AddRegion("I_equal_4_region2","I_equal_4_region2", &SR1l_I_equal4_350lessMETless450);
    AddRegion("I_equal_4_region3","I_equal_4_region3", &SR1l_I_equal4_450lessMETless550);
    AddRegion("I_equal_4_region4","I_equal_4_region4", &SR1l_I_equal4_550lessMETlessInf);




    //try of new definition of the I regions go into the definition to change them or to make them better 
    AddRegion("I_equal_4_cut_Pt_300_region1","I_equal_4_cut_Pt_300_region1", &SR1l_I_equal4_cut_Pt_300_250lessMETless350);
    AddRegion("I_equal_4_cut_Pt_300_region2","I_equal_4_cut_Pt_300_region2", &SR1l_I_equal4_cut_Pt_300_350lessMETless450);
    AddRegion("I_equal_4_cut_Pt_300_region3","I_equal_4_cut_Pt_300_region3", &SR1l_I_equal4_cut_Pt_300_450lessMETless550);
    AddRegion("I_equal_4_cut_Pt_300_region4","I_equal_4_cut_Pt_300_region4", &SR1l_I_equal4_cut_Pt_300_550lessMETlessInf);


    AddRegion("I_equal_3_cut_Pt_300_region1","I_equal_3_cut_Pt_300_region1", &SR1l_I_equal3_cut_Pt_300_250lessMETless350);
    AddRegion("I_equal_3_cut_Pt_300_region2","I_equal_3_cut_Pt_300_region2", &SR1l_I_equal3_cut_Pt_300_350lessMETless450);
    AddRegion("I_equal_3_cut_Pt_300_region3","I_equal_3_cut_Pt_300_region3", &SR1l_I_equal3_cut_Pt_300_450lessMETless550);
    AddRegion("I_equal_3_cut_Pt_300_region4","I_equal_3_cut_Pt_300_region4", &SR1l_I_equal3_cut_Pt_300_550lessMETlessInf);

    AddRegion("I_equal_34_cut_Pt_300_region1","I_equal_34_cut_Pt_300_region1", &I_equal_34_cut_Pt_300_250_350);
    AddRegion("I_equal_34_cut_Pt_300_region2","I_equal_34_cut_Pt_300_region2", &I_equal_34_cut_Pt_300_350_450);
    AddRegion("I_equal_34_cut_Pt_300_region3","I_equal_34_cut_Pt_300_region3", &I_equal_34_cut_Pt_300_450_550);
    AddRegion("I_equal_34_cut_Pt_300_region4","I_equal_34_cut_Pt_300_region4", &I_equal_34_cut_Pt_300_550_Inf);


    
    AddRegion("I_equal_4_cut_Pt_300","I_equal_4_cut_Pt_300", &SR1l_I_equal_4_cut_Pt_300);
    AddRegion("I_equal_3_cut_Pt_300","I_equal_3_cut_Pt_300", &SR1l_I_equal_3_cut_Pt_300);
    
    AddRegion("I_equal_34_cut_Pt_300","I_equal_34_cut_Pt_300", &I_equal_34_cut_Pt_300); 
    
 


    AddRegion("I_equal_34_cut_Pt_300_region1_looseB","I_equal_34_cut_Pt_300_region1_looseB", &I_equal_34_cut_Pt_300_250_350_looseB);
    AddRegion("I_equal_34_cut_Pt_300_region2_looseB","I_equal_34_cut_Pt_300_region2_looseB", &I_equal_34_cut_Pt_300_350_450_looseB);
    AddRegion("I_equal_34_cut_Pt_300_region3_looseB","I_equal_34_cut_Pt_300_region3_looseB", &I_equal_34_cut_Pt_300_450_550_looseB);
    AddRegion("I_equal_34_cut_Pt_300_region4_looseB","I_equal_34_cut_Pt_300_region4_looseB", &I_equal_34_cut_Pt_300_550_Inf_looseB);

    AddRegion("I_equal_34_cut_Pt_300_looseB","I_equal_34_cut_Pt_300_looseB", &I_equal_34_cut_Pt_300_looseB); 





    AddRegion("I_equal_34_cut_Pt_300_modtop_Sup0","I_equal_34_cut_Pt_300_modtop_Sup0", &I_equal_34_cut_Pt_300_modtop_Sup0); 
   
    AddRegion("I_equal_34_cut_Pt_300_region1_modtop_Sup0","I_equal_34_cut_Pt_300_region1_modtop_Sup0", &I_equal_34_cut_Pt_300_250_350_modtop_Sup0);
    AddRegion("I_equal_34_cut_Pt_300_region2_modtop_Sup0","I_equal_34_cut_Pt_300_region2_modtop_Sup0", &I_equal_34_cut_Pt_300_350_450_modtop_Sup0);
    AddRegion("I_equal_34_cut_Pt_300_region3_modtop_Sup0","I_equal_34_cut_Pt_300_region3_modtop_Sup0", &I_equal_34_cut_Pt_300_450_550_modtop_Sup0);
    AddRegion("I_equal_34_cut_Pt_300_region4_modtop_Sup0","I_equal_34_cut_Pt_300_region4_modtop_Sup0", &I_equal_34_cut_Pt_300_550_Inf_modtop_Sup0);



    AddRegion("I_equal_34_cut_Pt_300_modtop_Inf0","I_equal_34_cut_Pt_300_modtop_Inf0", &I_equal_34_cut_Pt_300_modtop_Inf0); 
   
    AddRegion("I_equal_34_cut_Pt_300_region1_modtop_Inf0","I_equal_34_cut_Pt_300_region1_modtop_Inf0", &I_equal_34_cut_Pt_300_250_350_modtop_Inf0);
    AddRegion("I_equal_34_cut_Pt_300_region2_modtop_Inf0","I_equal_34_cut_Pt_300_region2_modtop_Inf0", &I_equal_34_cut_Pt_300_350_450_modtop_Inf0);
    AddRegion("I_equal_34_cut_Pt_300_region3_modtop_Inf0","I_equal_34_cut_Pt_300_region3_modtop_Inf0", &I_equal_34_cut_Pt_300_450_550_modtop_Inf0);
    AddRegion("I_equal_34_cut_Pt_300_region4_modtop_Inf0","I_equal_34_cut_Pt_300_region4_modtop_Inf0", &I_equal_34_cut_Pt_300_550_Inf_modtop_Inf0);
    
    
    
    AddRegion("I1_modtop_Sup0","I1_modtop_Sup0", &SR1l_I_250lessMETless350_topmod_Sup0);
    AddRegion("I2_modtop_Sup0","I2_modtop_Sup0", &SR1l_I_350lessMETless450_topmod_Sup0);
    AddRegion("I3_modtop_Sup0","I3_modtop_Sup0", &SR1l_I_450lessMETless550_topmod_Sup0);
    AddRegion("I4_modtop_Sup0","I4_modtop_Sup0", &SR1l_I_550lessMETlessInf_topmod_Sup0);


    AddRegion("I1_modtop_Inf0","I1_modtop_Inf0", &SR1l_I_250lessMETless350_topmod_Inf0);
    AddRegion("I2_modtop_Inf0","I2_modtop_Inf0", &SR1l_I_350lessMETless450_topmod_Inf0);
    AddRegion("I3_modtop_Inf0","I3_modtop_Inf0", &SR1l_I_450lessMETless550_topmod_Inf0);
    AddRegion("I4_modtop_Inf0","I4_modtop_Inf0", &SR1l_I_550lessMETlessInf_topmod_Inf0);


    AddRegion("I_topmod_Inf0","I_topmod_Inf0", &I_topmod_Inf0);
    AddRegion("I_topmod_Sup0","I_topmod_Sup0", &I_topmod_Sup0);

    AddRegion("I1_nocutB","I1_nocutB",&SR1l_I_250lessMETless350_nocutB);
    AddRegion("I2_nocutB","I2_nocutB",&SR1l_I_350lessMETless450_nocutB);
    AddRegion("I3_nocutB","I3_nocutB",&SR1l_I_450lessMETless550_nocutB);
    AddRegion("I4_nocutB","I4_nocutB",&SR1l_I_550lessMETlessInf_nocutB);

    AddRegion("I34_nocutB","I34_nocutB", &SR1l_I34_nucutB);
    
    AddRegion("I_cut_tight","I_cut_tight", &I_cut_tight);
    AddRegion("I_cut_med","I_cut_med", &I_cut_med);
    AddRegion("I_cut_loose","I_cut_loose", &I_cut_loose);

    AddRegion("I1_NocutPtLep","I1_NocutPtLep", &SR1l_I_250lessMETless350_NocutPt);
    AddRegion("I2_NocutPtLep","I2_NocutPtLep", &SR1l_I_350lessMETless450_NocutPt);
    AddRegion("I3_NocutPtLep","I3_NocutPtLep", &SR1l_I_450lessMETless550_NocutPt);
    AddRegion("I4_NocutPtLep","I4_NocutPtLep", &SR1l_I_550lessMETlessInf_NocutPt);




//  fillYieldsVector();// @MJ@ TODO probably not needed when I do not need to zero negative
                                                                                                                       

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);
  


    ////old luminosity
     SetLumi(35.9);

    ///test with the new luminosity
    // SetLumi(120);
    //SetLumi(150);
   
     Add2DHisto("pt_lead_jet","MET");
    

     //just a try to have different plots such as the number of jets in finction of the stop mass, if I remember well it should work but you need to check the 2D plots
     /*  Add2DHisto("nJets","StopMass");
     Add2DHisto("nJets","NeutralinoMass");
     Add2DHisto("MET","StopMass");
     Add2DHisto("MET","NeutralinoMass");
     Add2DHisto("Mlb","StopMass");
     Add2DHisto("Mlb","NeutralinoMass");
     Add2DHisto("dphi","StopMass");
     Add2DHisto("dphi","NeutralinoMass");
     Add2DHisto("topnessMod","StopMass");
     Add2DHisto("topnessMod","NeutralinoMass");
     Add2DHisto("MT","StopMass");
     Add2DHisto("MT","NeutralinoMass");
     Add2DHisto("lep1_pt","StopMass");
     Add2DHisto("lep1_pt","NeutralinoMass");
     Add2DHisto("lep1_pt","dphi");*/
     

     //try for 3D plots but it don't work
     /*  Add3DHisto("NeutralinoMass","StopMass","nJets");
	 Add3DHisto("NeutralinoMass","StopMass","nBJets");
	 Add3DHisto("NeutralinoMass","StopMass","MET");
	 Add3DHisto("NeutralinoMass","StopMass","MT2W");
²	 Add3DHisto("NeutralinoMass","StopMass","MT");
	 Add3DHisto("NeutralinoMass","StopMass","topnessMod");
	 Add3DHisto("NeutralinoMass","StopMass","dphi");
	 Add3DHisto("NeutralinoMass","StopMass","Mlb");
	 Add3DHisto("NeutralinoMass","StopMass","lep1_pt");*/
     Create1DHistos();
     Add2DHisto("StopMass","NeutralinoMass");
     WriteXMLConfig(); 
     


    
}



void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{


    counter++;
    nthentry++;
  

    checkNegativeYields = false;
    if(nthentry == myEvent.nentries)
    {
      //cout << "checking histogram for negative values" << endl;
        nthentry =0;
        //checkNegativeYields = true; //@MJ@ TODO be aware of this, I can use multiprocessing now but I can endup with incorrect MC yields!!!
    }


    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);
 
 
 
    //weight on signal
    TFile *file = NULL;
    double weightSignal = -13;
    double weightSignal_Tight = -13;
     if(currentProcessType == "signal")     
      {
	//cout<<"j'ai passé le premier test"<<endl;
	if(currentDataset != storedDataset)// && h2 == NULL) //@MJ@ TODO this can work only with one signal dataset!!!
	  {
	    storedDataset = currentDataset;
	    TString fName =  babyTuplePath+currentDataset+".root";
	    file = new TFile(fName);
	    h2 = (TH2D*)file->Get("histNEvts")->Clone();
	    h3 = (TH3D*)file->Get("h_counterSMS")->Clone();
	  }
	if (h2 == NULL) throw std::runtime_error("The histogram used for CS was not filled!");

	//	for the mass of the stop we want in order to have only one or more couple of mass to study
		if (  (myEvent.mass_stop == 350 && myEvent.mass_lsp == 100))
	  /* || (myEvent.mass_stop == 500 && myEvent.mass_lsp == 325)
	      || (myEvent.mass_stop == 500 && myEvent.mass_lsp == 275)
	      || (myEvent.mass_stop == 525 && myEvent.mass_lsp == 325)
	      || (myEvent.mass_stop == 525 && myEvent.mass_lsp == 300)
	      || (myEvent.mass_stop == 525 && myEvent.mass_lsp == 275)
	      || (myEvent.mass_stop == 475 && myEvent.mass_lsp == 300)
	      || (myEvent.mass_stop == 475 && myEvent.mass_lsp == 325)
	      || (myEvent.mass_stop == 475 && myEvent.mass_lsp == 275))*/
	  /* ||
	     (myEvent.mass_stop == 300 && myEvent.mass_lsp == 100) ||
	     (myEvent.mass_stop == 500 && myEvent.mass_lsp == 100))
	     
								    */
	
	 {
	//  cout<<"et le deuxieme aussi 1!!!!"<<endl;
	    
	    int totalNrOfEvents = (int)h2->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp));
	    
	    float xsec = stopXSEC->GetBinContent(stopXSEC->FindBin(myEvent.mass_stop));
	    weightSignal =  1000* xsec * GetLumi() /totalNrOfEvents; //@MJ@ TODO improve my method
	    weightSignal *=  myEvent.weight_lepSF*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,27))  ) ;
	    //	  weightSignal *=  myEvent.weight_vetoLepSF*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,30))  ) ;
	    weightSignal_Tight = weightSignal;
	    weightSignal *=  myEvent.weight_btagsf*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,14))  ) ;
	    weightSignal_Tight *=  myEvent.weight_tightbtagsf*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,37))  ) ;
	    // AddToHisto(myEvent.mass_stop, myEvent.mass_lsp, weightSignal*myEvent.ngoodjets); 
	       }
	
	
	        	
	       else	
	  {
	    weightSignal=0;
	    weightSignal_Tight=0;
	    }
	   



 













	  /*myEvent.pfmet
  myEvent.MT2W
    myEvent.mt_met_lep
    myEvent.ngoodjets
myEvent.ngoodbtags
myEvent.topnessMod
myEvent.dphi_ak4pfjets_met
myEvent.Mlb
myEvent.lep1_pt	
ak4pfjets_pt
myEvent.ak4pfjets_leadbtag_p4.Pt()*/ 
		//here it is to plot the efficiency for every couple of mass you have to put the event with the cut on the moyenne variable and in the moyenne 2 wariable you put the one without cut see the ../sonicScrewdriver/interface/BabyScrewdriver.h for the plots
	    if ( currentProcessClass=="2Chargino" && (myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2) == true && myEvent.ak4pfjets_passMEDbtag->at(0) == false)//&&  SR1l() && myEvent.pfmet>=250 )//()
	    {
	      // cout<<"ccccx"<<endl;
	      /*  if ()
		{
		  moyenne->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal); 	     
		  }*/
	      // moyenne->Fill(myEvent.ak4pfjets_pt,myEvent.pfmet,weightSignal);
	      //somme +=weightSignal;
	      //     if( myEvent.lep1_mc_motherid==1000024 || myEvent.lep1_mc_motherid==(-1000024))//myEvent.lep1_mc_motherid==24 || myEvent.lep1_mc_motherid==(-24)) // // ||
	      //{
	      //  moyenne->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal );
	      moyenne->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal*(myEvent.ak4pfjets_pt-myEvent.ak4pfjets_pt2) ); 
		  //  cout<<weightSignal<<endl;
		  //}
	      //cout<<myEvent.ak4pfjets_pt<<endl;
	      // cout<<myEvent.ak4pfjets_leadbtag_p4.Pt()<<endl;
	      //cout<<myEvent.lep1_pt<<endl;
	    moyenne2->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal);
	    }
	     
      //	moyenne = &moyenne3;//	moyenne2->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal);
	

	//	x++;
	// cout<<"et j'aditionne le x "<<x<<endl;*/
      
	    //same but it could be used only for one point also
if ( currentProcessClass=="NoChargino" && (myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2) == true && myEvent.ak4pfjets_passMEDbtag->at(0) == false)//&&  SR1l() && myEvent.pfmet>=250 )//(myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2) == true && myEvent.ak4pfjets_passMEDbtag->at(0) == false)
      {
	// cout<<"ccccx"<<endl;
	      /*  if ()
		{
		  moyenne->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal); 	     
		  }*/
	      // moyenne->Fill(myEvent.ak4pfjets_pt,myEvent.pfmet,weightSignal);
	      //somme +=weightSignal;
	      //     if( myEvent.lep1_mc_motherid==1000024 || myEvent.lep1_mc_motherid==(-1000024))//myEvent.lep1_mc_motherid==24 || myEvent.lep1_mc_motherid==(-24)) // // ||
	      //{
		  moyenne3->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal*(myEvent.ak4pfjets_pt-myEvent.ak4pfjets_pt2) ); 
		  //  cout<<weightSignal<<endl;
		  //}
	      //cout<<myEvent.ak4pfjets_pt<<endl;
	      // cout<<myEvent.ak4pfjets_leadbtag_p4.Pt()<<endl;
	      //cout<<myEvent.lep1_pt<<endl;
	    moyenne4->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal);
	    }
	     
      //	moyenne = &moyenne3;//	moyenne2->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal);
	

	//	x++;
	// cout<<"et j'aditionne le x "<<x<<endl;*/
      


	    if ( currentProcessClass=="mixt" && (myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2) == true && myEvent.ak4pfjets_passMEDbtag->at(0) == false)// &&SR1l() && myEvent.pfmet>=250 )//(myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2) == true && myEvent.ak4pfjets_passMEDbtag->at(0) == false)
	    {
	      // cout<<"ccccx"<<endl;
	      /*  if ()
		{
		  moyenne->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal); 	     
		  }*/
	      // moyenne->Fill(myEvent.ak4pfjets_pt,myEvent.pfmet,weightSignal);
	      //somme +=weightSignal;
	      //     if( myEvent.lep1_mc_motherid==1000024 || myEvent.lep1_mc_motherid==(-1000024))//myEvent.lep1_mc_motherid==24 || myEvent.lep1_mc_motherid==(-24)) // // ||
	      //{
		  moyenne5->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal*(myEvent.ak4pfjets_pt-myEvent.ak4pfjets_pt2 ) ); 
		  // cout<<weightSignal<<endl;
		  //}
	      //cout<<myEvent.ak4pfjets_pt<<endl;
	      // cout<<myEvent.ak4pfjets_leadbtag_p4.Pt()<<endl;
		  moyenne6->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal);
	      //cout<<myEvent.lep1_pt<<endl;
	    }
	     
      //	moyenne = &moyenne3;//	moyenne2->Fill(myEvent.mass_stop, myEvent.mass_lsp, weightSignal);
	

	//	x++;
	// cout<<"et j'aditionne le x "<<x<<endl;*/
      }
    
    
     








  // test 
    if( (currentProcessClass == "bkgLostLepton") || (currentProcessClass == "bkgOneLepFromTop")) {
      if (myEvent.is1lepFromTop) currentProcessClass = "bkgOneLepFromTop";
      else currentProcessClass = "bkgLostLepton";
    }
    
    
    if(currentDataset != storedDataset && currentProcessType == "background") //@MJ@ TODO this can work only with one signal dataset!!!
      {
        storedDataset = currentDataset;
        scale1fbS2 = 1;
	
        if(currentDataset == "ttbar_singleLeptFromTbar_madgraph_pythia8_25ns")
	  {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns.root"; 
	    scale1fbS2= 1.13617e+07/(1.13617e+07 + 4.63189e+07);
	    
	  }
        else if(currentDataset == "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns")
	  {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromTbar_madgraph_pythia8_25ns.root"; 
            scale1fbS2= 4.63189e+07/(1.13617e+07 + 4.63189e+07);
	    
	  }
        else if(currentDataset == "ttbar_singleLeptFromT_madgraph_pythia8_25ns")
	  {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns.root"; 
            scale1fbS2=1.16509e+07/(1.16509e+07 + 4.08199e+07);
	    
	  }
        else if(currentDataset == "ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns")
	  {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromT_madgraph_pythia8_25ns.root"; 
            scale1fbS2=4.08199e+07/(1.16509e+07 + 4.08199e+07);
	  }
        else if(currentDataset == "ttbar_diLept_madgraph_pythia8_25ns")
	  {
            TString fBkgName =  babyTuplePath+"ttbar_diLept_madgraph_pythia8_ext1_25ns.root"; 
            scale1fbS2=5.77109e+06/(5.77109e+06 + 2.34556e+07);
	  }
        else if(currentDataset == "ttbar_diLept_madgraph_pythia8_ext1_25ns")
	  {
            TString fBkgName =  babyTuplePath+"ttbar_diLept_madgraph_pythia8_25ns.root"; 
            scale1fbS2=2.34556e+07/(5.77109e+06 + 2.34556e+07);
	  }
        else
	  {
            scale1fbS2 = 1; 
	  }
      }
    
    
    double weight = 0; 
    




    //here it was some test that I made about the origin of the lepton in the mixt scenario (T2tb)
    // if( !myEvent.is_data && counter<100 ) {
    
    
    /*  ofstream fichier1("mother.dat", ios::out | ios::app);
	if(currentProcessType == "signal")
	{ //   cout << "counter " << counter << endl;
	
	// match leading lepton first
	     int genLep_matchedTo_selLep__idx = -1;
	double matched_lep_dr = 0.1;
	for(int iGen=0; iGen<(int) myEvent.genleps_p4.size(); iGen++) {
	if( abs(myEvent.genleps_id[iGen]) != abs(myEvent.lep1_pdgid) ) continue;
	if( !myEvent.genleps_isLastCopy[iGen] ) continue;
 	if( !myEvent.genleps_fromHardProcessFinalState[iGen] &&
	!myEvent.genleps_fromHardProcessDecayed[iGen]       ) continue;
	if( ROOT::Math::VectorUtil::DeltaR(myEvent.genleps_p4[iGen], myEvent.lep1_p4) < matched_lep_dr ) {
	genLep_matchedTo_selLep__idx  = iGen;
	break;
	}
	}
    */
    
    /*      if( myEvent.lep1_mc_motherid!=24 && myEvent.lep1_mc_motherid!=(-24) && myEvent.lep1_mc_motherid!=1000024 && myEvent.lep1_mc_motherid!=(-1000024))
	    {
	    fichier1<< "lepton mother1 "<< myEvent.lep1_mc_motherid;//<<" "<<myEvent.genleps_id[1]<<" "<<myEvent.genleps_id[2] ;//<<" from "<<myEvent.lep1_mc_gmotherid;
            // If matched selected lepton, 
	    // if( genLep_matchedTo_selLep__idx>0 ) {
	    // cout << " mother2 " << myEvent.genleps_motherid[genLep_matchedTo_selLep__idx] << " from    " << myEvent. genleps_gmotherid[genLep_matchedTo_selLep__idx] ;  
		  //	  }
            fichier1 << endl;
	}
    }
 
  fichier1.close();*/



    if(currentProcessType == "signal")
    {
      weight = weightSignal;
      // x++;
      //cout<<x<<endl;
      //weight=1.0;
     
  
    }
    else if (currentProcessType == "data") 
    {
        weight = 1.0;
    }
    else 
      { 
	weight = getWeight(currentProcessType, GetLumi(), scale1fbS2,true); 
	if (myEvent.ak4pfjets_passMEDbtag->size()>0)
	  {
	    if ( (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 )== true)
	      { 
		if ( SR1l() && myEvent.pfmet>=250)
		  {
		    ratio+=weight;
		    // cout<<ratio<<endl;
		  }
		ratio2+=weight;
	      }
	  }
      }
	


// Check WNJets genPt
   if( currentDataset == "W1JetsToLNu_madgraph_pythia8_25ns" || currentDataset == "W2JetsToLNu_madgraph_pythia8_25ns" 
     || currentDataset == "W3JetsToLNu_madgraph_pythia8_25ns" || currentDataset == "W4JetsToLNu_madgraph_pythia8_25ns") {
     if( myEvent.nupt>200.0 ) { weight = 0.;}
   }
   if (currentDataset == "W1JetsToLNu_nupT200_madgraph_pythia8_25ns" || currentDataset == "W2JetsToLNu_nupT200_madgraph_pythia8_25ns"
     || currentDataset == "W3JetsToLNu_nupT200_madgraph_pythia8_25ns" || currentDataset == "W4JetsToLNu_nupT200_madgraph_pythia8_25ns") {
     if (myEvent.nupt<200.0) { weight = 0.; }
   }



   //test in order the have multiple mass on the same plot in the schedule plots
   // if  (myEvent.mass_stop == 800 && myEvent.mass_lsp == 600) 
   AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields);
   //  if (currentProcessType == "background")   AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields);

    
   // if (myEvent.mass_stop == 800 && myEvent.mass_lsp == 600 && (currentProcessClass=="2Chargino_600_400")) AutoFillProcessClass("2Chargino_800_600",weight,checkNegativeYields);
    // if (myEvent.mass_stop == 600 && myEvent.mass_lsp == 400 && (currentProcessClass=="2Chargino")) AutoFillProcessClass("2chargino_600_400",weight,checkNegativeYields);
   // if (myEvent.mass_stop == 300 && myEvent.mass_lsp == 100 && (currentProcessClass=="2Chargino_600_400")) AutoFillProcessClass("2Chargino_300_100",weight,checkNegativeYields); 
   // if (myEvent.mass_stop == 800 && myEvent.mass_lsp == 600 && (currentProcessClass=="NoChargino_600_400")) AutoFillProcessClass("NoChargino_800_600",weight,checkNegativeYields);
   //if (myEvent.mass_stop == 600 && myEvent.mass_lsp == 400 && (currentProcessClass=="NoChargino")) AutoFillProcessClass("Nochargino_600_400",weight,checkNegativeYields);
   //  if (myEvent.mass_stop == 300 && myEvent.mass_lsp == 100 && (currentProcessClass=="NoChargino_600_400")) AutoFillProcessClass("NoChargino_300_100",weight,checkNegativeYields);
   // if (myEvent.mass_stop == 800 && myEvent.mass_lsp == 600 && (currentProcessClass=="mixt_600_400")) AutoFillProcessClass("mixt_800_600",weight,checkNegativeYields);
   // if (myEvent.mass_stop == 600 && myEvent.mass_lsp == 400 && (currentProcessClass=="mixt")) AutoFillProcessClass("mixt_600_400",weight,checkNegativeYields);
   //  if (myEvent.mass_stop == 300 && myEvent.mass_lsp == 100 && (currentProcessClass=="mixt_600_400")) AutoFillProcessClass("mixt_300_100",weight,checkNegativeYields);   
   /*  if(counter % 10000 == 0)
       {
	//cout << counter << endl;
	}*/
    
}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
  // ######################
  //  Plot configuration and production
  // ######################
  
  // Schedule plots
  //

  // cout<<"exemple"<<endl;
  SchedulePlots("1DSuperimposed");
  //cout<<"exemple"<<endl;
  SchedulePlots("1DStack");
  SchedulePlots("1DDataMCComparison");
  SchedulePlots("2D");
  // SchedulePlots("2DFrom3DProjection");
  
  // Config plots
  
  SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
  SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 13 TeV");
  
  SetGlobalBoolOption("Plot", "exportPdf", false);
  SetGlobalBoolOption("Plot", "exportEps", false);
  SetGlobalBoolOption("Plot", "exportPng", false);
  
  // Make and write the plots
   
         cout << endl;
  cout << "   > Making plots..." << endl; 
   MakePlots();
  cout << "   > Saving plots..." << endl;
  WritePlots("./plotsTest/");
  
  
  // ######################
  //  Tables and other stuff uncomment the totyield and the CombineCardMaker in order to write the table that will be used either with the combine or to plot the number of event per region
  // ######################
  
  
  /*    vector<string> totYield2 = {"H2"}; 
	vector<string> totYield1 = {"H1"}; 
	
	
	vector<string> totYield4 = {"B2"}; 
	vector<string> totYield3 = {"E1"}; 
  */
  
  // vector<string> totYield = {"baseline", "test"};I1_Pt
    // vector<string> totYield = {"I1","I2","I3","I4"};//,"I1_3jets_Pt","I2_3jets_Pt","I3_3jets_Pt","I4_3jets_Pt"};
    // vector<string> totYield1 = {"I_equal_34_region1","I_equal_34_region2","I_equal_34_region3","I_equal_34_region4"};//vector<string> totYield = {"I1","I2","I3","I4","I1_test_mother","I2_test_mother","I3_test_mother","I4_test_mother"};

    /*  CombineCardMaker card5(this, totYield , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card5.Print("card_I_mixt.tab", 4);
  card5.ProduceCard("card_I_mixt.tab" ,"./cards14");
  CombineCardMaker card6(this, totYield1 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card6.Print("card_I34_mixt.tab", 4);
  card6.ProduceCard("card_I34_mixt.tab" ,"./cards15");*/
   	







  // vector<string> totYield = {"A","B","C","D","E","F","G","H","I","C_mod","D_mod","E_mod","F_mod","G_mod","H_mod","I_mod"};
  //  vector<string> totYield = {"I","I_mod_top_inf0","I_mod_top_sup0","I_Mlb_inf175","I_Mlb_sup175"};
    //--------------------------------
    //study for compressed regions 
    //----------------------------------------------
  
  // vector<string> totYield3 = {"A1","A2","A3","A4","B1","B2","B3","C1","C2","C3","C4","C5","D1","D2","D3","D4","E1","E2","E3","F1","F2","G1","G2","G3","G4","H1","H2","I1","I2","I3","I4"};
      /*  CombineCardMaker card(this, totYield , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
	  card.Print("card_test_selection_2Chargino.tab", 4);*/
  //    CombineCardMaker card3(this, totYield , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  // card3.Print("card_test_Pt_I_3jets_mixt.tab", 4);
 /* CombineCardMaker card4(this, totYield , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
 card4.Print("card_test_selection_NoChargino.tab", 4);*/
  //,"I14","I24","I34","I44","I13","I23","I33","I43"};

	 /*  CombineCardMaker card2(this, totYield1 , "lepChannel", "2Chargino_800_600", "StopMass", "NeutralinoMass");
	     card2.Print("card_pour_oral2.tab", 4);*/
	 
  /*    vector<string> totYield5 = {"I108","I208","I308","I408"};
	
	vector<string> totYield6 = {"A105","A205","A305","A405","B105","B205","B305","C105","C205","C305","C405","C505","D105","D205","D305","D405","E105","E205","E305","F105","F205","G105","G205","G305","G405","H105","H205"};
	
	vector<string> totYield7= {"A1","A2","A3","A4","B1NoNTight","B2NoNTight","B3NoNTight","C1","C2","C3","C4","C5","D1NoNTight","D2NoNTight","D3NoNTight","D4NoNTight","E1","E2","E3","F1NoNTight","F2NoNTight","G1","G2","G3","G4","H1NoNTight","H2NoNTight"};
	
	vector<string> totYield8 = {"A105","A205","A305","A405","B1NoNTight05","B2NoNTight05","B3NoNTight05","C105","C205","C305","C405","C505","D1NoNTight05","D2NoNTight05","D3NoNTight05","D4NoNTight05","E105","E205","E305","F1NoNTight05","F2NoNTight05","G105","G205","G305","G405","H1NoNTight05","H2NoNTight05"};*/
  
  
  
  
  
  /*  vector<string> totYield3 = {"A1","A2","A3HIGHT","A4HIGHT","B1","B2HIGHT","B3HIGHT","C1","C2","C3","C4HIGHT","C5HIGHT","D1","D2","D3HIGHT","D4HIGHT","E1","E2HIGHT","E3HIGHT","F1HIGHT","F2HIGHT","G1","G2","G3HIGHT","G4HIGHT","H1HIGHT","H2HIGHT","I1","I2","I3HIGHT","I4HIGHT","B1_topMod","B2_topMod","B3_topMod"};
  // CombineCardMaker card(this, totYield2 , "lepChannel", "2chargino", "StopMass", "NeutralinoMass");
  // card.Print("card.tab", 4);
  /* vector<string> totYield = {"A1A","A2A","A3A","A4A","B1A","B2A","B3A","C1A","C2A","C3A","C4A","C5A","D1A","D2A","D3A","D4A","E1A","E2A","E3A","F1A","F2A","G1A","G2A","G3A","G4A","H1A","H2A","I1A","I2A","I3A","I4A"};
     vector<string> totYield1 = {"A1B","A2B","A3B","A4B","B1B","B2B","B3B","C1B","C2B","C3B","C4B","C5B","D1B","D2B","D3B","D4B","E1B","E2B","E3B","F1B","F2B","G1B","G2B","G3B","G4B","H1B","H2B","I1B","I2B","I3B","I4B"};
     vector<string> totYield3 = {"A1C","A2C","A3C","A4C","B1C","B2C","B3C","C1C","C2C","C3C","C4C","C5C","D1C","D2C","D3C","D4C","E1C","E2C","E3C","F1C","F2C","G1C","G2C","G3C","G4C","H1C","H2C","I1C","I2C","I3C","I4C"};
     //vector<string> totYield4 = {"A1","A2","A3","A4","B1","B2","B3","C1","C2","C3","C4","C5","D1","D2","D3","D4","E1","E2","E3","F1","F2","G1","G2","G3","G4","H9","H10","I1","I2","I3","I4"};*/
  //vector<string> totYield = {"SR1l_C_250lessMETless350" , "SR1l_C_350lessMETless450" , "SR1l_C_450lessMETless550" , "SR1l_C_550lessMETless650" , "SR1l_C_650lessMETlessInf" };
  /* 
  TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+ "table.tab", 4);
     TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+ "sbottomSignal.tex", 4);
  */
  
   //outfile in order to use it with resultPlot (here for 2Chargino) 
  /*  CombineCardMaker card(this, totYield , "lepChannel", "2Chargino_600_400", "StopMass", "NeutralinoMass");
  card.Print("card_donnee.tab", 4);
  CombineCardMaker card1(this, totYield , "lepChannel", "NoChargino_600_400", "StopMass", "NeutralinoMass");
  card1.Print("card_donnee1.tab", 4);
  CombineCardMaker card2(this, totYield , "lepChannel", "mixt_600_400", "StopMass", "NeutralinoMass");
  card2.Print("card_donnee2.tab", 4);
  CombineCardMaker card3(this, totYield , "lepChannel", "2Chargino_800_600", "StopMass", "NeutralinoMass");
  card3.Print("card_donnee3.tab", 4);
  CombineCardMaker card4(this, totYield , "lepChannel", "NoChargino_800_600", "StopMass", "NeutralinoMass");
  card4.Print("card_donnee4.tab", 4);
  CombineCardMaker card5(this, totYield , "lepChannel", "mixt_800_600", "StopMass", "NeutralinoMass");
  card5.Print("card_donnee5.tab", 4);*/
  	
    





    
  
  //  vector<string> totYield = {"A1","A2","A3","A4","B1","B2","B3","C1","C2","C3","C4","C5","D1","D2","D3","D4","E1","E2","E3","F1","F2","G1","G2","G3","G4","H1","H2","I108","I208","I308","I408"};	 
   //   vector<string> totYield1 = {"A1","A2","A3","A4","B1","B2","B3","C1CMPR","C2CMPR","C3CMPR","C4CMPR","C5CMPR","D1CMPR","D2CMPR","D3CMPR","D4CMPR","E1CMPR","E2CMPR","E3CMPR","F1CMPR","F2CMPR","G1CMPR","G2CMPR","G3CMPR","G4CMPR","H1CMPR","H2CMPR","I1","I2","I3","I4"};     	
  // vector<string> totYield2 = {"A1","A2","A3","A4","B1","B2","B3","C1","C2","C3","C4","C5","D1","D2","D3","D4","E1","E2","E3","F1","F2","G1","G2","G3","G4","H1","H2"};
  // vector<string> totYield3 = {"Baseline","A1","A2","A3","A4","B1","B2","B3","C1","C2","C3","C4","C5","D1","D2","D3","D4","E1","E2","E3","F1","F2","G1","G2","G3","G4","H1","H2","I1","I2","I3","I4"};     	
    // vector<string> totYield2 = {"I1","I2","I3","I4"};
    //  vector<string> totYield3 = {"I1_NocutPtLep","I2_NocutPtLep","I3_NocutPtLep","I4_NocutPtLep"};
  
  /* vector<string> totYield4 = {"A1","A2","A3","A4","B1","B2","B3","C1CMPR","C2CMPR","C3CMPR","C4CMPR","C5CMPR","D1CMPR","D2CMPR","D3CMPR","D4CMPR","E1CMPR","E2CMPR","E3CMPR","F1CMPR","F2CMPR","G1CMPR","G2CMPR","G3CMPR","G4CMPR","H1CMPR","H2CMPR"};
  
   */																																							
    /*  CombineCardMaker card(this, totYield , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
   card.Print("card_new_incertenties_R1_2Chargino_multiPoint.tab", 4);
   card.ProduceCard("card_new_incertenties_R1_2Chargino_multiPoint.tab" ,"./cards1");
   CombineCardMaker card1(this, totYield1 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
   card1.Print("card_new_incertenties_R0_2Chargino_multiPoint.tab", 4);
   card1.ProduceCard("card_new_incertenties_R0_2Chargino_multiPoint.tab" ,"./cards2");   
   CombineCardMaker card2(this, totYield2 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
   card2.Print("card_new_incertenties_I_2Chargino_multiPoint.tab", 4);
   card2.ProduceCard("card_new_incertenties_I_2Chargino_multiPoint.tab" ,"./cards7");    
   CombineCardMaker card3(this, totYield3 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
   card3.Print("card_new_incertenties_I_NocutPtLep_2Chargino_multiPoint.tab", 4);
   card3.ProduceCard("card_new_incertenties_I_NocutPtLep_2Chargino_multiPoint.tab" ,"./cards8");*/
   /*	CombineCardMaker card3(this, totYield3 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
   //	card3.Print("card_baseline_R0_2Chargino.tab", 4);
	CombineCardMaker card4(this, totYield3 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
	card4.Print("card_baseline_R0_mixt.tab", 4);
	//	CombineCardMaker card5(this, totYield3 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
	//	card5.Print("card_baseline_R0_NoChargino.tab", 4);
	card3.ProduceCard("card_I_2Chargino_150.tab" ,"./cards3");
		CombineCardMaker card4(this, totYield4 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
	card4.Print("card_A_to_H_Mod_2Chargino_120.tab", 4);
		card4.ProduceCard("card_A_to_H_Mod_2Chargino_150.tab" ,"./cards4");
	
										*/
  /*   CombineCardMaker card5(this, totYield , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
   card5.Print("card_new_incertenties_R1_mixt_multiPoint.tab", 4);
   card5.ProduceCard("card_new_incertenties_R1_mixt_multiPoint.tab" ,"./cards3");
   CombineCardMaker card6(this, totYield1 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
   card6.Print("card_new_incertenties_R0_mixt_multiPoint.tab", 4);
   card6.ProduceCard("card_new_incertenties_R0_mixt_multiPoint.tab" ,"./cards4");
   CombineCardMaker card7(this, totYield2 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
   card7.Print("card_new_incertenties_I_mixt_multiPoint.tab", 4);
   card7.ProduceCard("card_new_incertenties_I_mixt_multiPoint.tab" ,"./cards9");
   CombineCardMaker card8(this, totYield3 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
   card8.Print("card_new_incertenties_I_NocutPtLep_mixt_multiPoint.tab", 4);
   card8.ProduceCard("card_new_incertenties_I_NocutPtLep_mixt_multiPoint.tab" ,"./cards10");*/
  /*	CombineCardMaker card9(this, totYield4 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
	card9.Print("card_A_to_H_Mod_mixt_150.tab", 4);
		card9.ProduceCard("card_A_to_H_Mod_mixt_150.tab" ,"./cards9");
	
	*/
  /*	CombineCardMaker card10(this, totYield , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
	card10.Print("card_new_incertenties_R1_NoChargino_multiPoint.tab", 4);
	card10.ProduceCard("card_new_incertenties_R1_NoChargino_multiPoint.tab" ,"./cards5");
		CombineCardMaker card11(this, totYield1 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
	card11.Print("card_new_incertenties_R0_NoChargino_multiPoint.tab", 4);
	card11.ProduceCard("card_new_incertenties_R0_NoChargino_multiPoint.tab" ,"./cards6");
	CombineCardMaker card12(this, totYield2 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
	card12.Print("card_new_incertenties_I_NoChargino_multiPoint.tab", 4);
	card12.ProduceCard("card_new_incertenties_I_NoChargino_multiPoint.tab" ,"./cards11");
  	CombineCardMaker card13(this, totYield3 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
	card13.Print("card_new_incertenties_I_NocutPtLep_NoChargino_multiPoint.tab", 4);
	card13.ProduceCard("card_new_incertenties_I_NocutPtLep_NoChargino_multiPoint.tab" ,"./cards12");*/
		/*	CombineCardMaker card14(this, totYield4 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
	card14.Print("card_A_to_H_Mod_NoChargino_150.tab", 4);
		card14.ProduceCard("card_A_to_H_Mod_NoChargino_150.tab" ,"./cards14");
	
  */ 

  /* vector<string> totYield1 = {"I_equal_4_cut_Pt_300_region1","I_equal_4_cut_Pt_300_region2","I_equal_4_cut_Pt_300_region3","I_equal_4_cut_Pt_300_region4"};
     vector<string> totYield2 = {"I_equal_34_cut_Pt_300_region1","I_equal_34_cut_Pt_300_region2","I_equal_34_cut_Pt_300_region3","I_equal_34_cut_Pt_300_region4"};
     vector<string> totYield3 = {"I_equal_3_cut_Pt_300_region1","I_equal_3_cut_Pt_300_region2","I_equal_3_cut_Pt_300_region3","I_equal_3_cut_Pt_300_region4"};*/
  /*  vector<string> totYield1= {"C1","C2","C3","C4","C5"};
  vector<string> totYield2= {"C1CMPR","C2CMPR","C3CMPR","C4CMPR","C5CMPR"};
  vector<string> totYield3 = {"I1","I2","I3","I4"};
  vector<string> totYield4 = {"I108","I208","I308","I408"};
  vector<string> totYield5= {"C1","C2","C3","C4","C5","I108","I208","I308","I408"};
  vector<string> totYield6= {"C1CMPR","C2CMPR","C3CMPR","C4CMPR","C5CMPR","I1","I2","I3","I4"};

  CombineCardMaker card1(this, totYield1 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
  card1.Print("card_C_500_300_2Chargino.tab", 4);
  card1.ProduceCard("card_C_500_300_2Chargino.tab" ,"./cards1");
  CombineCardMaker card2(this, totYield2 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
  card2.Print("card_C_CMPR_500_300_2Chargino.tab", 4);
  card2.ProduceCard("card_C_CMPR_500_300_2Chargino.tab" ,"./cards2");
  CombineCardMaker card3(this, totYield3 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
  card3.Print("card_I_500_300_2Chargino.tab", 4);
  card3.ProduceCard("card_I_500_300_2Chargino.tab" ,"./cards3");
  CombineCardMaker card4(this, totYield4 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
  card4.Print("card_I_CMPR_500_300_2Chargino.tab", 4);
  card4.ProduceCard("card_I_CMPR_500_300_2Chargino.tab" ,"./cards4");
  CombineCardMaker card5(this, totYield5 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
  card5.Print("card_R1_500_300_2Chargino.tab", 4);
  card5.ProduceCard("card_R1_500_300_2Chargino.tab" ,"./cards5");
  CombineCardMaker card6(this, totYield6 , "lepChannel", "2Chargino", "StopMass", "NeutralinoMass");
  card6.Print("card_R0_500_300_2Chargino.tab", 4);
  card6.ProduceCard("card_R0_500_300_2Chargino.tab" ,"./cards6");
	
  CombineCardMaker card7(this, totYield1 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
  card7.Print("card_C_500_300_NoChargino.tab", 4);
  card7.ProduceCard("card_C_500_300_NoChargino.tab" ,"./cards7");
  CombineCardMaker card8(this, totYield2 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
  card8.Print("card_C_CMPR_500_300_NoChargino.tab", 4);
  card8.ProduceCard("card_C_CMPR_500_300_NoChargino.tab" ,"./cards8");
  CombineCardMaker card9(this, totYield3 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
  card9.Print("card_I_500_300_NoChargino.tab", 4);
  card9.ProduceCard("card_I_500_300_NoChargino.tab" ,"./cards9");
  CombineCardMaker card10(this, totYield4 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
  card10.Print("card_I_CMPR_500_300_NoChargino.tab", 4);
  card10.ProduceCard("card_I_CMPR_500_300_NoChargino.tab" ,"./cards10");
  CombineCardMaker card11(this, totYield5 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
  card11.Print("card_R1_500_300_NoChargino.tab", 4);
  card11.ProduceCard("card_R1_500_300_NoChargino.tab" ,"./cards11");
  CombineCardMaker card12(this, totYield6 , "lepChannel", "NoChargino", "StopMass", "NeutralinoMass");
  card12.Print("card_R0_500_300_NoChargino.tab", 4);
  card12.ProduceCard("card_R0_500_300_NoChargino.tab" ,"./cards12");
	
  CombineCardMaker card13(this, totYield1 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card13.Print("card_C_500_300_mixt.tab", 4);
  card13.ProduceCard("card_C_500_300_mixt.tab" ,"./cards13");
  CombineCardMaker card14(this, totYield2 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card14.Print("card_C_CMPR_500_300_mixt.tab", 4);
  card14.ProduceCard("card_C_CMPR_500_300_mixt.tab" ,"./cards14");
  CombineCardMaker card15(this, totYield3 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card15.Print("card_I_500_300_mixt.tab", 4);
  card15.ProduceCard("card_I_500_300_mixt.tab" ,"./cards15");
  CombineCardMaker card16(this, totYield4 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card16.Print("card_I_CMPR_500_300_mixt.tab", 4);
  card16.ProduceCard("card_I_CMPR_500_300_mixt.tab" ,"./cards16");
  CombineCardMaker card17(this, totYield5 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card17.Print("card_R1_500_300_mixt.tab", 4);
  card17.ProduceCard("card_R1_500_300_mixt.tab" ,"./cards17");
  CombineCardMaker card18(this, totYield6 , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
  card18.Print("card_R0_500_300_mixt.tab", 4);
  card18.ProduceCard("card_R0_500_300_mixt.tab" ,"./cards18");
  */

  /*    CombineCardMaker card3(this, totYield3 , "lepChannel", "2chargino", "StopMass", "NeutralinoMass");
	 card3.Print("card3.tab", 4);
	 
	 Same but for mixed signal
	 CombineCardMaker card2(this, totYield , "lepChannel", "mixt", "StopMass", "NeutralinoMass");
	 card2.Print("card_mixt.tab", 4);
	 card2.ProduceCard("card_mixt.tab" ,"cards");
  */ 
  






  //if you want to use the resulplot program you need to write the table with the name of the region by uncommenting the lines bellow and by putting the right totyield
    /*	   ofstream sigfile("test_I.txt");
	if (sigfile.is_open())
	{
	for(uint32_t r=0; r<totYield.size(); r++)
	{
	sigfile << totYield.at(r) << endl;
        }
	sigfile.close();
	} 
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   1"<<endl;    
   

	   ofstream sigfile1("test_I34.txt");
	if (sigfile1.is_open())
	{
	for(uint32_t r=0; r<totYield1.size(); r++)
	{
	sigfile1 << totYield1.at(r) << endl;
        }
	sigfile1.close();
	} 
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   1"<<endl;    
    */
	/*

	   ofstream sigfile2("test_I_34_cutPt_300.txt");
	if (sigfile2.is_open())
	{
	for(uint32_t r=0; r<totYield2.size(); r++)
	{
	sigfile2 << totYield2.at(r) << endl;
        }
	sigfile2.close();
	} 
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   1"<<endl;    
   


	   ofstream sigfile3("test_I_3_cutPt_300.txt");
	if (sigfile3.is_open())
	{
	for(uint32_t r=0; r<totYield3.size(); r++)
	{
	sigfile3 << totYield3.at(r) << endl;
        }
	sigfile3.close();
	} 
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   1"<<endl;    
  */ 
	/*  ofstream sigfile1("signalRegMor_pour_oral2.txt");
	if (sigfile1.is_open())
	{
	for(uint32_t r=0; r<totYield1.size(); r++)
	{
	sigfile1 << totYield1.at(r) << endl;
        }
	sigfile1.close();
	} 
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   1"<<endl;    
  
    ofstream sigfile2("signalRegMor_pour_oral3.txt");
	if (sigfile2.is_open())
	{
	for(uint32_t r=0; r<totYield3.size(); r++)
	{
	sigfile2 << totYield3.at(r) << endl;
        }
	sigfile2.close();
	} 
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   1"<<endl;   */ 
  
	/*	
	ofstream sigfile1("signalRegMor_R0.txt");
	if (sigfile1.is_open())
	{
	for(uint32_t r=0; r<totYield1.size(); r++)
	{
	sigfile1 << totYield1.at(r) << endl;
        }
	sigfile1.close();
	}
	
	
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   2"<<endl;    
	
	ofstream sigfile2("signalRegMor_A_to_H.txt");
	if (sigfile2.is_open())
	{
	for(uint32_t r=0; r<totYield2.size(); r++)
	{
	sigfile2 << totYield2.at(r) << endl;
        }
	sigfile2.close();
	} 
     
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   3"<<endl;    
	
	ofstream sigfile3("signalRegMor_I.txt");
	if (sigfile3.is_open())
	{
	for(uint32_t r=0; r<totYield3.size(); r++)
	{
	sigfile3 << totYield3.at(r) << endl;
        }
	sigfile3.close();
	}
	
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   4"<<endl;    
	
	ofstream sigfile4("signalRegMor__A_to_H_Mod.txt");
	if (sigfile4.is_open())
	{
	for(uint32_t r=0; r<totYield4.size(); r++)
	{
	sigfile4 << totYield4.at(r) << endl;
        }
	sigfile4.close();
	}
	
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   5"<<endl;    
	
	ofstream sigfile5("signalRegMor5.txt");
	if (sigfile5.is_open())
	{
	for(uint32_t r=0; r<totYield5.size(); r++)
	{
	sigfile5 << totYield5.at(r) << endl;
        }
	sigfile5.close();
	}
	
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   6"<<endl;    
	ofstream sigfile6("signalRegMor6.txt");
	if (sigfile6.is_open())
	{
	for(uint32_t r=0; r<totYield6.size(); r++)
	{
	sigfile6 << totYield6.at(r) << endl;
        }
	sigfile6.close();
	}
	
	
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   7"<<endl;    
	
	ofstream sigfile7("signalRegMor7.txt");
	if (sigfile7.is_open())
	{
	for(uint32_t r=0; r<totYield7.size(); r++)
	{
	sigfile7 << totYield7.at(r) << endl;
        }
	sigfile7.close();
	}
	
	cout<<"EEEEEEEEEEEEEEEEEEEEENNNNNNNNNNNNNNNNNNDDDDDDDDDDDDDDDDDDDDDDDDDDD  SSSSSIIIIIIGGGGGGGG   8"<<endl;    
	
	ofstream sigfile8("signalRegMor8.txt");
	if (sigfile8.is_open())
	{
	for(uint32_t r=0; r<totYield8.size(); r++)
	{
 	sigfile8 << totYield8.at(r) << endl;
        }
	sigfile8.close();
	}*/
  // moyenne->Divide(moyenne2);
  /*  TCanvas can("limits","limits");
      moyenne->Draw("colz");
      can.SaveAs("test.root");*/
  // cout<<x<<endl;
    
   cout << "end of processing" << endl;

  
}


double getWeight(string currentProcessType, float lumi, double s1fb2, bool mediumbtag)
{
  double nEvents =  myEvent.wNormalization.at(22);
  double all_weights = lumi*  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) );
  //    all_weights*= myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31));
  if (mediumbtag) all_weights *= myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ;
  else all_weights *= myEvent.weight_tightbtagsf*( nEvents / myEvent.wNormalization.at(37) ) ;
  // --> with the sytematics, the weights change but also the Normalization bins !!!! To complexify in future...
  
  if(s1fb2 != 1) all_weights *= s1fb2;
  
  if(currentProcessType == "signal")
    {
      all_weights = lumi; 
    }
  
  return all_weights;
}
