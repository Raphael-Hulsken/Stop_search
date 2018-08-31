#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
#include <cmath>
#include <exception>
#include <ctime>

#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THashList.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TLegend.h"

using namespace std;

//this function give the value for every point even the one were we don't have data. It fills the gap
TGraph2D* GetInterpolatingGraph(TH2F *hold){
  float binsize = hold->GetXaxis()->GetBinWidth(1)/2.;
  TString name = hold->GetName();
  name.ReplaceAll("Org","");
  TGraph2D *g = new TGraph2D(hold);
  g->SetNpx(int(g->GetXmax()-g->GetXmin())/binsize);
  g->SetNpy(int(g->GetYmax()-g->GetYmin())/binsize);
  return g;
}

//This one create the limit for a given value here the value is 1
TGraph* GetContour(TGraph2D *g, TString name, TGraph *gempty){
  TGraph *gnew;
  TH2D *temp = (TH2D*)g->GetHistogram();//need this for list to work?
  //g->Draw("alp");//need this for list to work?
  TList *glist = (TList*)g->GetContourList(1.0);
  if(glist == nullptr) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  if(glist->GetSize()==0) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  int max_points = -1;
  int nn = glist->GetSize();
  for(int i = 0; i<glist->GetSize(); ++i){
    TGraph *gtemp = (TGraph*)glist->At(i);
    int Npoints = gtemp->GetN();
    if(Npoints>max_points){
      gnew = (TGraph*)gtemp->Clone(name);
      max_points = Npoints;
      if(gnew->GetN()>0){
	double x,y;
	gnew->GetPoint(0,x,y);
	if(y<13.) gnew->SetPoint(0,x,0);
	gnew->GetPoint(gnew->GetN()-1,x,y);
	if(y<13.) { gnew->SetPoint(gnew->GetN(),x,0); }//they way ROOT does the interpolation the LSP = 0 strip becomes empty
      }
    }
  }
  return gnew;
}


//This one create the limit for a given value here the value is 2.07, the second function is usefull when you want two type of limit
TGraph* GetContour2(TGraph2D *g, TString name, TGraph *gempty){
  TGraph *gnew;
  TH2D *temp = (TH2D*)g->GetHistogram();//need this for list to work?
  //g->Draw("alp");//need this for list to work?
  TList *glist = (TList*)g->GetContourList(2.07);
  if(glist == nullptr) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  if(glist->GetSize()==0) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  int max_points = -1;
  int nn = glist->GetSize();
  for(int i = 0; i<glist->GetSize(); ++i){
    TGraph *gtemp = (TGraph*)glist->At(i);
    int Npoints = gtemp->GetN();
    if(Npoints>max_points){
      gnew = (TGraph*)gtemp->Clone(name);
      max_points = Npoints;
      if(gnew->GetN()>0){
	double x,y;
	gnew->GetPoint(0,x,y);
	if(y<13.) gnew->SetPoint(0,x,0);
	gnew->GetPoint(gnew->GetN()-1,x,y);
	if(y<13.) { gnew->SetPoint(gnew->GetN(),x,0); }//they way ROOT does the interpolation the LSP = 0 strip becomes empty
      }
    }
  }
  return gnew;
}



//if ever you want a third one (especially if you want to plots the significance limit)
TGraph* GetContour3(TGraph2D *g, TString name, TGraph *gempty){
  TGraph *gnew;
  TH2D *temp = (TH2D*)g->GetHistogram();//need this for list to work?
  //g->Draw("alp");//need this for list to work?
  TList *glist = (TList*)g->GetContourList(5.0);
  if(glist == nullptr) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  if(glist->GetSize()==0) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  int max_points = -1;
  int nn = glist->GetSize();
  for(int i = 0; i<glist->GetSize(); ++i){
    TGraph *gtemp = (TGraph*)glist->At(i);
    int Npoints = gtemp->GetN();
    if(Npoints>max_points){
      gnew = (TGraph*)gtemp->Clone(name);
      max_points = Npoints;
      if(gnew->GetN()>0){
	double x,y;
	gnew->GetPoint(0,x,y);
	if(y<13.) gnew->SetPoint(0,x,0);
	gnew->GetPoint(gnew->GetN()-1,x,y);
	if(y<13.) { gnew->SetPoint(gnew->GetN(),x,0); }//they way ROOT does the interpolation the LSP = 0 strip becomes empty
      }
    }
  }
  return gnew;
}



/*TGraph* GetContour2(TGraph2D *g, TString name, TGraph *gempty){
  TGraph *gnew;
  TH2D *temp = (TH2D*)g->GetHistogram();//need this for list to work?
  //g->Draw("alp");//need this for list to work?
  TList *glist = (TList*)g->GetContourList(1.0);
  if(glist == nullptr) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  if(glist->GetSize()==0) {gnew = (TGraph*)gempty->Clone(name); return gnew; }
  int max_points = -1;
  int nn = glist->GetSize();
  int noti = -1;
  //cout << "number of entries in list " << nn << endl;
  for(int i = 0; i<glist->GetSize(); ++i){
    TGraph *gtemp = (TGraph*)glist->At(i);
    int Npoints = gtemp->GetN();
    if(Npoints>max_points){
      gnew = (TGraph*)gtemp->Clone(name);
      max_points = Npoints;
      noti = i;
      if(gnew->GetN()>0){
	double x,y;
	gnew->GetPoint(0,x,y);
	if(y<13.) gnew->SetPoint(0,x,0);
	gnew->GetPoint(gnew->GetN()-1,x,y);
	if(y<13.) { gnew->SetPoint(gnew->GetN(),x,0);}//they way ROOT does the interpolation the LSP = 0 strip becomes empty
      }
    }
  }
  max_points = -1;
  for(int i = 0; i<glist->GetSize(); ++i){
    if(i==noti&&glist->GetSize()>1) continue;
    TGraph *gtemp = (TGraph*)glist->At(i);
    int Npoints = gtemp->GetN();
    if(Npoints>max_points){
      gnew = (TGraph*)gtemp->Clone(name);
      max_points = Npoints;
      if(gnew->GetN()>0){
	double x,y;
	gnew->GetPoint(0,x,y);
	if(y<13.) gnew->SetPoint(0,x,0);
	gnew->GetPoint(gnew->GetN()-1,x,y);
	if(y<13.) { gnew->SetPoint(gnew->GetN(),x,0);}//they way ROOT does the interpolation the LSP = 0 strip becomes empty
      }
    }
  }
  return gnew;
  }*/

int main(int argc, char *argv[])
{

    if(argc != 5)
        throw std::runtime_error("Bad number of arguments!");

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kMint);
    gROOT->ForceStyle();

    TString file1 = argv[1];
    TString file2 = argv[2];
      TString file3 = argv[3];
        TString file4 = argv[4];
	//TString file5 = argv[5];
    
    TGraph *gEmpty = new TGraph();

    TFile* f1 = NULL;
    TTree* t1 = NULL;
    f1 = TFile::Open(file1);
    if(f1==NULL)
        throw std::runtime_error("File 1 address not set");
    t1 = dynamic_cast<TTree*>(f1->Get("limit"));
    if(t1==NULL)
        throw std::runtime_error("Tree 1 address not set");


    TFile* f2 = NULL;
    TTree* t2 = NULL;
    f2 = TFile::Open(file2);
    if(f2==NULL)
        throw std::runtime_error("File 2 address not set");
    t2 = dynamic_cast<TTree*>(f2->Get("limit"));
    if(t2==NULL)
        throw std::runtime_error("Tree 2 address not set");

    
    
         TFile* f3 = NULL;
    TTree* t3 = NULL;
    f3 = TFile::Open(file3);
    if(f3==NULL)
        throw std::runtime_error("File 3 address not set");
    t3 = dynamic_cast<TTree*>(f3->Get("limit"));
    if(t3==NULL)
        throw std::runtime_error("Tree 3 address not set");


    
    TFile* f4 = NULL;
    TTree* t4 = NULL;
    f4 = TFile::Open(file4);
    if(f4==NULL)
        throw std::runtime_error("File 4 address not set");
    t4 = dynamic_cast<TTree*>(f4->Get("limit"));
    if(t4==NULL)
        throw std::runtime_error("Tree 4 address not set");


    /*

    TFile* f5 = NULL;
    TTree* t5 = NULL;
    f5 = TFile::Open(file5);
    if(f5==NULL)
        throw std::runtime_error("File 5 address not set");
    t5 = dynamic_cast<TTree*>(f5->Get("limit"));
    if(t5==NULL)
        throw std::runtime_error("Tree 4 address not set");

    */

    double limit = -13;
    float mSbot = -13;
    float mNeutr = -13;
    t1->SetBranchAddress("limit", &limit );
    t1->SetBranchAddress("spart_mass_point", &mSbot );
    t1->SetBranchAddress("neutr_mass_point", &mNeutr );

   

    t2->SetBranchAddress("limit", &limit );
    t2->SetBranchAddress("spart_mass_point", &mSbot );
    t2->SetBranchAddress("neutr_mass_point", &mNeutr );

    
        
    t3->SetBranchAddress("limit", &limit );
    t3->SetBranchAddress("spart_mass_point", &mSbot );
    t3->SetBranchAddress("neutr_mass_point", &mNeutr );

    
    
    t4->SetBranchAddress("limit", &limit );
    t4->SetBranchAddress("spart_mass_point", &mSbot );
    t4->SetBranchAddress("neutr_mass_point", &mNeutr );
    /*
    t5->SetBranchAddress("limit", &limit );
    t5->SetBranchAddress("spart_mass_point", &mSbot );
    t5->SetBranchAddress("neutr_mass_point", &mNeutr );
    */
    TH2F* limitPlot1 = new TH2F("ratio section efficace exclue/section efficace theorique", "ratio section efficace exclue/section efficace theorique", 150, 0 , 1500, 130, 0, 1300); //@MJ@ TODO do smething more general
    // TH2F* limitPlot1 = new TH2F("ratio section efficace d'exclusion/section efficace theorique", "ratio section efficace d'exclusion/section efficace theorique", 40, 0 , 400, 3, 20, 50); //@MJ@ TODO do smething more general
    TH2F* limitPlot2 = new TH2F("limitPlot1", "limitPlot1", 150, 0 , 1500, 130, 0, 1300); //@MJ@ TODO do smething more general
    limitPlot1->GetXaxis()->SetTitle("masse du stop (GeV)");   
         TH2F* limitPlot3 = new TH2F("limitPlot2", "limitPlot2", 150, 0 , 1500, 130, 0, 1300); //@MJ@ TODO do smething more general
        TH2F* limitPlot4 = new TH2F("limitPlot3", "limitPlot3", 150, 0 , 1500, 130, 0, 1300); //@MJ@ TODO do smething more general
	//  TH2F* limitPlot5 = new TH2F("limitPlot4", "limitPlot4", 150, 0 , 1500, 130, 0, 1300); //@MJ@ TODO do smething more general
    
    Int_t nentries1 = (Int_t)t1->GetEntries();
    Int_t nentries2 = (Int_t)t2->GetEntries();
    cout<<nentries1<<endl;
     Int_t nentries3 = (Int_t)t3->GetEntries();
         Int_t nentries4 = (Int_t)t4->GetEntries();
	 // Int_t nentries5 = (Int_t)t5->GetEntries();
    
    
 


    TCanvas can("ratio section efficace d'exclusion/section efficace theorique pour A a H","ratio section efficace d'exclusion/section efficace theorique pour A a H");
    //limitPlot->Draw("colz");
    /*  h =can.DrawFrame(0.,0.,1.,1.)
	h->SetXTitle("x")      */     

 
   
    
    //take the e==2 or (e-2)%6==0 in order to have the 95% exclusion
    for (Int_t e=0; e<nentries1; e++)
      {
	t1->GetEntry(e);
	
		if(e==2 || (e-2)%6==0)
	{
	  //   cout << "entry " << e << " limit " <<limit << " for point: sbottom " << mSbot << " neutralino " << mNeutr << endl;
	    //if(limit < 1)
	    limitPlot1->Fill(mSbot,mNeutr,limit);
	    
	     }
      }

    // cout<<file1[file1.Sizeof()-7]<<endl;

     for (Int_t e=0; e<nentries2; e++)
     {
       	t2->GetEntry(e);
	
		if(e==2 || (e-2)%6==0)
		 {
	//	cout << "entry " << e << " limit " <<limit << " for point: sbottom " << mSbot << " neutralino " << mNeutr << endl;
	//if(limit < 1)
	limitPlot2->Fill(mSbot,mNeutr,limit);
	/*	if (mSbot==500 && mNeutr ==300)
		cout<<"la signif est :"<<limit<<endl;	    */
	  }
    	  
     }

         for (Int_t e=0; e<nentries3; e++)
      {
	t3->GetEntry(e);
	
	if(e==2 || (e-2)%6==0)
	  {
	    //   cout << "entry " << e << " limit " <<limit << " for point: sbottom " << mSbot << " neutralino " << mNeutr << endl;
	    //if(limit < 1)
	    limitPlot3->Fill(mSbot,mNeutr,limit);
	    
	  }
      }
    
    
	 
    for (Int_t e=0; e<nentries4; e++)
      {
	t4->GetEntry(e);
	
	if(e==2 || (e-2)%6==0)
	  {
	    //  cout << "entry " << e << " limit " <<limit << " for point: sbottom " << mSbot << " neutralino " << mNeutr << endl;
	    //if(limit < 1)
	    limitPlot4->Fill(mSbot,mNeutr,limit);
	    
	  }
      }

    /*

    for (Int_t e=0; e<nentries5; e++)
      {
	t5->GetEntry(e);
	
	if(e==2 || (e-2)%6==0)
	  {
	    cout << "entry " << e << " limit " <<limit << " for point: sbottom " << mSbot << " neutralino " << mNeutr << endl;
	    //if(limit < 1)
	    limitPlot5->Fill(mSbot,mNeutr,limit);
	    
	  }
	  }
    */
    // limitPlot1->Draw("colz");
    
   
    
    TGraph2D *g2Exp1 = (TGraph2D*)GetInterpolatingGraph(limitPlot1);
    //    g2Exp1->SetTitle("Donnee","masse du stop (GeV)","masse du neutralino (GeV)","colz");   
    TGraph *gExp_c1 = (TGraph*)GetContour(g2Exp1, "gExp", gEmpty);
    

   
    TGraph2D *g2Exp2 = (TGraph2D*)GetInterpolatingGraph(limitPlot2);
    
    TGraph *gExp_c2 = (TGraph*)GetContour(g2Exp2, "gExp", gEmpty);

    //   TGraph *gExp_c2a = (TGraph*)GetContour3(g2Exp2, "gExp", gEmpty);
    

       TGraph2D *g2Exp3 = (TGraph2D*)GetInterpolatingGraph(limitPlot3);
    
    TGraph *gExp_c3 = (TGraph*)GetContour(g2Exp3, "gExp1", gEmpty);
    
    
    TGraph2D *g2Exp4 = (TGraph2D*)GetInterpolatingGraph(limitPlot4);
    
    TGraph *gExp_c4 = (TGraph*)GetContour(g2Exp4, "gExp3", gEmpty);
    /*
    TGraph2D *g2Exp5 = (TGraph2D*)GetInterpolatingGraph(limitPlot5);
    
    TGraph *gExp_c5 = (TGraph*)GetContour(g2Exp5, "gExp4", gEmpty);
    */
  
    gExp_c1->GetXaxis()->SetTitle("Masse du stop (GeV)");
    gExp_c1->GetYaxis()->SetTitle("Masse du neutralino (GeV)");

    g2Exp1->Draw("colz");
    gExp_c1->SetLineColor(kGreen-8);
    // gExp_c1->SetLineWidth(4);
    gExp_c1->SetLineWidth(4);
    gExp_c1->Draw("same");
     

      gExp_c3->SetLineColor(kGreen);
 gExp_c3->SetLineWidth(4);
  gExp_c3->Draw("same");
    gExp_c4->SetLineColor(kRed);
 gExp_c4->SetLineWidth(4);
  gExp_c4->Draw("same");
    gExp_c2->SetLineColor(kRed-2);
 gExp_c2->SetLineWidth(4);
 //Exp_c2a->SetLineColor(2);
 //Exp_c2a->SetLineWidth(4);
    //gExp_c2->SetLineWidth(4);
  gExp_c2->Draw("same"); 
  //Exp_c2a->Draw("same");

    /*  gExp_c5->SetLineColor(7);
      gExp_c5->SetLineWidth(4);
   gExp_c5->Draw("same");
    */
    

   auto legend = new TLegend(0.1,0.7,0.48,0.9);
   /*  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(gExp_c1, "limit pour l'analyse CMS");
     legend->AddEntry(gExp_c2, "decouverte 3 sigma nouvelle luminosite");
        legend->AddEntry(gExp_c2a, "decouverte 5 sigma nouvelle luminosite");
   //  legend->AddEntry(gExp_c1,"Limite pour la premiere etude");
   // legend->AddEntry(gExp_c2,"Limite pour l'etude avec la region compressee");
   */

   //can.GetCSetTitle("test");
   legend->AddEntry(gExp_c1, "limit for A to H NoChargino");// L=35.9 fb^{-1}" );
   legend->AddEntry(gExp_c2, "limit for I NoChargino");//A to H pm 1 sigma");// L=120 fb^{-1}");
   legend->AddEntry(gExp_c3, "limit for I 2Chargino");// L=120 fb^{-1}");
   legend->AddEntry(gExp_c4, "limit for I mixt");
 //  legend->AddEntry(gExp_c5, "limit for I");
   legend->Draw();

   //TF1* funck = new TF1 ("fonction","x-250",200,1400); 
    //funck->Draw("same");

   
can.SaveAs("limitPlot.root"); 
    
    
    //TGraph *gExp_c2 = (TGraph*)GetContour2(g2Exp, "gExp_2_c2", gEmpty);
    
   // gExp_c2->Draw("same");
    

    
} 

