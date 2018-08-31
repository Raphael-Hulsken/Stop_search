#include <exception>
#include <iostream>

#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TStyle.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"
#include "interface/Figure.h"


//./resultPlot card.tab signalRegMor.txt 
using namespace std;
using namespace theDoctor;

double max(double a ,double b)
{
  if (a>=b)
    {
      return a;
    }
  else 
    {
      return b;
    }
}


int main(int argc, char *argv[])
{

  //first I take the Stack and the superimposed plots

  TFile* f = new TFile("./plotsTest/1DSuperimposed_1point.root");//,"READ");
  TFile* f1 = new TFile("./plotsTest/1DStack_1point.root","READ");

  string a = "I_equal_34";

  string direc = "lepChannel/"+a;
  TDirectory* dir = f->GetDirectory(direc.c_str());
  TDirectory* dir1 = f1->GetDirectory(direc.c_str());
  
  string x= "pt_lead_jet";

  TCanvas* can = (TCanvas*)dir ->Get(x.c_str());
  TCanvas* can1 = (TCanvas*)dir1 ->Get(x.c_str());
  TCanvas* can2 = (TCanvas*)dir1 ->Get("yield");
  /*TH1D* hmixt = (TH1D*)can -> GetPrimitive("v:nJets|p:mixt|r:SR1l|c:lepChannel|t:1DEntries");
  TH1D* h2Chargino = (TH1D*)can -> GetPrimitive("v:nJets|p:2Chargino|r:SR1l|c:lepChannel|t:1DEntries");
  TH1D* hNoChargino = (TH1D*)can -> GetPrimitive("v:nJets|p:NoChargino|r:SR1l|c:lepChannel|t:1DEntries");
  */

  //I take it for the signal and the background
  string y="v:"+x+"|p:mixt|r:"+a+"|c:lepChannel|t:1DEntries";
  string z="v:"+x+"|p:2Chargino|r:"+a+"|c:lepChannel|t:1DEntries";
  string k="v:"+x+"|p:NoChargino|r:"+a+"|c:lepChannel|t:1DEntries";
  string bkgtot="v:"+x+"|p:bkgTotal|r:"+a+"|c:lepChannel|t:1DEntries";

   TH1D* hnbmixt = (TH1D*)can1 -> GetPrimitive(y.c_str());
  TH1D* hnb2Chargino = (TH1D*)can1 -> GetPrimitive(z.c_str());
  TH1D* hnbNoChargino = (TH1D*)can1 -> GetPrimitive(k.c_str());
  // TH1D* hnbBkgTot = (TH1D*)can1 -> GetPrimitive("v:nJets|p:bkgTotal|r:SR1l|c:lepChannel|t:1DEntries");
    THStack* hnbBkg= (THStack*) can1-> GetPrimitive("");
    /* can1->GetListOfPrimitives()->Print();
       can->GetListOfPrimitives()->Print();*/

    TH1D* hnbmixt0 = (TH1D*)can -> GetPrimitive(y.c_str());
    TH1D* hnb2Chargino0 = (TH1D*)can -> GetPrimitive(z.c_str());
    TH1D* hnbNoChargino0 = (TH1D*)can -> GetPrimitive(k.c_str());
    TH1D* hnbBkgTotal0 = (TH1D*)can -> GetPrimitive(bkgtot.c_str());
    

    string e ="v:yield|p:mixt|r:"+a+"|c:lepChannel|t:1DEntries";
    string L ="v:yield|p:2Chargino|r:"+a+"|c:lepChannel|t:1DEntries";
    string g ="v:yield|p:NoChargino|r:"+a+"|c:lepChannel|t:1DEntries";
    TH1D* hnbmixt2 = (TH1D*)can2 -> GetPrimitive(e.c_str());
    TH1D* hnb2Chargino2 = (TH1D*)can2 -> GetPrimitive(L.c_str());
    TH1D* hnbNoChargino2 = (TH1D*)can2 -> GetPrimitive(g.c_str());
    THStack* hnbBkg2 = (THStack*) can2-> GetPrimitive("");
   




    /*
  TH1D* mixt0 = hnbmixt0;
  TH1D* TChargino0 = hnb2Chargino0;
  TH1D* NoChargino0 = hnbNoChargino0;
  

   TH1D* mixt = hnbmixt;
  TH1D* TChargino = hnb2Chargino;
  TH1D* NoChargino = hnbNoChargino;
  // TH1D* bkg = new TH1D("bkg","bkg",hnbBkg->GetNbinsX(),1,10);
 
  TH1D* mixt2 = hnbmixt2;
  TH1D* TChargino2 = hnb2Chargino2;
  TH1D* NoChargino2 = hnbNoChargino2;
  



  //  TH1D* bkg = (TH1D*) hnbBkg->GetHistogram();
  /*hnbBkg->SetHistogram(bkg);*/
    /*  TH1D* bkg = (TH1D*)  hnbBkg->GetStack()->Last();;
  TH1D* bkg2 = (TH1D*)  hnbBkg2->GetStack()->Last();;
  
 for(int b=0; b< bkg->GetNbinsX(); b++)
   {
     bkg->SetBinError(b,0);
     cout<<bkg->GetBinError(b)<<endl;;
     }




    mixt->Add(bkg,-1); 
  TChargino->Add(bkg,-1);
  NoChargino->Add(bkg,-1);
 


  mixt2->Add(bkg2,-1); 
  TChargino2->Add(bkg2,-1);
  NoChargino2->Add(bkg2,-1);
 
  cout<<mixt0->GetNbinsX()<<"    "<<mixt2->GetNbinsX()<<endl;
  mixt0->Scale(mixt2->GetBinContent(1));
  TChargino0->Scale(TChargino2->GetBinContent(1));
  NoChargino0->Scale(NoChargino2->GetBinContent(1));
 
  

  TCanvas* c1=new TCanvas("c1","transparent pad",200,10,700,500);
  TPad* pad1 = new TPad("pad 1","",0,0,1,1);
 

  
  TChargino->GetXaxis()->SetTitle(x.c_str());
  TChargino->GetYaxis()->SetTitle("Entries");


  TChargino->SetLineColor(kOrange-3);
  NoChargino->SetLineColor(kGreen-3);
  mixt->SetLineColor(kBlue-3);



  TChargino2->SetLineColor(kOrange-3);
  NoChargino2->SetLineColor(kGreen-3);
  mixt2->SetLineColor(kBlue-3);

  TChargino0->GetXaxis()->SetTitle(x.c_str());
  TChargino0->GetYaxis()->SetTitle("Entries");

			
  TChargino0->SetLineColor(kOrange-3);
  NoChargino0->SetLineColor(kGreen-3);
  mixt0->SetLineColor(kBlue-3);

  pad1->Draw();
  pad1->cd();
  
  

  TChargino->Draw("e");
   NoChargino->Draw("same e");
   mixt->Draw("same e");
   //  bkg->Draw("e");

   
    auto legend = new TLegend(0.1,0.7,0.2,0.9);
   legend->SetHeader("Masse 500/300 for Baseline","C"); // option "C" allows to center the header
   legend->AddEntry(TChargino, "T2bW");
   legend->AddEntry(NoChargino, "T2tt");
   legend->AddEntry(mixt, "T2tb");
   //   legend->AddEntry(ratio10, "R0 NoChargino");
   
   legend->Draw();
   // bkg->Draw("");

   c1->Update();
   c1->SaveAs("renormalisation_test.root"); //@MJ@ TODO svae to root file
   
   TCanvas* c2=new TCanvas("c1","transparent pad",200,10,700,500);
   TPad* pad2 = new TPad("pad 1","",0,0,1,1);
   
   
   pad2->Draw();
   pad2->cd();

   double kax,Max;
   kax = max (TChargino0->GetBinContent(TChargino0->GetMaximumBin()),NoChargino0->GetBinContent(NoChargino0->GetMaximumBin()));
   Max = max (kax,mixt0->GetBinContent(mixt0->GetMaximumBin()));
     TChargino0->SetMaximum(Max+1);   
   TChargino0->Draw("e");
   NoChargino0->Draw("same e");
   mixt0->Draw("same e");
   

   auto legend1 = new TLegend(0.1,0.7,0.2,0.9);
   legend1->SetHeader("Masse 500/300 for Baseline","C"); // option "C" allows to center the header
   legend1->AddEntry(TChargino0, "T2bW");
   legend1->AddEntry(NoChargino0, "T2tt");
   legend1->AddEntry(mixt0, "T2tb");
   //   legend->AddEntry(ratio10, "R0 NoChargino");
   
   legend1->Draw();
   // bkg->Draw("");


   c2->Update();
   c2->SaveAs("renormalisation_test2.root"); //@MJ@ TODO svae to root file
   
*/

    cout<<hnbBkgTotal0->GetNbinsX()<<" "<<hnbmixt0->GetNbinsX()<<endl;
      double re=0, ke=0;
      double erorS=0, erorB=0, erorInter=0, erorTot=0;
      TH1D* fin = new TH1D("efficiency of the cut","efficiency of the cut",hnbBkgTotal0->GetNbinsX(),0,hnbBkgTotal0->GetNbinsX()*50);
      


      //this one is if you want a cut that is > the cut 
      /*     for(int te=1 ; te<= hnbBkgTotal0->GetNbinsX(); te++)
     {
       re+=hnbmixt0->GetBinContent(te);
       ke+= hnbBkgTotal0->GetBinContent(te) ;
       erorS = sqrt(erorS*erorS + hnbmixt0->GetBinError(te)/re* hnbmixt0->GetBinError(te)/re);
       // cout<<"pour le signal l'erreur est de "<<hnbmixt0->GetBinErroro(te)<<" est la valeur de "<<re<<" ce qui donne une erreur relative de "<<erorS<<endl;
       erorB = sqrt(erorS*erorS + hnbBkgTotal0->GetBinError(te)/ke * hnbBkgTotal0->GetBinError(te)/ke) ;
       //  cout<<"pour le bruit de fond l'erreur est de "<<hnbBkgTotal0->GetBinError(te)<<" est la valeur de "<<ke<<" ce qui donne une erreur relative de "<<erorB<<endl;
       erorInter = sqrt(erorS*erorS/ke + re*re*erorB*erorB/(4*ke*ke*ke));
       // erorTot = sqrt(erorTot*erorTot + erorInter*erorInter);
       cout<<"re="<<re<<", ke="<<ke<<", erorS="<<erorS<<", erorB="<<erorB<<", erorInter="<<erorInter<<", erorTot="<<erorTot<<endl;
       //       cout<<re<<" "<<ke<<endl;
       //cout<<"pour le bin "<<te*0.1<<" on a le ratio "<<(re)/sqrt(ke)<<endl;
       //cout<<"l'erreur est "<<erorTot<<endl;
       if (ke <  1)
       {
	 //cout<<ke<<endl;
	 fin->SetBinContent(te,(re)/sqrt(ke));
	 fin->SetBinError(te,erorInter*(re)/sqrt(ke));
       }
     }
     cout<<e<<" "<<ke<<endl;
      */


      //this one is if you want a cut that is < the cut 
        for(int te= hnbBkgTotal0->GetNbinsX() ; te>0 ; te--)
     {
       re+=hnbmixt0->GetBinContent(te);
       ke+= hnbBkgTotal0->GetBinContent(te) ;
       erorS = sqrt(erorS*erorS + hnbmixt0->GetBinError(te)/re* hnbmixt0->GetBinError(te)/re);
       // cout<<"pour le signal l'erreur est de "<<hnbmixt0->GetBinErroro(te)<<" est la valeur de "<<re<<" ce qui donne une erreur relative de "<<erorS<<endl;
       erorB = sqrt(erorS*erorS + hnbBkgTotal0->GetBinError(te)/ke * hnbBkgTotal0->GetBinError(te)/ke) ;
       //  cout<<"pour le bruit de fond l'erreur est de "<<hnbBkgTotal0->GetBinError(te)<<" est la valeur de "<<ke<<" ce qui donne une erreur relative de "<<erorB<<endl;
       erorInter = sqrt(erorS*erorS/ke + re*re*erorB*erorB/(4*ke*ke*ke));
       // erorTot = sqrt(erorTot*erorTot + erorInter*erorInter);
       cout<<"re="<<re<<", ke="<<ke<<", erorS="<<erorS<<", erorB="<<erorB<<", erorInter="<<erorInter<<", erorTot="<<erorTot<<endl;
       //       cout<<re<<" "<<ke<<endl;
       //cout<<"pour le bin "<<te*0.1<<" on a le ratio "<<(re)/sqrt(ke)<<endl;
       //cout<<"l'erreur est "<<erorTot<<endl;
       if (ke <  1)
       {
	 //cout<<ke<<endl;
	 fin->SetBinContent(te,(re)/sqrt(ke));
	 fin->SetBinError(te,erorInter*(re)/sqrt(ke));
       }
       }
       cout<<e<<" "<<ke<<endl;

   TCanvas* c3=new TCanvas("c1","transparent pad",200,10,700,500);
   TPad* pad3 = new TPad("pad 1","",0,0,1,1);
   pad3->Draw();
   pad3->cd();
   
	 TH1D* ratio2 =  new TH1D ("ratio2","ratio2",hnbBkgTotal0->GetNbinsX(),0,hnbBkgTotal0->GetNbinsX()*50);
              for(uint32_t b=0; b< hnbBkgTotal0->GetNbinsX(); b++)
              {
		ratio2->SetBinContent(b+1, 1);
	      }


   fin->GetXaxis()->SetTitle(x.c_str());
   fin->GetYaxis()->SetTitle("ratio S/sqrt(B)");
   fin->SetLineColor(kRed);
   fin->Draw("TEXT45");
   ratio2->SetLineColor(kBlack);
   ratio2->Draw("same");
   TLine * loose = new TLine(0.9535,0.3,0.9535,1.2);
   loose->SetLineColor(kBlue);
   //  loose->Draw("same") ;  
   
   TLine * med = new TLine(0.8484,0.3,0.8484,1.2);
   med->SetLineColor(kGreen);
   //  med->Draw("same") ;  
   

      TLine * tight = new TLine(0.5426,0.3,0.5426,1.2);
   tight->SetLineColor(kRed);
   //  tight->Draw("same") ;  
   
   TLine * ptlep = new TLine(150,0.3,150,1.2);
   ptlep->SetLineColor(kBlue);
   //   ptlep->Draw("same") ;  
   
c3->Update();
   c3->SaveAs("renormalisation_test3.root"); //@MJ@ TODO svae to root file
     
}
