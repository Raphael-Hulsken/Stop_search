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


int main(int argc, char *argv[]){


 
  THStack *hs =new THStack("hs","Stacked 1D");
        gStyle->SetOptStat(0);
        gStyle->SetLegendBorderSize(0);
        gROOT->ForceStyle();


        string inputTab = "card_bkg.tab";
	/*	string inputTab2 = argv[2];
		string inputTab3 = argv[3];*/
        TString inputFile = "signalRegMor_bkg.txt";
	//	TString inputFile2 = argv[5];




        vector<string> regions;
        vector<string> datasets = {"bkgOneLepFromTop","bkgLostLepton","bkgOneLepFromW","bkgZnunu","totalSM"};

//name for the root plot
        vector<string> datasetsLeg = {"bkgOneLepFromTop","bkgLostLepton","bkgOneLepFromW","bkgZnunu","totalSM"};
      	

string line;
        ifstream regfile(inputFile);
        if (regfile.is_open())
        {
            while ( getline (regfile,line) )
            {
                regions.push_back(line);
            }
            regfile.close();
        }

   vector<string> regLabels = regions;
        const std::string ss = "lessMETless";
        const std::string t = "-";

        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( ss, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, ss.size(), t );
                n += t.size();
            }
        }

       TLegend *leg = new TLegend(0.6,0.65,0.85,0.85);

       //theDoctor::SonicScrewdriver sonic;
          Table tab(inputTab);
	  /*  Table tab2(inputTab2);
	      Table tab3(inputTab3);*/

       vector<TH1F*> histo;
       for(uint32_t h =0; h<datasets.size(); h++)
       {
           histo.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), regions.size(), 0, regions.size()));
	   }
       for(uint32_t r=0; r<regions.size();r++)
       {
           for(uint32_t d=0; d<datasets.size();d++)
           {
               theDoctor::Figure res = tab.Get(regions.at(r),datasets.at(d));
               histo.at(d)->SetBinContent(r+1,res.value());
               histo.at(d)->SetBinError(r+1,res.error());
               histo.at(d)->GetXaxis()->SetBinLabel(r+1,regLabels.at(r).c_str());
 
           }
       }

         
       vector<double> table_bkg_lep_top;
       string ligne2;
       ifstream fichier2("./bkg_lep_from_top.tab");			
       while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
	 
	 {
	   
	   table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
	   
        }
       fichier2.close();

	vector<double> table_bkg_lep_W;
	string ligne3;
	ifstream fichier3("./bkg_lep_from_W.tab");			
	while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
	  
	  {
	    
	    table_bkg_lep_W.push_back(stod(ligne3));  // on l'affiche
	    
	  }
       fichier3.close();
	
	vector<double> table_bkg_lost_lep;
	string ligne4;
	ifstream fichier4("./bkg_Lost_lepton.tab");			
	while(getline(fichier4, ligne4))  // tant que l'on peut mettre la ligne dans "contenu"
	  
	  {
	    
	    table_bkg_lost_lep.push_back(stod(ligne4));  // on l'affiche
	    
	  }
       fichier4.close();

	vector<double> table_bkg_Z_nunu;
	string ligne5;
	ifstream fichier5("./bkg_Z_nunu.tab");			
	while(getline(fichier5, ligne5))  // tant que l'on peut mettre la ligne dans "contenu"
	  
	  {
	    
	    table_bkg_Z_nunu.push_back(stod(ligne5));  // on l'affiche
	    
	  }
       fichier5.close();





       vector<double> table_incertitude_lep_top;
       string ligne6;
       ifstream fichier6("./incertitude_lep_from_top.tab");			
       while(getline(fichier6, ligne6))  // tant que l'on peut mettre la ligne dans "contenu"
	 
	 {
	   
	   table_incertitude_lep_top.push_back(stod(ligne6));  // on l'affiche
	   
        }
       fichier6.close();

	vector<double> table_incertitude_lep_W;
	string ligne7;
	ifstream fichier7("./incertitude_lep_from_W.tab");			
	while(getline(fichier7, ligne7))  // tant que l'on peut mettre la ligne dans "contenu"
	  
	  {
	    
	    table_incertitude_lep_W.push_back(stod(ligne7));  // on l'affiche
	    
	  }
       fichier7.close();
	
	vector<double> table_incertitude_lost_lep;
	string ligne8;
	ifstream fichier8("./incertitude_lost_lepton.tab");			
	while(getline(fichier8, ligne8))  // tant que l'on peut mettre la ligne dans "contenu"
	  
	  {
	    
	    table_incertitude_lost_lep.push_back(stod(ligne8));  // on l'affiche
	    
	  }
       fichier8.close();

	vector<double> table_incertitude_Z_nunu;
	string ligne9;
	ifstream fichier9("./incertitude_Z_nunu.tab");			
	while(getline(fichier9, ligne9))  // tant que l'on peut mettre la ligne dans "contenu"
	  
	  {
	    
	    table_incertitude_Z_nunu.push_back(stod(ligne9));  // on l'affiche
	    
	  }
       fichier9.close();



       TH1D* histo_lost_lepton = new TH1D ("lost_lepton","lost_lepton",31,0,31);
       TH1D* histo_lep_top = new TH1D ("lepton_from_top","lepton_from_top",31,0,31);
       TH1D* histo_Z_nunu = new TH1D ("Z_nunu","Z_nunu",31,0,31);
       TH1D* histo_lep_W = new TH1D ("lep_from_W","lep_from_W",31,0,31);
       

 for(uint32_t r=0; r<31;r++)
       {
         
               histo_lost_lepton->SetBinContent(r+1, table_bkg_lost_lep[r]);
               histo_lost_lepton->SetBinError(r+1,table_incertitude_lost_lep[r]*table_bkg_lost_lep[r]);
              
	       histo_lep_top->SetBinContent(r+1, table_bkg_lep_top[r]);
               histo_lep_top->SetBinError(r+1,table_incertitude_lep_top[r]*table_bkg_lep_top[r]);
	       
	       histo_lep_W->SetBinContent(r+1, table_bkg_lep_W[r]);
               histo_lep_W->SetBinError(r+1,table_incertitude_lep_W[r]*table_bkg_lep_W[r]);
	       
	       histo_Z_nunu->SetBinContent(r+1, table_bkg_Z_nunu[r]);
               histo_Z_nunu->SetBinError(r+1,table_incertitude_Z_nunu[r]*table_bkg_Z_nunu[r]);
              
       }       


       TH1D* histo_total_bkg = new TH1D ("total_bkg_paper","total_bkg_paper",31,0,31);

       histo_total_bkg->Add(histo_lost_lepton);
       histo_total_bkg->Add(histo_lep_top);
       histo_total_bkg->Add(histo_lep_W);
       histo_total_bkg->Add(histo_Z_nunu);

       TCanvas *can = new TCanvas("can","can");
       TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0);
       TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,1.0);
       pad1->Draw();        
       pad2->Draw();        


       pad1->cd();


	 pad1->SetLogy();
	 histo.at(3)->SetFillColor(0);
       histo.at(3)->SetLineColor(kGreen-3);
       histo.at(3)->Draw("e");

	 histo_Z_nunu->SetFillColor(0);
       histo_Z_nunu->SetLineColor(kBlue-3);
       histo_Z_nunu->Draw("same e");
	

 TH1D* ratio = (TH1D*) histo.at(4)->Clone() ;
	 ratio->Divide(histo_total_bkg);	
	
	 
	 auto legend = new TLegend(0.1,0.7,0.2,0.9);
	 //legend->SetHeader("Masse 500/300","C"); // option "C" allows to center the header
	 legend->AddEntry(histo.at(0), "bkg simu");
	 legend->AddEntry(histo_Z_nunu, "bkg paper");
	 //legend->AddEntry(histo.at(0), "ratio");
	 // legend->AddEntry(ratio9, "R0 T2tb");
	 //   legend->AddEntry(ratio10, "R0 T2tt");
	 
	 legend->Draw();
	 pad2->cd();
	 ratio->SetTitle("ratio bkg Total simu/paper");
	 ratio->SetLineColor(kRed-3);

	 TH1D* ratio2 =  new TH1D ("ratio2","ratio2",31,0,31);
              for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
              {
                  ratio2->SetBinContent(b+1, 1);
}
 ratio->Draw("same")  ;
	 

	      ratio2->Draw("same");
	

	 
	 
	 can->Update();
	 can->SaveAs("bkg.root"); //@MJ@ TODO svae to root file

	 double value=0, value1=0;

	 for (int i=0 ; i< 28 ; i++)
	   {
	     value += histo_total_bkg->GetBinContent(i);
	     value1 += histo.at(4)->GetBinContent(i);
	   }
	 
	 cout<<value1/value<<endl;

 
value=0;
 value1=0;
 
for (int i=28 ; i< 32 ; i++)
	   {
	     value += histo_total_bkg->GetBinContent(i);
	     value1 += histo.at(4)->GetBinContent(i);
	   }
 cout<<value1/value<<endl;

}
