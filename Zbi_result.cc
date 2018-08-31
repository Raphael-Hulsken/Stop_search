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


//------------------
//Function to compute the significance
//------------------
double Zbi(double n_sig, double n_b, double rel_uncert_b = 0.4 )
{
  //if (n_sig==0){return 0;}
  // else 
  //{
  double n_on = n_sig+n_b;
  //cout<<"n_on ="<<n_on<<endl;
  double mu_b_hat = n_b;
  //cout<<"mu_b_hat ="<<mu_b_hat<<endl;
  double sigma_b = rel_uncert_b*n_b;
  //cout<<"sigma_b ="<<sigma_b<<endl;
  double tau = mu_b_hat/(sigma_b*sigma_b);
  // cout<<"tau ="<<tau<<endl;
  double n_off = tau*mu_b_hat;
  //cout<<"n_off ="<<n_off<<endl;
  double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);//TMath::Beta(n_on,1+n_off);
  //  cout<<"P_Bi ="<<P_Bi<<endl;
  double Z_Bi = sqrt(2)*TMath::ErfInverse(1 - 2*P_Bi);
  //  cout<<"Z_Bi ="<<Z_Bi<<endl;
  if (Z_Bi<=0)
    {
      return 0;
    }
  else 
    { 
      return Z_Bi;
    }
  //}
}



int main(int argc, char *argv[]){


  //small code that is written the same way as the resultPlots but here we study only one set of region and it compute some specific value for Zbi
   
  if(argc != 3)
    throw std::runtime_error("Bad number of arguments!");
  
  string inputTab = argv[1];
   TString inputFile = argv[2];
  
  vector<string> regions;
  vector<string> datasets = {"bkgOneLepFromTop","bkgLostLepton","bkgOneLepFromW","bkgZnunu", "(600,400)","(800,600)","(200,0)","(300,100)","data1","totalSM"};
  
  
  
  //input for the error per region
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

  Table tab(inputTab);
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
	  // histo.at(d)->GetXaxis()->SetBinLabel(r+1,regLabels.at(r).c_str());
	  
	}
    }


  cout<<"for 600,400 : in the regions I1, Zbi="<<Zbi(histo.at(4)->GetBinContent(1),histo.at(9)->GetBinContent(1))<<endl;
  cout<<"for 600,400 : in the regions I14, Zbi="<<Zbi(histo.at(4)->GetBinContent(5),histo.at(9)->GetBinContent(5))<<endl;
  cout<<"for 600,400 : in the regions I2, Zbi="<<Zbi(histo.at(4)->GetBinContent(2),histo.at(9)->GetBinContent(2))<<endl;
  cout<<"for 600,400 : in the regions I24, Zbi="<<Zbi(histo.at(4)->GetBinContent(6),histo.at(9)->GetBinContent(6))<<endl;
  cout<<"for 600,400 : in the regions I3, Zbi="<<Zbi(histo.at(4)->GetBinContent(3),histo.at(9)->GetBinContent(3))<<endl;
  cout<<"for 600,400 : in the regions I34, Zbi="<<Zbi(histo.at(4)->GetBinContent(7),histo.at(9)->GetBinContent(7))<<endl;
  cout<<"for 600,400 : in the regions I4, Zbi="<<Zbi(histo.at(4)->GetBinContent(4),histo.at(9)->GetBinContent(4))<<endl;
  cout<<"for 600,400 : in the regions I44, Zbi="<<Zbi(histo.at(4)->GetBinContent(8),histo.at(9)->GetBinContent(8))<<endl;


  cout<<"for 800,650 : in the regions I1, Zbi="<<Zbi(histo.at(5)->GetBinContent(1),histo.at(9)->GetBinContent(1))<<endl;
  cout<<"for 800,650 : in the regions I14, Zbi="<<Zbi(histo.at(5)->GetBinContent(5),histo.at(9)->GetBinContent(5))<<endl;
  cout<<"for 800,650 : in the regions I2, Zbi="<<Zbi(histo.at(5)->GetBinContent(2),histo.at(9)->GetBinContent(2))<<endl;
  cout<<"for 800,650 : in the regions I24, Zbi="<<Zbi(histo.at(5)->GetBinContent(6),histo.at(9)->GetBinContent(6))<<endl;
  cout<<"for 800,650 : in the regions I3, Zbi="<<Zbi(histo.at(5)->GetBinContent(3),histo.at(9)->GetBinContent(3))<<endl;
  cout<<"for 800,650 : in the regions I34, Zbi="<<Zbi(histo.at(5)->GetBinContent(7),histo.at(9)->GetBinContent(7))<<endl;
  cout<<"for 800,650 : in the regions I4, Zbi="<<Zbi(histo.at(5)->GetBinContent(4),histo.at(9)->GetBinContent(4))<<endl;
  cout<<"for 800,650 : in the regions I44, Zbi="<<Zbi(histo.at(5)->GetBinContent(8),histo.at(9)->GetBinContent(8))<<endl;

 cout<<"for 200,0 : in the regions I1, Zbi="<<Zbi(histo.at(6)->GetBinContent(1),histo.at(9)->GetBinContent(1))<<endl;
  cout<<"for 200,0 : in the regions I14, Zbi="<<Zbi(histo.at(6)->GetBinContent(5),histo.at(9)->GetBinContent(5))<<endl;
  cout<<"for 200,0 : in the regions I2, Zbi="<<Zbi(histo.at(6)->GetBinContent(2),histo.at(9)->GetBinContent(2))<<endl;
  cout<<"for 200,0 : in the regions I24, Zbi="<<Zbi(histo.at(6)->GetBinContent(6),histo.at(9)->GetBinContent(6))<<endl;
  cout<<"for 200,0 : in the regions I3, Zbi="<<Zbi(histo.at(6)->GetBinContent(3),histo.at(9)->GetBinContent(3))<<endl;
  cout<<"for 200,0 : in the regions I34, Zbi="<<Zbi(histo.at(6)->GetBinContent(7),histo.at(9)->GetBinContent(7))<<endl;
  cout<<"for 200,0 : in the regions I4, Zbi="<<Zbi(histo.at(6)->GetBinContent(4),histo.at(9)->GetBinContent(4))<<endl;
  cout<<"for 200,0 : in the regions I44, Zbi="<<Zbi(histo.at(6)->GetBinContent(8),histo.at(9)->GetBinContent(8))<<endl;

 cout<<"for 300,150 : in the regions I1, Zbi="<<Zbi(histo.at(7)->GetBinContent(1),histo.at(9)->GetBinContent(1))<<endl;
  cout<<"for 300,150 : in the regions I14, Zbi="<<Zbi(histo.at(7)->GetBinContent(5),histo.at(9)->GetBinContent(5))<<endl;
  cout<<"for 300,150 : in the regions I2, Zbi="<<Zbi(histo.at(7)->GetBinContent(2),histo.at(9)->GetBinContent(2))<<endl;
  cout<<"for 300,150 : in the regions I24, Zbi="<<Zbi(histo.at(7)->GetBinContent(6),histo.at(9)->GetBinContent(6))<<endl;
  cout<<"for 300,150 : in the regions I3, Zbi="<<Zbi(histo.at(7)->GetBinContent(3),histo.at(9)->GetBinContent(3))<<endl;
  cout<<"for 300,150 : in the regions I34, Zbi="<<Zbi(histo.at(7)->GetBinContent(7),histo.at(9)->GetBinContent(7))<<endl;
  cout<<"for 300,150 : in the regions I4, Zbi="<<Zbi(histo.at(7)->GetBinContent(4),histo.at(9)->GetBinContent(4))<<endl;
  cout<<"for 300,150 : in the regions I44, Zbi="<<Zbi(histo.at(7)->GetBinContent(8),histo.at(9)->GetBinContent(8))<<endl;
  cout<<histo.at(9)->GetBinContent(1)<<"          "<<histo.at(7)->GetBinContent(1)<<endl;





}
