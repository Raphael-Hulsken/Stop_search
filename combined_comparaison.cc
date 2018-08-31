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


  //code that work the same as resulsPlot.cc but here the output is the number of event in the regions, I used it to compute the number of event in R1 and R0 to be sure that there are the same and that I do not count some event twice


  // THStack *hs =new THStack("hs","Stacked 1D");
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);
  gROOT->ForceStyle();
  
  if(argc != 3)
    throw std::runtime_error("Bad number of arguments!");
  
  string inputTab = argv[1];
  // string inputTab2 = argv[2];
  TString inputFile = argv[2];
  /* TString inputFile2 = argv[4];
  string inputTab3 = argv[5];
  TString inputFile3 = argv[6];
 string inputTab4 = argv[7];
 TString inputFile4 = argv[8];

  */
  vector<string> regions; 
  /* vector<string> regions2;
 vector<string> regions3;
 vector<string> regions4;*/
  vector<string> datasets = {"bkgOneLepFromTop","bkgLostLepton","bkgOneLepFromW","bkgZnunu", "(600,325)","(300,50)","(900,50)","data1","totalSM"};
  
  
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
    


  /*
string line2;
  ifstream regfile2(inputFile2);
  if (regfile2.is_open())
    {
      while ( getline (regfile2,line2) )
	{
	  regions2.push_back(line2);
	}
      regfile2.close();
    }



  string line3;
  ifstream regfile3(inputFile3);
  if (regfile3.is_open())
    {
      while ( getline (regfile3,line3) )
	{
	  regions3.push_back(line3);
	}
      regfile3.close();
    }


  string line4;
  ifstream regfile4(inputFile4);
  if (regfile4.is_open())
    {
      while ( getline (regfile4,line4) )
	{
	  regions4.push_back(line4);
	}
      regfile4.close();
    }

  */
  
  //theDoctor::SonicScrewdriver sonic;
  Table tab(inputTab);
  /*  Table tab2(inputTab2);
  Table tab3(inputTab3);
  Table tab4(inputTab4);*/
  //  double y=0;
  for(uint32_t r=0; r<regions.size();r++)
    {
       theDoctor::Figure res = tab.Get(regions.at(r),datasets.at(8));
       cout<<"in the regions "<< regions.at(r)<<" there are "<<res.value()<<" background"<<endl;
       // theDoctor::Figure res2 = tab2.Get(regions2.at(r),datasets.at(8));
       // cout<<"between the region "<<regions.at(r)<<" and "<<regions2.at(r)<<" we have a difference of "<<res.value()-res2.value()<<endl;
       // y+=res.value()-res2.value();
    }


   //cout<<"l'ecart total est de "<<y<<endl;
   /*
   double z=0;
   double k=0;
   // cout<<regions.size()<<" "<<regions2.size()<<endl;
   // if (regions2.size()>regions.size())
   //{
   for(uint32_t r=regions.size()-4; r<regions.size();r++)
	 {
	   theDoctor::Figure res = tab.Get(regions.at(r),datasets.at(8));
	   theDoctor::Figure res2 = tab2.Get(regions2.at(r),datasets.at(8));
	   //  cout<<"in the compressed region "<<regions2.at(r)<<" there are "<<res2.value()<<" events, and in the compressed regions "<<regions.at(r)<<" there are "<<res.value()<<" events"<<endl;
	   z+=res2.value();
	   k+=res.value();
	 }
       //}


   /* cout<<z<<" is the number of events in the total compressed regions"<<endl;
   cout<<k<<" is the number of events in the total compressed regions without repetition"<<endl;
   cout<<"the difference is "<<z-(k+y)<<endl;*/


   //double l=0;
   /*   // double m=0;

  for(uint32_t r=0; r<regions.size();r++)
	 {
	   theDoctor::Figure res = tab2.Get(regions2.at(r),datasets.at(0))-tab.Get(regions.at(r),datasets.at(0));
	   theDoctor::Figure res1 = tab2.Get(regions2.at(r),datasets.at(1))-tab.Get(regions.at(r),datasets.at(1));
	   theDoctor::Figure res2 = tab2.Get(regions2.at(r),datasets.at(2))-tab.Get(regions.at(r),datasets.at(2));
	   theDoctor::Figure res3 = tab2.Get(regions2.at(r),datasets.at(3))-tab.Get(regions.at(r),datasets.at(3));
	   theDoctor::Figure res4 = tab2.Get(regions2.at(r),datasets.at(4))-tab.Get(regions.at(r),datasets.at(4));
	   theDoctor::Figure res5 = tab2.Get(regions2.at(r),datasets.at(5))-tab.Get(regions.at(r),datasets.at(5));
	   theDoctor::Figure res6 = tab2.Get(regions2.at(r),datasets.at(6))-tab.Get(regions.at(r),datasets.at(6));
	  
	   //cout<<"in the compressed region "<<regions2.at(r)<<" there are "<<res3.value()<<" events, and in the compressed regions "<<regions.at(r)<<" there are "<<res.value()<<" events"<<endl;
	   cout<<"The difference between the region "<<regions2.at(r)<<" and the region "<<regions.at(r)<<" is "<<res.value()<<" for the bkg from top, "<<res1.value()<<" for the bkg from lost lepton, "<<res2.value()<<" for the bkg from W, "<<res3.value() <<" fro the bkg nunu, "<<res4.value() <<" for the 600-325, "<<res5.value() <<" for the 300-50, and "<<res6.value()<<" for the 900-50"<<endl; 
	   // l+=res3.value();
	   // m+=res.value();
	 }


  //  cout<<"with selection on compressed region "<<m<<" with selection on the others regions "<<l<<endl;
  /*
  double p=0;
  double n=0;
  double w=0;
  double q=0;
  double x=0;
  double test=0;
  double rest=0;
  for (int i=4; i<7; i++)
    {
      n=0;
      w=0;
      test=0;
      rest=0;
     for(uint32_t r=0; r<regions.size();r++)
       {
	 theDoctor::Figure res = tab.Get(regions.at(r),datasets.at(i));
	 theDoctor::Figure res2 = tab2.Get(regions2.at(r),datasets.at(i));
	 theDoctor::Figure res3 = tab3.Get(regions3.at(r),datasets.at(i));
	 theDoctor::Figure res4 = tab4.Get(regions4.at(r),datasets.at(i));
	 //	 cout<<"for the events "<<datasets.at(i)<<" in the compressed region "<<regions2.at(r)<<" there are "<<res2.value()<<" events, and in the compressed regions "<<regions.at(r)<<" there are "<<res.value()<<" events and in the regions "<<regions4.at(r)<<" there are "<<res4.value()<<" events, with a total "<<res3.value()<<" in the region"<<regions3.at(r)<<endl;
	 //	 cout<<"diff ="<<res3.value()-res.value()-res2.value()-res4.value()<<endl;
	 n+=res2.value();
	 w+=res.value();
	 test+=res3.value();
	 rest+=res4.value();
       }
     //   cout<<"for the event "<<datasets.at(i)<<" we have "<<w<<" leptons events that come from chargino and "<<n<<" leptons event that come from the top and "<<rest<<" lepton event that come from neither of the two for a total of "<<test<<endl;
     //   cout<<"total diff ="<<test-n-w-rest<<endl;
     //q+=n;
     //x+=w;
     //p+=n-w;*/
     //cout<<n-w<<" with selection on compressed region "<<w<<" with selection on the other regions "<<n<<endl;
  // }
 
 // cout<<"p : "<<p<<"   l-q : "<<l-q<<"  m-x : "<<m-x<<endl;






  /* double t=0;
 double e=0;



 for(uint32_t r=0; r<regions.size()-4;r++)
	{
	  theDoctor::Figure res = tab.Get(regions.at(r),datasets.at(8));
	  theDoctor::Figure res3 = tab2.Get(regions2.at(r),datasets.at(8));
	  //	  cout<<"in the compressed region "<<regions2.at(r)<<" there are "<<res3.value()<<" events, and in the compressed regions "<<regions.at(r)<<" there are "<<res.value()<<" events"<<endl;
	  t+=res3.value();
	  e+=res.value();
	}
 //cout <<t<<" "<<e<<endl;
 /*
 for (int i=0;i<5;i++)
   {
     bool k=(i==3);
     cout<<k<<edl;
   }
 */

}
