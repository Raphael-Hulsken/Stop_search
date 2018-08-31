#include <exception>
#include <iostream>
#include "stdlib.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TStyle.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"
#include "interface/Figure.h"

using namespace std;
using namespace theDoctor;


int main()
{
  ofstream fichier1("result_combine.dat");
  int y =0;
  fichier1<<"card & 300-50 & 600-325 & 600-350 & 900-50 \\"<<endl;
  string commande ="cd ./cards";
  string commande2;
  string commande3; 


  fichier1<<"region0";

  for (int j=0; j<4 ;j++)
    {
      
      
      if (j==0)
	{
	  commande3="combine -M ProfileLikelihood --significance cards/300_50.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	}
      
      if (j==1)
	{
	  commande3="combine -M ProfileLikelihood --significance cards/600_325.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	}
      
      
      if (j==2)
	{
	  commande3="combine -M ProfileLikelihood --significance cards/600_350.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	}
      
      if (j==3)
	{
	  commande3="combine -M ProfileLikelihood --significance cards/900_50.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	}
      
      std::system(commande3.c_str());
      
      // cout <<j<<endl;
      ifstream fichier("essaie.txt");
      string ligne;
      cout<<commande3<<endl;
      
      while (getline(fichier, ligne))
	{
	  y++;
	  if (y==3)
	    {
	      size_t found = ligne.find(" ");
	      string val_neutr = ligne.substr(found+1);
	      fichier1<<" & "<<val_neutr;
	      cout<<val_neutr<<endl;
	    }
	}
      y=0;
      fichier.close();
    }

  fichier1<<" \\"<<endl;

  for (int i=1; i<18 ; i++)
    {
      string str = std::to_string(i);
      fichier1<<"region"<<i;
      /* commande2=commande+str+"/";
	 std::system(commande2.c_str());*/
      //cout<<*foo<<endl;
      cout<<commande2.c_str()<<endl;
      for (int j=0; j<4 ;j++)
	{
	 
	  
	  if (j==0)
	    {
	      commande3="combine -M ProfileLikelihood --significance cards"+str+"/300_50.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	    }
	  
	  if (j==1)
	    {
	       commande3="combine -M ProfileLikelihood --significance cards"+str+"/600_325.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	    }


	  if (j==2)
	    {
	       commande3="combine -M ProfileLikelihood --significance cards"+str+"/600_350.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	    }

	  if (j==3)
	    {
	        commande3="combine -M ProfileLikelihood --significance cards"+str+"/900_50.tab -t -1 --expectSignal=1 --toysFreq -n 1_1 >> essaie.txt";
	    }

	  std::system(commande3.c_str());
      
	  // cout <<j<<endl;
	  ifstream fichier("essaie.txt");
	   string ligne;
	  cout<<commande3<<endl;
	  
	   while (getline(fichier, ligne))
	    {
	      y++;
	      if (y==3)
		{
		  size_t found = ligne.find(" ");
		  string val_neutr = ligne.substr(found+1);
		  fichier1<<" & "<<val_neutr;
		  cout<<val_neutr<<endl;
		}
	    }
	   y=0;

	   fichier.close();
	  
	}
      fichier1<<" \\"<<endl;
      
      
      }

  fichier1.close();
  
  
}
