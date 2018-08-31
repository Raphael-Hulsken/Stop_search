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
double Zbi(double n_sig, double n_b, double rel_uncert_b = 0.3 )
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

//estimation of the significance for discovery 
double SDisc(double s, double b, double error=0.03)
{
  return s/sqrt(b+error*error*b*b);
}

//estimation of the significance for exclusion
double SExclu(double s, double b, double error=0.3)
{
  return s/sqrt(s+b+error*error*b*b);
}


//usefull for the computation of the Zbi error 
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



//--------------------
//Function to create the table in latex form, I write it in order to make it on multiple line othewise we can not read it but it should be change for less number of regions
//---------------------
void PrintLatex(vector<string> regions, vector<string> datasetLeg ,vector<TH1F*> histo,std::ostream& output,vector<double> table)
{
	output << left;
    // Begin tabular
    output << "\\documentclass{article}" << endl;
    output << "\\usepackage{graphicx}" << endl;
    output << "\\usepackage{pdflscape}" << endl;
    output << "\\usepackage{xcolor}" << endl;
    output << "\\begin{document}" << endl;
    output << "\\begin{landscape}" << endl;
    output<<"\\LARGE"<<endl;

    output << "\\noindent\\hrulefill" << endl;
    output << "\\smallskip\\noindent" << endl;




    //////////////Separation that could be used if the table is to big
        output << "\\resizebox{\\linewidth}{!}{%" << endl;
    output << "\\begin{tabular}{|l|";
    for (unsigned int i = 0 ; i < 10 ; i++) output << "c";
    output << "|}" << endl;

    // Line before header
	output << "\\hline" << endl;

	// Header
	output << "&" << endl;
	for (unsigned int i = 0 ; i < 10 ; i++)
	{
		output << "\\textbf{" << regions[i] << "} ";
		if (i < 10 - 1) output << "\t&";
        else output << "\t\\\\";
	    output << endl;
	}
	
    // Line after header
	output << "\\hline" << endl;

	// Rows
	for (unsigned int i = 4 ; i < 7 ; i++)
	{
		output << "\\textbf{$" << datasetLeg[i] << "$} ";
		output << "\t & ";
		for (unsigned int j = 0 ; j < 10 ; j++)
		{
		  //in order to co;pute the error
		  double Zbi1;
		  double Zbi2;
		  Zbi1=sqrt(pow(Zbi(histo.at(i)->GetBinContent(j+1)+histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j]),2)+pow(Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)-histo.at(8)->GetBinError(j+1), table[j]),2));
		  Zbi2=sqrt(pow(Zbi(histo.at(i)->GetBinContent(j+1)-histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j]),2)+pow(Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)+histo.at(8)->GetBinError(j+1), table[j]),2));
		  //to choose the important results
		  if (Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])>=0.5)//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)), table[j])>=0.5)
		    {
		      output <<"\\textcolor{green}{"<< Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j]) << " $\\pm$ "<<max(abs(Zbi1 -Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])),abs(Zbi2- Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])))<<"}";//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)), table[j]) << "}";// and "<<SDisc(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1))<<" and "<<SExclu(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)) <<"  ";
		    }
		  else 
		    {
		      output << Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j]) << " $\\pm$ "<<max(abs(Zbi1 -Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])),abs(Zbi2- Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])))<<"";//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)),table[j]) << "";// and "<<SDisc(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1))<<" and "<<SExclu(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)) <<" ";
		    }
		  if (j < 10 - 1) output << "\t & ";
		  else output << "\t \\\\";
		}
		output << endl;
	}
	
    // Line after content
    output << "\\hline" << endl;
    
    // End tabular
    output << "\\end{tabular}}\\\\ " << endl;






output << "\\resizebox{\\linewidth}{!}{%" << endl;
 output << "\\begin{tabular}{|l|";
 for (unsigned int i = 10 ; i < 20 ; i++) output << "c";
    output << "|}" << endl;

    // Line before header
	output << "\\hline" << endl;

	// Header
	output << "&" << endl;
	for (unsigned int i = 10 ; i < 20 ; i++)
	{
		output << "\\textbf{" << regions[i] << "} ";
		if (i < 20 - 1) output << "\t&";
        else output << "\t\\\\";
	    output << endl;
	}
	
    // Line after header
	output << "\\hline" << endl;

	// Rows
	for (unsigned int i = 4 ; i < 7 ; i++)
	{
		output << "\\textbf{$" << datasetLeg[i] << "$} ";
		output << "\t & ";
		for (unsigned int j = 10 ; j < 20 ; j++)
		{
		  //in order to co;pute the error
		  double Zbi1;
		  double Zbi2;
		  Zbi1=sqrt(pow(Zbi(histo.at(i)->GetBinContent(j+1)+histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j]),2)+pow(Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)-histo.at(8)->GetBinError(j+1), table[j]),2));
		  Zbi2=sqrt(pow(Zbi(histo.at(i)->GetBinContent(j+1)-histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j]),2)+pow(Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)+histo.at(8)->GetBinError(j+1), table[j]),2));
		  //to choose the important results
		  if (Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])>=0.5)//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)), table[j])>=0.5)
		    {
		      output <<"\\textcolor{green}{"<< Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j]) << " $\\pm$ "<<max(abs(Zbi1 -Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])),abs(Zbi2- Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])))<<"}";//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)), table[j]) << "}";// and "<<SDisc(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1))<<" and "<<SExclu(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)) <<"  ";
		    }
		  else 
		    {
		      output << Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j]) << " $\\pm$ "<<max(abs(Zbi1 -Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])),abs(Zbi2- Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])))<<"";//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)),table[j]) << "";// and "<<SDisc(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1))<<" and "<<SExclu(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)) <<" ";
		    }
		  if (j < 20 - 1) output << "\t & ";
            else output << "\t \\\\";
		}
		output << endl;
	}
	
    // Line after content
    output << "\\hline" << endl;
    
    // End tabular
    output << "\\end{tabular}}\\\\ " << endl;





    




output << "\\resizebox{\\linewidth}{!}{%" << endl;
 output << "\\begin{tabular}{|l|";
 for (unsigned int i = 20 ; i < regions.size() ; i++) output << "c";
    output << "|}" << endl;

    // Line before header
	output << "\\hline" << endl;

	// Header
	output << "&" << endl;
	for (unsigned int i = 20 ; i < regions.size() ; i++)
	{
		output << "\\textbf{" << regions[i] << "} ";
		if (i < regions.size() - 1) output << "\t&";
        else output << "\t\\\\";
	    output << endl;
	}
	
    // Line after header
	output << "\\hline" << endl;

	// Rows
	for (unsigned int i = 4 ; i < 7 ; i++)
	{
		output << "\\textbf{$" << datasetLeg[i] << "$} ";
		output << "\t & ";
		for (unsigned int j = 20 ; j < regions.size() ; j++)
		{
		  //in order to co;pute the error
		  double Zbi1;
		  double Zbi2;
		  Zbi1=sqrt(pow(Zbi(histo.at(i)->GetBinContent(j+1)+histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j]),2)+pow(Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)-histo.at(8)->GetBinError(j+1), table[j]),2));
		  Zbi2=sqrt(pow(Zbi(histo.at(i)->GetBinContent(j+1)-histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j]),2)+pow(Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)+histo.at(8)->GetBinError(j+1), table[j]),2));
		  //to choose the important results
		  if (Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])>=0.5)//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)), table[j])>=0.5)
		    {
		      output <<"\\textcolor{green}{"<< Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j]) << " $\\pm$ "<<max(abs(Zbi1 -Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])),abs(Zbi2- Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])))<<"}";//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)), table[j]) << "}";// and "<<SDisc(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1))<<" and "<<SExclu(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)) <<"  ";
		    }
		  else 
		    {
		      output << Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j]) << " $\\pm$ "<<max(abs(Zbi1 -Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])),abs(Zbi2- Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])))<<"";//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)),table[j]) << "";// and "<<SDisc(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1))<<" and "<<SExclu(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)) <<" ";
		    }
		  if (j < regions.size() - 1) output << "\t & ";
            else output << "\t \\\\";
		}
		output << endl;
	}
	
    // Line after content
    output << "\\hline" << endl;
    
    // End tabular
    output << "\\end{tabular}}" << endl;







    output << "\\end{landscape}" << endl;
    output << "\\end{document}" << endl;


}
















































int main(int argc, char *argv[]){


  //to add the error
  vector<double> table;
    string ligne;
    /* ifstream fichier("file_error.dat");
  while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"

        {

	  table.push_back(stod(ligne));  // on l'affiche

	  }*/
  
  /*for (int i=0; i<table.size();i++)
    {cout<<table[i]<<endl;}
  */
 
  THStack *hs =new THStack("hs","Stacked 1D");
        gStyle->SetOptStat(0);
        gStyle->SetLegendBorderSize(0);
        gROOT->ForceStyle();

	//need five arguments, the way it is writen right know it needs 3 table for data and 2 table for the name of the regions so it could work for differents set of regions or for differents scenario depending on the table you put in
        if(argc != 6)
            throw std::runtime_error("Bad number of arguments!");

        string inputTab = argv[1];
	string inputTab2 = argv[2];
	string inputTab3 = argv[3];
        TString inputFile = argv[4];
		TString inputFile2 = argv[5];

        vector<string> regions;
 vector<string> regions2;

 //here we choose the signal that we want to annalyse
        //vector<string> datasets = {"bkgOneLepFromTop","bkgLostLepton","bkgOneLepFromW","bkgZnunu", "(600,400)","(300,100)","(900,300)","data1","totalSM"};
	vector<string> datasets = {"(500,300)","(300,100)","(500,100)","data1","totalSM"};
//,"(800,500)","(1200,200)","(200,25)","data1","totalSM"};
//"(900,300)", "(400,300)", "(1000,100)"}; //@MJ@ TODO change the datasets to meaningful ones


//name for the root plot
        //vector<string> datasetsLeg = {"bkgOneLepFromTop","bkgLostLepton","bkgOneLepFromW","bkgZnunu","#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(600,400)", "#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(300,100)", "#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(800,600)","data","totalSM"};
	vector<string> datasetsLeg = {"#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(500,300)", "#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(300,100)", "#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(900,300)","data","totalSM"};
				      // "#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(800,500)", "#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(1200,200)", "#tilde{t}#rightarrowb#tilde{#chi}^{#pm}_{1}(200,25)","data","totalSM"}; //@MJ@ TODO change the datasets to meaningful ones

	//name for the latex table 
	vector<string> datasetsLeg2 = {"bkgOneLepFromTop","bkgLostLepton","bkgOneLepFromW","bkgZnunu","\\tilde{t} \\rightarrow b \\tilde{\\chi}^{\\pm}_{1}(600,325)",              "\\tilde{t}  \\rightarrow b \\tilde{\\chi}^{\\pm}_{1}(300,50)", "\\tilde{t} \\rightarrow b \\tilde{\\chi}^{\\pm}_{1}(900,50)","data","totalSM"};
				       // "\\tilde{t} \\rightarrow b \\tilde{\\chi}^{\\pm}_{1}(800,500)", "\\tilde{t} \\rightarrow b \\tilde{\\chi}^{\\pm}_{1}(1200,200)", "\\tilde{t} \\rightarrow b \\tilde{\\chi}^{\\pm}_{1}(200,25)","data","totalSM"};


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

	/*    const std::string uu = "SR1l_";
        const std::string v ="";

        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( uu, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, uu.size(), v );
                n += v.size();
            }
        }

        const std::string ww = "_";
        const std::string x = ";";

        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( ww, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, ww.size(), x );
                n += x.size();
            }
	    }*/

        /*for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            regLabels.at(reg) = regLabels.at(reg)+")";
        }*/


       TLegend *leg = new TLegend(0.6,0.65,0.85,0.85);

       //theDoctor::SonicScrewdriver sonic;
          Table tab(inputTab);
	   Table tab2(inputTab2);
	  Table tab3(inputTab3);


	  //implementation of the histogram with the data of the table
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

           vector<TH1F*> histo2;
       for(uint32_t h =0; h<datasets.size(); h++)
	 {
           histo2.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), regions.size(), 0, regions.size()));
	 }
       for(uint32_t r=0; r<regions2.size();r++)
	 {
           for(uint32_t d=0; d<datasets.size();d++)
	     {
               theDoctor::Figure res = tab2.Get(regions2.at(r),datasets.at(d));
               histo2.at(d)->SetBinContent(r+1,res.value());
               histo2.at(d)->SetBinError(r+1,res.error());
               histo2.at(d)->GetXaxis()->SetBinLabel(r+1,regLabels.at(r).c_str());
	       
	     }
	 }
       

           vector<TH1F*> histo3;
       for(uint32_t h =0; h<datasets.size(); h++)
	 {
           histo3.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), regions.size(), 0, regions.size()));
	 }
       for(uint32_t r=0; r<regions.size();r++)
	 {
           for(uint32_t d=0; d<datasets.size();d++)
	     {
               theDoctor::Figure res = tab3.Get(regions.at(r),datasets.at(d));
               histo3.at(d)->SetBinContent(r+1,res.value());
               histo3.at(d)->SetBinError(r+1,res.error());
               histo3.at(d)->GetXaxis()->SetBinLabel(r+1,regLabels.at(r).c_str());
	       
	     }
	     }



       /*for(uint32_t r=0; r<regions.size();r++)
       {
	 table.push_back(histo.at(8)->GetBinError(r+1)/histo.at(8)->GetBinContent(r+1));
	 }
       */





       //here we plot the number of event per regions
       vector<uint16_t> colors = {810, 603, 616};//, 400,632,860,880,1111};
      if(colors.size() < (histo.size()-4))
          throw std::runtime_error("Not enough colors for all histograms");

       TCanvas *can = new TCanvas("can","can");
       TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.0,1.0,1.0);
       // TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2);
         pad1->Draw();
	 // pad2->Draw();
       
        //can->SetLogy();

       pad1->cd();
       //pad1->GetFrame()->SetFillColor(0);
       pad1->SetFillStyle(0);
       //pad1->GetFrame()->SetFillStyle(0);
       pad1->SetLogy();

       histo.at(0)->SetTitle("");
       //histo.at(c)->SetLineColor(kGray+2);
       histo.at(0)->SetLineColor(kOrange-3);
       histo.at(0)->SetLineWidth(3);
       //  histo.at(0)->SetFillColor(kOrange-3);
       // histo.at(0)->SetFillStyle(1001);
       // hs->Add(histo.at(0));
       histo.at(0)->LabelsOption("v");
         histo.at(0)->Draw("e");

       histo2.at(0)->SetLineColor(kBlue-3);
       histo2.at(0)->SetLineWidth(3);
           histo2.at(0)->Draw("same e");

       histo3.at(0)->SetLineColor(kGreen-3);
       histo3.at(0)->SetLineWidth(3);
                histo3.at(0)->Draw("same e");
       
//float max = histo.at(c)->GetMaximum();
       //histo.at(c)->SetMaximum(1.5*max);
       //histo.at(c)->GetXaxis()->LabelsOption("v");
       //histo.at(c)->GetXaxis()->SetTitleOffset(1.2);
       // histo.at(c)->Draw("hist e"); //SET titles (per bin?!)
       histo.at(0)->GetXaxis()->SetTitle("");
       histo.at(0)->GetYaxis()->SetTitle("Events");


       /*histo2.at(0)->SetLineColor(kBlue-3);
       histo2.at(0)->SetLineWidth(3);
       histo2.at(0)->Draw("same e");
       */for(uint32_t r=0; r<regions.size();r++)
  {
    cout<<histo2.at(0)->GetBinContent(r)<<"    "<<histo.at(0)->GetBinContent(r)<<endl;
  } 
  auto legend1 = new TLegend(0.1,0.8,0.2,0.9);
  //  legend1->SetHeader("Masse 500/300","C"); // option "C" allows to center the header
     legend1->AddEntry(histo.at(0), "R0 T2bW");
	        legend1->AddEntry(histo2.at(0), "R0 T2tb");
			legend1->AddEntry(histo3.at(0), "R0 T2tt");
			//legend->AddEntry(ratio10, "mixt");
  
			 legend1->Draw();
       


			 //here if you want to plot the background
  /*  for(uint32_t c=0; c<histo.size(); c++)
       {
           if(c==0)
           {
              histo.at(c)->SetTitle("");
              //histo.at(c)->SetLineColor(kGray+2);
	       histo.at(c)->SetLineColor(kOrange-3);
              histo.at(c)->SetLineWidth(3);
              histo.at(c)->SetFillColor(kOrange-3);
              histo.at(c)->SetFillStyle(1001);
	      hs->Add(histo.at(c));
              //float max = histo.at(c)->GetMaximum();
              //histo.at(c)->SetMaximum(1.5*max);
              //histo.at(c)->GetXaxis()->LabelsOption("v");
              //histo.at(c)->GetXaxis()->SetTitleOffset(1.2);
	      // histo.at(c)->Draw("hist e"); //SET titles (per bin?!)
              histo.at(c)->GetXaxis()->SetTitle("");
              histo.at(c)->GetYaxis()->SetTitle("Events");
           }
	    else if (c==1)
	     {

	       //histo.at(c)->SetLineColor(kGray+2);
	       histo.at(c)->SetLineColor(kCyan+3);
	       histo.at(c)->SetLineWidth(3);
	       histo.at(c)->SetFillColor(kCyan+3);
	       histo.at(c)->SetFillStyle(1001);
	       //float max = histo.at(c)->GetMaximum();
	       //histo.at(c)->SetMaximum(1.5*max);
	       //histo.at(c)->GetXaxis()->LabelsOption("v");
	       //histo.at(c)->GetXaxis()->SetTitleOffset(1.2);
	       //histo.at(c)->Draw("hist e"); //SET titles (per bin?!)
	       //histo.at(c)->GetXaxis()->SetTitle("");
	       //histo.at(c)->GetYaxis()->SetTitle("Events");
	       hs->Add(histo.at(c));
	       //hs->Draw("");
	       }
	    else if (c==2)
	     {

	       //histo.at(c)->SetLineColor(kGray+2);
	       histo.at(c)->SetLineColor(kGreen);
	       histo.at(c)->SetLineWidth(3);
	       histo.at(c)->SetFillColor(kGreen);
	       histo.at(c)->SetFillStyle(1001);
	       //float max = histo.at(c)->GetMaximum();
	       //histo.at(c)->SetMaximum(1.5*max);
	       //histo.at(c)->GetXaxis()->LabelsOption("v");
	       //histo.at(c)->GetXaxis()->SetTitleOffset(1.2);
	       // histo.at(c)->Draw("hist e"); //SET titles (per bin?!)
	       //histo.at(c)->GetXaxis()->SetTitle("");
	       //histo.at(c)->GetYaxis()->SetTitle("Events");
	       hs->Add(histo.at(c));
	       //hs->Draw("");
	     }
	    else if (c==3)
	     {

	       //histo.at(c)->SetLineColor(kGray+2);
	       histo.at(c)->SetLineColor(kPink);
	       histo.at(c)->SetLineWidth(3);
	       histo.at(c)->SetFillColor(kPink);
	       histo.at(c)->SetFillStyle(1001);
	       //float max = histo.at(c)->GetMaximum();
	       //histo.at(c)->SetMaximum(1.5*max);
	       //histo.at(c)->GetXaxis()->LabelsOption("v");
	       //histo.at(c)->GetXaxis()->SetTitleOffset(1.2);
	       // histo.at(c)->Draw("hist e"); //SET titles (per bin?!)
	       //histo.at(c)->GetXaxis()->SetTitle("");
	       //histo.at(c)->GetYaxis()->SetTitle("Events");
	       hs->Add(histo.at(c));
	       hs->Draw("hist e");
	       }
	     else if (c==8)
	      {
		//histo.at(c)->SetLineColor(kGray+2);
		//histo.at(c)->SetLineColor(kPink);
		//histo.at(c)->SetLineWidth(3);
		//histo.at(c)->SetFillColor(kPink);
		//histo.at(c)->SetFillStyle(1001);
	      }
           else
           {               
              //histo.at(c)->SetLineColor(c+2);
              histo.at(c)->SetLineColor(colors.at(c-4));
              histo.at(c)->SetLineWidth(3);
              histo.at(c)->SetLineStyle(2);
		if (c==7)
	      {
		histo.at(c)->Draw("same e");
	      }
	      else 
		{
		  histo.at(c)->Draw("same hist e");
		}
		}
              leg->AddEntry(histo.at(c), (TString)datasetsLeg.at(c));
       }*/
      
// leg->Draw("l");

	      // pad2->cd();
              TH1F* ratio = (TH1F*) histo.at(0)->Clone();
	      //  TH1F* ratio2 = (TH1F*) histo.at(0)->Clone();
	      /*   for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
              {
		ratio->SetBinContent(b+1,histo.at(7)->GetBinContent(b+1)/(histo.at(8)->GetBinContent(b+1)));// );
              }
	      for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
	      {
		ratio2->SetBinContent(b+1,1);// );
	      }
              ratio->SetFillColor(0);
              ratio->SetLineColor(kBlack);
              ratio->SetMinimum(0);
              ratio->SetMaximum(2);
              ratio->GetXaxis()->SetTitle("MET");
              ratio->GetYaxis()->SetTitle("data/MC");
              ratio->Draw("e");
	      ratio2->SetFillColor(0);
              ratio2->SetLineColor(kBlack);
              ratio2->SetMinimum(0);
              ratio2->SetMaximum(2);
              ratio2->GetXaxis()->SetTitle("MET");
              ratio2->GetYaxis()->SetTitle("data/MC");
              ratio2->Draw("same hist");
	
	      */

       /*cout << "chi2 test 1: " << histo.at(0)->Chi2Test(histo.at(1), "WW P") << endl;

       TH1D* h1 = new TH1D("ttZ1", "ttZ1", regions.size()-6, 0, regions.size()-6);
       TH1D* h2 = new TH1D("ttZ2", "ttZ2", regions.size()-6, 0, regions.size()-6);

       uint32_t bin = 1;
       for(uint32_t b = 0; b<regions.size(); b++)
       {
           if(b+1==6 || b+1==7 || b+1==12 || b+1==15 || b+1==16 || b+1==27)
               continue;
           double value1 = histo.at(0)->GetBinContent(b+1); 
           double error1 = histo.at(0)->GetBinError(b+1); 
           h1->SetBinContent(bin,value1);
           h1->SetBinError(bin,error1);
           double value2 = histo.at(1)->GetBinContent(b+1); 
           double error2 = histo.at(1)->GetBinError(b+1); 
           h2->SetBinContent(bin,value2);
           h2->SetBinError(bin,error2);
           bin++;
       }

       cout << "chi2 test 2: " << h1->Chi2Test(h2, "WW P") << endl;*/
	      //  can->Update();
       can->SaveAs("plotTest.root"); //@MJ@ TODO svae to root file




       //then I compute the significance with Zbi per regions the both number of event and Zbi can be put in the same plot by changing the definition of the canvac and the histogram

          TCanvas *can2 = new TCanvas("can2","can2");
       TPad *pad3 = new TPad("pad3", "The pad 30% of the height",0.0,0.0,1.0,1.0);
       pad3->SetFillColor(0);
       pad3->SetFillStyle(4000);
       pad3->SetFrameFillStyle(0);/*       TPad *pad4 = new TPad("pad4", "The pad 30% of the height",0.0,0.0,1.0,0.3);
       TPad *pad5 = new TPad("pad5", "The pad 30% of the height",0.0,0.6,1.0,0.9);
       */
          pad3->Draw();
       /* pad4->Draw();
       pad5->Draw();
       */
       	      double Zbi1;
	      double Zbi2;
	        pad3->cd();
	TH1F* ratio3 = (TH1F*) histo.at(0)->Clone();
	TH1F* ratio5 = (TH1F*) histo.at(0)->Clone();
	TH1F* ratio7 = (TH1F*) histo.at(0)->Clone();
	TH1F* ratio9 = (TH1F*) histo.at(0)->Clone();
	TH1F* ratio10 = (TH1F*) histo.at(0)->Clone();

 	       //	TH1F* ratio4 = (TH1F*) histo.at(0)->Clone();
              for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
              {
		ratio3->SetBinContent(b+1,Zbi(histo.at(0)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2));//(histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)), table[b]));// );
			Zbi1=sqrt(pow(Zbi(histo.at(0)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(0)->GetBinContent(b+1)+histo.at(0)->GetBinError(b+1),histo.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo.at(0)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(0)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1)-histo.at(4)->GetBinError(b+1), 0.2),2));
			Zbi2=sqrt(pow(Zbi(histo.at(0)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(0)->GetBinContent(b+1)-histo.at(0)->GetBinError(b+1),histo.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo.at(0)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(0)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1)+histo.at(4)->GetBinError(b+1), 0.2),2));
			ratio3->SetBinError (b+1,max(abs(Zbi1),abs(Zbi2)));
			/*  ratio5->SetBinContent(b+1,Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2));//,(histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)), table[b]));// );
	      	Zbi1=sqrt(pow(Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(1)->GetBinContent(b+1)+histo.at(1)->GetBinError(b+1),histo.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1)-histo.at(4)->GetBinError(b+1), 0.2),2));

	      	Zbi2=sqrt(pow(Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(1)->GetBinContent(b+1)-histo.at(1)->GetBinError(b+1),histo.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1)+histo.at(4)->GetBinError(b+1), 0.2),2));

		ratio5->SetBinError (b+1,max(abs(Zbi1),abs(Zbi2)));
	    	ratio7->SetBinContent(b+1,Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2));//,(histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)), table[b]));// );
			Zbi1=sqrt(pow(Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(2)->GetBinContent(b+1)+histo.at(2)->GetBinError(b+1),histo.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1)-histo.at(4)->GetBinError(b+1), 0.2),2));
			Zbi2=sqrt(pow(Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(2)->GetBinContent(b+1)-histo.at(2)->GetBinError(b+1),histo.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1)+histo.at(4)->GetBinError(b+1), 0.2),2));
			ratio7->SetBinError (b+1,max(abs(Zbi1),abs(Zbi2)));*/
	     


				ratio9->SetBinContent(b+1,Zbi(histo2.at(0)->GetBinContent(b+1),histo2.at(4)->GetBinContent(b+1), 0.2));//,(histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)), table[b]));// );
	      	Zbi1=sqrt(pow(Zbi(histo2.at(0)->GetBinContent(b+1),histo2.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo2.at(0)->GetBinContent(b+1)+histo2.at(0)->GetBinError(b+1),histo2.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo2.at(0)->GetBinContent(b+1),histo2.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo2.at(0)->GetBinContent(b+1),histo2.at(4)->GetBinContent(b+1)-histo2.at(4)->GetBinError(b+1), 0.2),2));

	      	Zbi2=sqrt(pow(Zbi(histo2.at(0)->GetBinContent(b+1),histo2.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo2.at(0)->GetBinContent(b+1)-histo2.at(0)->GetBinError(b+1),histo2.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo2.at(0)->GetBinContent(b+1),histo2.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo2.at(0)->GetBinContent(b+1),histo2.at(4)->GetBinContent(b+1)+histo2.at(4)->GetBinError(b+1), 0.2),2));

		ratio9->SetBinError (b+1,max(abs(Zbi1),abs(Zbi2)));
	    	

			ratio10->SetBinContent(b+1,Zbi(histo3.at(0)->GetBinContent(b+1),histo3.at(4)->GetBinContent(b+1), 0.2));//,(histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)), table[b]));// );
			Zbi1=sqrt(pow(Zbi(histo3.at(0)->GetBinContent(b+1),histo3.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo3.at(0)->GetBinContent(b+1)+histo3.at(0)->GetBinError(b+1),histo3.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo3.at(0)->GetBinContent(b+1),histo3.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo3.at(0)->GetBinContent(b+1),histo3.at(4)->GetBinContent(b+1)-histo3.at(4)->GetBinError(b+1), 0.2),2));
			Zbi2=sqrt(pow(Zbi(histo3.at(0)->GetBinContent(b+1),histo3.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo3.at(0)->GetBinContent(b+1)-histo3.at(0)->GetBinError(b+1),histo3.at(4)->GetBinContent(b+1), 0.2),2)+pow(Zbi(histo3.at(0)->GetBinContent(b+1),histo3.at(4)->GetBinContent(b+1), 0.2)-Zbi(histo3.at(0)->GetBinContent(b+1),histo3.at(4)->GetBinContent(b+1)+histo3.at(4)->GetBinError(b+1), 0.2),2));
			ratio10->SetBinError (b+1,max(abs(Zbi1),abs(Zbi2)));
	      }
	      /*for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
	      {
		ratio4->SetBinContent(b+1,1);// );
		}*/
	      /*   ratio5->SetFillColor(0);
              ratio5->SetLineColor(kBlack);
              ratio5->SetMinimum(0);
              ratio5->SetMaximum(ratio5->GetBinContent(ratio5->GetMaximumBin())+ratio5->GetBinError(ratio5->GetMaximumBin())+0.1);
              ratio5->GetXaxis()->SetTitle("MET");
              ratio5->GetYaxis()->SetTitle("Zbi");
              ratio5->Draw("e");*/
	      /*   ratio7->SetFillColor(0);
              ratio7->SetLineColor(kRed);
              ratio7->SetMinimum(0);
              ratio7->SetMaximum(ratio7->GetBinContent(ratio7->GetMaximumBin())+ratio7->GetBinError(ratio7->GetMaximumBin())+0.1);
              ratio7->GetXaxis()->SetTitle("MET");
              ratio7->GetYaxis()->SetTitle("Zbi");
	      ratio7->Draw("e");*/


	      ratio3->SetFillColor(0);
	       ratio3->SetLineColor(kOrange-3);
	       //ratio3->SetFillStyle(3001);
	       ratio3->SetFillStyle(0);
	      // ratio3->SetFillColor(kOrange-3);
	      //     ratio3->SetLineColor(kOrange-3);
	      ratio3->SetMinimum(0);
	      ratio3->SetMaximum(ratio3->GetBinContent(ratio3->GetMaximumBin())+ratio3->GetBinError(ratio3->GetMaximumBin())+0.2);
              ratio3->GetXaxis()->SetTitle("MET");
              ratio3->GetYaxis()->SetTitle("Zbi");
	      //  ratio3->Draw("same HIST e, Y+");
	     ratio3->Draw("e");

	        ratio9->SetFillColor(0);
	      ratio9->SetLineColor(kBlue-3);
	      ratio9->SetMinimum(0);
	      // ratio9->SetMaximum(10);
	      //  ratio9->SetMaximum(ratio9->GetBinContent(ratio9->GetMaximumBin())+ratio9->GetBinError(ratio9->GetMaximumBin())+0.1);
              ratio9->GetXaxis()->SetTitle("MET");
              ratio9->GetYaxis()->SetTitle("Zbi");
	        ratio9->Draw("same e");
		

	      ratio10->SetFillColor(0);
	      ratio10->SetLineColor(kGreen-3);
	      ratio10->SetMinimum(0);
	      // ratio9->SetMaximum(10);
	      ratio10->SetMaximum(ratio9->GetBinContent(ratio9->GetMaximumBin())+ratio9->GetBinError(ratio9->GetMaximumBin())+0.1);
              ratio10->GetXaxis()->SetTitle("MET");
              ratio10->GetYaxis()->SetTitle("Zbi");
	      //   ratio10->Draw("same e");
	      

	          auto legend = new TLegend(0.1,0.7,0.2,0.9);
	      legend->SetHeader("Masse 500/300","C"); // option "C" allows to center the header
	      legend->AddEntry(ratio3, "regions I");
	      // legend->AddEntry(histo.at(0), "Nb event R0 T2tt");
	      legend->AddEntry(ratio9, "new region I with 3 and 4 Jets");
	      //   legend->AddEntry(ratio10, "R0 T2tt");

	      legend->Draw();



	      /*
	      
	      ratio10->SetFillColor(0);
	      ratio10->SetLineColor(kGreen);
	      ratio10->SetMinimum(0);
	      ratio10->SetMaximum(ratio10->GetBinContent(ratio10->GetMaximumBin())+ratio10->GetBinError(ratio10->GetMaximumBin())+0.1);
              ratio10->GetXaxis()->SetTitle("MET");
              ratio10->GetYaxis()->SetTitle("Zbi");
	      ratio10->Draw("same e");

	      auto legend = new TLegend(0.1,0.7,0.48,0.9);
	      legend->SetHeader("Masse 500/300","C"); // option "C" allows to center the header
	      legend->AddEntry(ratio3, "2Chargino");
	      legend->AddEntry(ratio9, "NoChargino");
	      legend->AddEntry(ratio10, "mixt");
  
	      legend->Draw();



	      /* ratio4->SetFillColor(0);
              ratio4->SetLineColor(kBlack);
              ratio4->SetMinimum(0);
              ratio4->SetMaximum(2);
              ratio4->GetXaxis()->SetTitle("MET");
              ratio4->GetYaxis()->SetTitle("#tilde{b}#rightarrowt#tilde{#chi}^{#pm}_{1}(600,325)/MC");
              ratio4->Draw("same hist");*/


	      // pad4->cd();
 /*
 TH1F* ratio5 = (TH1F*) histo.at(0)->Clone();
   //   TH1F* ratio6 = (TH1F*) histo.at(0)->Clone();
              for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
              {
		ratio5->SetBinContent(b+1,Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), table[b]));//,(histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)), table[b]));// );
		//	Zbi1=sqrt(pow(Zbi(histo.at(5)->GetBinContent(b+1)+histo.at(5)->GetBinError(b+1),histo.at(8)->GetBinContent(b+1), table[b]),2)+pow(Zbi(histo.at(5)->GetBinContent(b+1),histo.at(8)->GetBinContent(b+1)-histo.at(8)->GetBinError(b+1), table[b]),2));

		//	Zbi2=sqrt(pow(Zbi(histo.at(5)->GetBinContent(b+1)-histo.at(5)->GetBinError(b+1),histo.at(8)->GetBinContent(b+1), table[b]),2)+pow(Zbi(histo.at(5)->GetBinContent(b+1),histo.at(8)->GetBinContent(b+1)+histo.at(8)->GetBinError(b+1), table[b]),2));

		ratio5->SetBinError (b+1,max(abs(Zbi1 -Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), table[b])),abs(Zbi2- Zbi(histo.at(1)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), table[b]))));
	}
	      /*for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
	      {
		ratio6->SetBinContent(b+1,1);// );
		}*/
	      /*           ratio5->SetFillColor(0);
              ratio5->SetLineColor(kBlack);
              ratio5->SetMinimum(0);
              ratio5->SetMaximum(ratio5->GetBinContent(ratio5->GetMaximumBin())+ratio5->GetBinError(ratio5->GetMaximumBin())+0.1);
              ratio5->GetXaxis()->SetTitle("MET");
              ratio5->GetYaxis()->SetTitle( (TString)datasetsLeg.at(1));
              ratio5->Draw("e");
	      /*ratio6->SetFillColor(0);
              ratio6->SetLineColor(kBlack);
              ratio6->SetMinimum(0);
              ratio6->SetMaximum(2);
              ratio6->GetXaxis()->SetTitle("MET");
              ratio6->GetYaxis()->SetTitle("#tilde{b}#rightarrowt#tilde{#chi}^{#pm}_{1}(800,350)/MC");
              ratio6->Draw("same hist");*/


 // pad5->cd();
 //	      TH1F* ratio7 = (TH1F*) histo.at(0)->Clone();
 	      //TH1F* ratio8 = (TH1F*) histo.at(0)->Clone();
 //           for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
 /*           {
		ratio7->SetBinContent(b+1,Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), table[b]));//,(histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)), table[b]));// );
		//	Zbi1=sqrt(pow(Zbi(histo.at(6)->GetBinContent(b+1)+histo.at(6)->GetBinError(b+1),histo.at(8)->GetBinContent(b+1), table[b]),2)+pow(Zbi(histo.at(6)->GetBinContent(b+1),histo.at(8)->GetBinContent(b+1)-histo.at(8)->GetBinError(b+1), table[b]),2));
		//	Zbi2=sqrt(pow(Zbi(histo.at(6)->GetBinContent(b+1)-histo.at(6)->GetBinError(b+1),histo.at(8)->GetBinContent(b+1), table[b]),2)+pow(Zbi(histo.at(6)->GetBinContent(b+1),histo.at(8)->GetBinContent(b+1)+histo.at(8)->GetBinError(b+1), table[b]),2));
				ratio7->SetBinError (b+1,max(abs(Zbi1 -Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), table[b])),abs(Zbi2- Zbi(histo.at(2)->GetBinContent(b+1),histo.at(4)->GetBinContent(b+1), table[b]))));
	     }
	      /* for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
	      {
		ratio8->SetBinContent(b+1,1);// );
		}*/
 /*      ratio7->SetFillColor(0);
              ratio7->SetLineColor(kBlack);
              ratio7->SetMinimum(0);
              ratio7->SetMaximum(ratio7->GetBinContent(ratio7->GetMaximumBin())+ratio7->GetBinError(ratio7->GetMaximumBin())+0.1);
              ratio7->GetXaxis()->SetTitle("MET");
              ratio7->GetYaxis()->SetTitle( (TString)datasetsLeg.at(2));
              ratio7->Draw("e");
	      /*ratio8->SetFillColor(0);
              ratio8->SetLineColor(kBlack);
              ratio8->SetMinimum(0);
              ratio8->SetMaximum(2);
              ratio8->GetXaxis()->SetTitle("MET");
              ratio8->GetYaxis()->SetTitle("#tilde{b}#rightarrowt#tilde{#chi}^{#pm}_{1}(900,50)/MC");
              ratio8->Draw("same hist");*/



	    
       
 can2->Update();
       can2->SaveAs("plotTest2.root"); //@MJ@ TODO svae to root file
       ofstream file("table_significance.tex");



       //if you want to have a table in the latex format see the definition of the region at the beggining of the file
       // PrintLatex(regions, datasetsLeg2 ,histo,file,table);




       //some test for the Zbi value
       /*  cout<<"Zbi for H1="<<Zbi(histo.at(6)->GetBinContent(26),histo.at(8)->GetBinContent(26),table[25])<<" ad Zbi for H2="<<Zbi(histo.at(6)->GetBinContent(27),histo.at(8)->GetBinContent(27),table[26])<<" Zbi for E1="<<Zbi(histo.at(6)->GetBinContent(17),histo.at(8)->GetBinContent(17),table[16])<<" Zbi for B1="<<Zbi(histo.at(6)->GetBinContent(6),histo.at(8)->GetBinContent(6),table[4])<<endl;
       cout<<"S/sqrtB H1="<<histo.at(6)->GetBinContent(26)/sqrt(histo.at(8)->GetBinContent(26))<<" S/sqrtB H2="<<histo.at(6)->GetBinContent(27)/sqrt(histo.at(8)->GetBinContent(27))<<" S/sqrtB E1="<<histo.at(6)->GetBinContent(17)/sqrt(histo.at(8)->GetBinContent(17))<<" S/sqrtB B1="<<histo.at(6)->GetBinContent(6)/sqrt(histo.at(8)->GetBinContent(6))<<endl;
       cout<<"S/sqrtB+sigma H1="<<histo.at(6)->GetBinContent(26)/sqrt(histo.at(8)->GetBinContent(26)+pow(histo.at(8)->GetBinContent(26)*table[25],2))<<" S/sqrtB+sigma H2="<<histo.at(6)->GetBinContent(27)/sqrt(histo.at(8)->GetBinContent(27)+pow(histo.at(8)->GetBinContent(27)*table[26],2))<<" S/sqrtB+sigma E1="<<histo.at(6)->GetBinContent(17)/sqrt(histo.at(8)->GetBinContent(17)+pow(histo.at(8)->GetBinContent(17)*table[16],2))<<" S/sqrtB+sigma B1="<<histo.at(6)->GetBinContent(6)/sqrt(histo.at(8)->GetBinContent(6)+pow(histo.at(8)->GetBinContent(6)*table[5],2))<<endl;*/

       /* cout<<"Zbi for H1="<<Zbi(histo.at(6)->GetBinContent(26),histo.at(8)->GetBinContent(26),0.2)<<" ad Zbi for H2="<<Zbi(histo.at(6)->GetBinContent(27),histo.at(8)->GetBinContent(27),0.2)<<" Zbi for E1="<<Zbi(histo.at(6)->GetBinContent(17),histo.at(8)->GetBinContent(17),0.2)<<" Zbi for B1="<<Zbi(histo.at(6)->GetBinContent(6),histo.at(8)->GetBinContent(6),0.2)<<endl;
       cout<<"S/sqrtB H1="<<histo.at(6)->GetBinContent(26)/sqrt(histo.at(8)->GetBinContent(26))<<" S/sqrtB H2="<<histo.at(6)->GetBinContent(27)/sqrt(histo.at(8)->GetBinContent(27))<<" S/sqrtB E1="<<histo.at(6)->GetBinContent(17)/sqrt(histo.at(8)->GetBinContent(17))<<" S/sqrtB B1="<<histo.at(6)->GetBinContent(6)/sqrt(histo.at(8)->GetBinContent(6))<<endl;
       cout<<"S/sqrtB+sigma H1="<<histo.at(6)->GetBinContent(26)/sqrt(histo.at(8)->GetBinContent(26)+pow(histo.at(8)->GetBinContent(26)*0.2,2))<<" S/sqrtB+sigma H2="<<histo.at(6)->GetBinContent(27)/sqrt(histo.at(8)->GetBinContent(27)+pow(histo.at(8)->GetBinContent(27)*0.2,2))<<" S/sqrtB+sigma E1="<<histo.at(6)->GetBinContent(17)/sqrt(histo.at(8)->GetBinContent(17)+pow(histo.at(8)->GetBinContent(17)*0.2,2))<<" S/sqrtB+sigma B1="<<histo.at(6)->GetBinContent(6)/sqrt(histo.at(8)->GetBinContent(6)+pow(histo.at(8)->GetBinContent(6)*0.2,2))<<endl;
       */
       //  cout<<"Zbi for H1="<<Zbi(histo.at(6)->GetBinContent(26)*10,histo.at(8)->GetBinContent(26)*10,table[25])<<endl;
       // cout<<"Zbi for H1="<<Zbi(histo.at(6)->GetBinContent(26)*10,histo.at(8)->GetBinContent(26)*10,table[25])<<" ad Zbi for H2="<<Zbi(histo.at(6)->GetBinContent(27)*10,histo.at(8)->GetBinContent(27)*10,table[26])<<" Zbi for E1="<<Zbi(histo.at(6)->GetBinContent(17)*10,histo.at(8)->GetBinContent(17)*10,table[16])<<" Zbi for B1="<<Zbi(histo.at(6)->GetBinContent(6)*10,histo.at(8)->GetBinContent(6)*10,table[4])<<endl;
       
       /*for (int j =0;j<100;j++)
	 {
	   cout<<j<<" "<<Zbi(j,1000)<<" "<<SDisc(j,1000)<<" "<<SExclu(j,1000)<<endl;
	   }*/
       

       /* for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
              {
		cout <<histo.at(0)->GetBinContent(b+1)+histo.at(1)->GetBinContent(b+1)+histo.at(2)->GetBinContent(b+1)+histo.at(3)->GetBinContent(b+1)-histo.at(8)->GetBinContent(b+1)<<endl;
		}*/




       /*


	for (unsigned int i = 4 ; i < 7 ; i++)
	{
	  for (unsigned int j = 0 ; j < regions.size() ; j++)
	    {
	      double Zbi1;
	      double Zbi2;
	      Zbi1=sqrt(Zbi(histo.at(i)->GetBinContent(j+1)+histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j])*Zbi(histo.at(i)->GetBinContent(j+1)+histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j])+Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)-histo.at(8)->GetBinError(j+1), table[j])*Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)-histo.at(8)->GetBinError(j+1), table[j]));
	      Zbi2=sqrt(Zbi(histo.at(i)->GetBinContent(j+1)-histo.at(i)->GetBinError(j+1),histo.at(8)->GetBinContent(j+1), table[j])*Zbi(histo.at(i)->GetBinContent(j+1)-histo.at(i)->GetBinError(j+1),histo.at(11)->GetBinContent(j+1), table[j])+Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)+histo.at(8)->GetBinError(j+1), table[j])*Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1)+histo.at(8)->GetBinError(j+1), table[j]));
	      cout <<"Zbi pour "<<i <<", et "<<j<<" :"<< Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j]) << " Zbi1: "<<Zbi1<<" Zbi2 : "<<Zbi2<<" "<<"error : "<<max(abs(Zbi1 -Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])),abs(Zbi2- Zbi(histo.at(i)->GetBinContent(j+1),histo.at(8)->GetBinContent(j+1), table[j])))<<""<<endl;//(histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)), table[j]) << "}";// and "<<SDisc(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1))<<" and "<<SExclu(histo.at(i)->GetBinContent(j+1),histo.at(0)->GetBinContent(j+1)+histo.at(1)->GetBinContent(j+1)+histo.at(2)->GetBinContent(j+1)+histo.at(3)->GetBinContent(j+1)) <<"  ";
	    }
	    }*/
}
