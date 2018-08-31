#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
//bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

bool CheckTrigger(bool data, string dataset){
       if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");
        
        bool trigged = false;

        if(data && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 )
        {
            cout << "1 lep data event " << endl;
            if( dataset.find("data_single_electron")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl)
                trigged = true;
            else if( dataset.find("data_single_muon")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
                trigged = true;
            else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET110_MHT110 || myEvent.HLT_MET120_MHT120 || myEvent.HLT_MET) ) 
            {
                if(abs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl ) //they are not doing that!!!
                    trigged = false;
                else if(abs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
                    trigged = false;
                else
                    trigged = true;
            }
            else
            {
                cout << "for data no trigger was found 1l " << endl;
                cout << "dataset " << dataset << " myEvent.HLT_SingleEl "<< myEvent.HLT_SingleEl << " myEvent.HLT_SingleMu " <<  myEvent.HLT_SingleMu << endl;
            }
            cout << "datset " << dataset << " triggered by 1lep trigger " << trigged << endl;

        }
        else if(data && (myEvent.ngoodleps + myEvent.nvetoleps)>2 )
        {
            if( dataset.find("data_double_eg")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==11 && abs(myEvent.lep2_pdgid)==11  && myEvent.HLT_DiEl)
                trigged = true;
            else if( dataset.find("data_double_mu")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==13 &&  abs(myEvent.lep2_pdgid)==13 && myEvent.HLT_DiMu)
                trigged = true;
            else if( dataset.find("data_muon_eg")!=std::string::npos &&  abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==24 && myEvent.HLT_MuE)
                trigged = true;
            else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET110_MHT110 || myEvent.HLT_MET120_MHT120 || myEvent.HLT_MET))
            {
                if(  abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==22 && !myEvent.HLT_DiEl )
                    trigged = false;
                else if( abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==26 && !myEvent.HLT_DiMu)
                    trigged = false;
                else if(abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==24 && !myEvent.HLT_MuE)
                    trigged = false;
                else
                    trigged = true;
            }
            else
            {
            } 

        }
        else
        {
            trigged = true;
        }
        return trigged;

}

bool passFilters()
{
    bool filtRes = false;
    if(!myEvent.is_data)
        filtRes = true;
    else if(myEvent.filt_met && 
            myEvent.filt_badChargedCandidateFilter && 
            myEvent.filt_jetWithBadMuon && 
            myEvent.filt_pfovercalomet && 
            myEvent.filt_badMuonFilter && 
            !myEvent.filt_duplicatemuons && 
            !myEvent.filt_badmuons && 
            myEvent.filt_nobadmuons )
        filtRes = true;
    else
        filtRes = false ;  

    //cout << "result of filters " << filtRes << endl;
    return filtRes;
        
}

bool passGoodVtx()
{
    bool selRes = false;
    if(myEvent.nvertex>=1)
        selRes= true;

   // cout << "good vertex returned " << selRes << endl;
    return selRes;
}

bool tightCSVV2()
{
    if(myEvent.ngoodbtags ==0 )
        return true;
    /*for (uint32_t j =0 ; j<myEvent.ak4pfjets_CSV->size(); j++)
    {
        if( myEvent.ak4pfjets_CSV->at(j) >0.935)
            return true;
    }
    return false;*/
    if(myEvent.ntightbtags> 0)
        return true;

    return false;
}


float Mlb()
{

    if(myEvent.ngoodbtags==0)
    {
        cout << "csv size " << myEvent.ak4pfjets_CSV->size()  << endl;
        if(myEvent.ak4pfjets_CSV->size()==0)
            return -13;
        if(myEvent.ngoodleps==0)
            return -13;

 
         ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >   lb =  myEvent.ak4pfjets_leadbtag_p4 + myEvent.lep1_p4;
         float Mlb2 = lb.M2();
         if(Mlb2<0)
             throw std::runtime_error("Square root of Mlb is smaller than zero");
         return sqrt(Mlb2);
    }
    else
        return myEvent.Mlb;
}

bool dilepSel()
{

    if ( ( (myEvent.ngoodleps + myEvent.nvetoleps)>2) &&
         ( (abs(myEvent.lep1_pdgid)==13 && myEvent.lep1_passMediumID) || 
	   (abs(myEvent.lep1_pdgid)==11 && fabs(myEvent.lep1_eta)<1.4442 && myEvent.lep1_passMediumID ) ) &&
         ( (abs(myEvent.lep2_pdgid)==13 && myEvent.lep2_passMediumID) || 
	   (abs(myEvent.lep2_pdgid)==11 && fabs(myEvent.lep2_eta)<1.4442 && myEvent.lep2_passMediumID ) ) )
    {
        myEvent.pfmet = myEvent.pfmet_rl;
        myEvent.pfmet_phi = myEvent.pfmet_phi_rl;
        myEvent.MT2W = myEvent.MT2W_rl;
        myEvent.dphi_ak4pfjets_met = myEvent.dphi_ak4pfjets_met_rl;
        myEvent.mt_met_lep = myEvent.mt_met_lep_rl;
        myEvent.topnessMod = myEvent.topnessMod_rl;
        myEvent.lep1_dphiMET = myEvent.lep1_dphiMET_rl;
        myEvent.lep2_dphiMET = myEvent.lep2_dphiMET_rl;
        return true;
    }

    else
        return false;   
}






//---------------------
//Normal selection
//----------------------
bool baseline(){ return myEvent.pfmet>=0 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0 && passGoodVtx() && passFilters();}
bool SR1l() { return (myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && baseline()  ); }
bool SR1l_A_250lessMETless350() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_A_350lessMETless450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_A_450lessMETless600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_A_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_B_250lessMETless450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}
bool SR1l_C_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_C_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_C_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_C_550lessMETless650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l_C_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650);}
bool SR1l_D_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}
bool SR1l_E_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_E_350lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l_E_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_G_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_G_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_G_450lessMETless600() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_G_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_H_250lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_I_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}




//---------------------------
//New regions to find 900_50 2Chargino 
//---------------
//here modified middle MET for the H regions 
bool SR1l_H_250lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<550);}
bool SR1l_H_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}

bool SR1l_H_250lessMETless650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<650);}
bool SR1l_H_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=650);}

bool SR1l_H_250lessMETless500() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<500);}
bool SR1l_H_500lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=500);}

bool SR1l_H_250lessMETless600() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<600);}
bool SR1l_H_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}


//same modifications but with modified topness put at 12
bool SR1l_H_250lessMETless450_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInf_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}

bool SR1l_H_250lessMETless550_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<550);}
bool SR1l_H_550lessMETlessInf_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}

bool SR1l_H_250lessMETless650_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<650);}
bool SR1l_H_650lessMETlessInf_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=650);}

bool SR1l_H_250lessMETless500_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<500);}
bool SR1l_H_500lessMETlessInf_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=500);}

bool SR1l_H_250lessMETless600_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<600);}
bool SR1l_H_600lessMETlessInf_12() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=12 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}





//--------------------
//New regions to find 600_325 2Chargino
//------------------------
/*bool SR1l_E_250lessMETless350B() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350) && myEvent.ngoodbtags==2;}
bool SR1l_E_350lessMETless550B() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<550) && myEvent.ngoodbtags==2;}
bool SR1l_G_250lessMETless350B() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350) && myEvent.ngoodbtags==2;}
bool SR1l_G_350lessMETless450B() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450) && myEvent.ngoodbtags==2;}
bool SR1l_I_350lessMETless450B() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450) && myEvent.ngoodbtags==2;}
bool SR1l_I_450lessMETless550B() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550) && myEvent.ngoodbtags==2;}
*/




//--------------------
//New regions to find 300_50 2Chargino
//------------------------
bool SR1l_C_150lessMETless250() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=150 && myEvent.pfmet<250);}
bool SR1l_D_150lessMETless250() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=150 && myEvent.pfmet<250);}
bool SR1l_E_150lessMETless250() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=150 && myEvent.pfmet<250);}
bool SR1l_I_150lessMETless250() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=150 && myEvent.pfmet<250);}





//---------------------
//Modified regions in order to include the compressed one
//---------------------
bool SR1l_C_250lessMETless350CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_C_350lessMETless450CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_C_450lessMETless550CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<550) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true)) ;}
bool SR1l_C_550lessMETless650CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550 && myEvent.pfmet<650) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_C_650lessMETlessInfCMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_D_250lessMETless350CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<350) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_D_350lessMETless450CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=350 && myEvent.pfmet<450) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_D_450lessMETless550CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<550) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_D_550lessMETlessInfCMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_E_250lessMETless350CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_E_350lessMETless550CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<550) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_E_550lessMETlessInfCMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_F_250lessMETless450CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_F_450lessMETlessInfCMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_G_250lessMETless350CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_G_350lessMETless450CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_G_450lessMETless600CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_G_600lessMETlessInfCMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_H_250lessMETless450CMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
bool SR1l_H_450lessMETlessInfCMPR() { return ((SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450) && (myEvent.ngoodjets<5 || myEvent.lep1_pt>=150 || myEvent.lep1_dphiMET>=2 || myEvent.ak4pfjets_passMEDbtag->at(0) == true));}
		




//to compute the total I we need to find the one that appear in I that don't go into the others
bool SR1l_I_250lessMETless350_0_8() { return (SR1l() && myEvent.ngoodjets>=5 && (myEvent.dphi_ak4pfjets_met<0.8 || ( myEvent.Mlb>175 && myEvent.ntightbtags<1) ) && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450_0_8() { return (SR1l() && myEvent.ngoodjets>=5 && (myEvent.dphi_ak4pfjets_met<0.8 || ( myEvent.Mlb>175 && myEvent.ntightbtags<1) ) && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550_0_8() { return (SR1l() && myEvent.ngoodjets>=5 && (myEvent.dphi_ak4pfjets_met<0.8 || ( myEvent.Mlb>175 && myEvent.ntightbtags<1) ) && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInf_0_8() { return (SR1l() && myEvent.ngoodjets>=5 && (myEvent.dphi_ak4pfjets_met<0.8  || ( myEvent.Mlb>175 && myEvent.ntightbtags<1) )&& myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}








//----------------------
//Selection for the mother of the lepton (mother Chargino 1000024)
//-----------------------------
bool SR1lA() { return (myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && baseline() &&( myEvent.lep1_mc_motherid==1000024 || myEvent.lep1_mc_motherid==(-1000024))) ; }
bool SR1l_A_250lessMETless350A() { return (SR1lA() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_A_350lessMETless450A() { return (SR1lA() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_A_450lessMETless600A() { return (SR1lA() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_A_600lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_B_250lessMETless450A() { return (SR1lA() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600A() { return (SR1lA() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}
bool SR1l_C_250lessMETless350A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_C_350lessMETless450A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_C_450lessMETless550A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_C_550lessMETless650A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l_C_650lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650);}
bool SR1l_D_250lessMETless350A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}
bool SR1l_E_250lessMETless350A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_E_350lessMETless550A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l_E_550lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_G_250lessMETless350A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_G_350lessMETless450A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_G_450lessMETless600A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_G_600lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_H_250lessMETless450A() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_I_250lessMETless350A() { return (SR1lA() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450A() { return (SR1lA() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550A() { return (SR1lA() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInfA() { return (SR1lA() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}

















//--------------
//lepton come from top
//--------------
bool SR1lB() { return (myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && baseline() &&( myEvent.lep1_mc_motherid==24 || myEvent.lep1_mc_motherid==(-24))) ; }
bool SR1l_A_250lessMETless350B() { return (SR1lB() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_A_350lessMETless450B() { return (SR1lB() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_A_450lessMETless600B() { return (SR1lB() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_A_600lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_B_250lessMETless450B() { return (SR1lB() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600B() { return (SR1lB() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}
bool SR1l_C_250lessMETless350B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_C_350lessMETless450B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_C_450lessMETless550B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_C_550lessMETless650B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l_C_650lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650);}
bool SR1l_D_250lessMETless350B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}
bool SR1l_E_250lessMETless350B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_E_350lessMETless550B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l_E_550lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_G_250lessMETless350B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_G_350lessMETless450B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_G_450lessMETless600B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_G_600lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_H_250lessMETless450B() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_I_250lessMETless350B() { return (SR1lB() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450B() { return (SR1lB() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550B() { return (SR1lB() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInfB() { return (SR1lB() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}






//----------------
//lepton come from another one
//---------------
bool SR1lC() { return (myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && baseline() &&( myEvent.lep1_mc_motherid!=24 && myEvent.lep1_mc_motherid!=(-24) && myEvent.lep1_mc_motherid!=1000024 && myEvent.lep1_mc_motherid!=(-1000024)) ); }
bool SR1l_A_250lessMETless350C() { return (SR1lC() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_A_350lessMETless450C() { return (SR1lC() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_A_450lessMETless600C() { return (SR1lC() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_A_600lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_B_250lessMETless450C() { return (SR1lC() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600C() { return (SR1lC() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}
bool SR1l_C_250lessMETless350C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_C_350lessMETless450C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_C_450lessMETless550C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_C_550lessMETless650C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l_C_650lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650);}
bool SR1l_D_250lessMETless350C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}
bool SR1l_E_250lessMETless350C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_E_350lessMETless550C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l_E_550lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_G_250lessMETless350C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_G_350lessMETless450C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_G_450lessMETless600C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_G_600lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_H_250lessMETless450C() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_I_250lessMETless350C() { return (SR1lC() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450C() { return (SR1lC() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550C() { return (SR1lC() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInfC() { return (SR1lC() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}







//lower cut for Delta Phi ->0.5
bool SR1l_A_250lessMETless350_0_5() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_A_350lessMETless450_0_5() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_A_450lessMETless600_0_5() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_A_600lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=600);}
bool SR1l_B_250lessMETless450_0_5() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600_0_5() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}
bool SR1l_C_250lessMETless350_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_C_350lessMETless450_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_C_450lessMETless550_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_C_550lessMETless650_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l_C_650lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=650);}
bool SR1l_D_250lessMETless350_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}
bool SR1l_E_250lessMETless350_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_E_350lessMETless550_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l_E_550lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_G_250lessMETless350_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_G_350lessMETless450_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_G_450lessMETless600_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_G_600lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=600);}
bool SR1l_H_250lessMETless450_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInf_0_5() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}







//new regions with no cut on ntight
bool SR1l_B_250lessMETless450_NoNTight() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600_NoNTight() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInf_NoNTight() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_D_250lessMETless350_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 &&  myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 &&  myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInf_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 &&  myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 &&  myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInf_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 &&  myEvent.pfmet>=450);}
bool SR1l_H_250lessMETless450_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 &&  myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInf_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450);}














//lower cut in Delta Phi ->0.5 and no cut in Ntight
bool SR1l_B_250lessMETless450_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInf_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=600);}
bool SR1l_D_250lessMETless350_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 &&  myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 &&  myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInf_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 &&  myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 &&  myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInf_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 &&  myEvent.pfmet>=450);}
bool SR1l_H_250lessMETless450_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 &&  myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInf_0_5_NoNTight() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.pfmet>=450);}



//Big regions 

bool SR1l_A() {return (SR1l_A_250lessMETless350() || SR1l_A_350lessMETless450() || SR1l_A_450lessMETless600() || SR1l_A_600lessMETlessInf() );}
bool SR1l_B() {return (SR1l_B_250lessMETless450() || SR1l_B_450lessMETless600() || SR1l_B_600lessMETlessInf() );}
bool SR1l_C() {return (SR1l_C_250lessMETless350() || SR1l_C_350lessMETless450() || SR1l_C_450lessMETless550() || SR1l_C_550lessMETless650() || SR1l_C_650lessMETlessInf() );}
bool SR1l_D() {return (SR1l_D_250lessMETless350() || SR1l_D_350lessMETless450() || SR1l_D_450lessMETless550() || SR1l_D_550lessMETlessInf() );}
bool SR1l_E() {return (SR1l_E_250lessMETless350() || SR1l_E_350lessMETless550() || SR1l_E_550lessMETlessInf() );}
bool SR1l_F() {return (SR1l_F_250lessMETless450() || SR1l_F_450lessMETlessInf() );}
bool SR1l_G() {return (SR1l_G_250lessMETless350() || SR1l_G_350lessMETless450() || SR1l_G_450lessMETless600() || SR1l_G_600lessMETlessInf() );}
bool SR1l_H() {return (SR1l_H_250lessMETless450() || SR1l_H_450lessMETlessInf() );}
bool SR1l_I() {return (SR1l_I_250lessMETless350() || SR1l_I_350lessMETless450() || SR1l_I_450lessMETless550() || SR1l_I_550lessMETlessInf() );}





//---------------------------
//Change higher MET cut
//---------------------------
bool SR1l_A_450lessMETless700() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<700);}
bool SR1l_A_700lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=700);}
bool SR1l_B_450lessMETless650() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<650);}
bool SR1l_B_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=650);}
bool SR1l_C_550lessMETless750() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550 && myEvent.pfmet<750);}
bool SR1l_C_750lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=750);}
bool SR1l_D_450lessMETless700() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<700);}
bool SR1l_D_700lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=700);}
bool SR1l_E_350lessMETless650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<650);}
bool SR1l_E_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650);}
bool SR1l_F_250lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<550);}
bool SR1l_F_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}
bool SR1l_G_450lessMETless650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<650);}
bool SR1l_G_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650);}
//bool SR1l_H_250lessMETless650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<650);}
//bool SR1l_H_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=650);}
bool SR1l_I_450lessMETless700() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<700);}
bool SR1l_I_700lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=700);}


bool SR1l_less175Mlb() {return (SR1l() && myEvent.Mlb<=175  && myEvent.ngoodjets<=3 &&  myEvent.ntightbtags>=1 && myEvent.dphi_ak4pfjets_met>=0.8);}

bool SR1l_more175Mlb() {return (SR1l() && myEvent.Mlb>175  && myEvent.ngoodjets<=3 &&  myEvent.ntightbtags>=1 && myEvent.dphi_ak4pfjets_met>=0.8);}

bool SR1l_J() { return (SR1l_less175Mlb() &&  myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false);}


//bool I { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false);}

bool I_no_pt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false);}

bool I_no_dphiMET() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150  && myEvent.ak4pfjets_passMEDbtag->at(0) == false);}

bool I_no_jets() { return (SR1l() && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false);}


//bool I { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false);}


bool SR1l_B_topMod_250lessMETless450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>0  && myEvent.topnessMod<10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_topMod_450lessMETless600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>0  &&myEvent.topnessMod<10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_topMod_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>0  &&myEvent.topnessMod<10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}








bool SR1l_C_mod() {return (SR1l_C_250lessMETless350CMPR() || SR1l_C_350lessMETless450CMPR() || SR1l_C_450lessMETless550CMPR() || SR1l_C_550lessMETless650CMPR() || SR1l_C_650lessMETlessInfCMPR() );}
bool SR1l_D_mod() {return (SR1l_D_250lessMETless350CMPR() || SR1l_D_350lessMETless450CMPR() || SR1l_D_450lessMETless550CMPR() || SR1l_D_550lessMETlessInfCMPR() );}
bool SR1l_E_mod() {return (SR1l_E_250lessMETless350CMPR() || SR1l_E_350lessMETless550CMPR() || SR1l_E_550lessMETlessInfCMPR() );}
bool SR1l_F_mod() {return (SR1l_F_250lessMETless450CMPR() || SR1l_F_450lessMETlessInfCMPR() );}
bool SR1l_G_mod() {return (SR1l_G_250lessMETless350CMPR() || SR1l_G_350lessMETless450CMPR() || SR1l_G_450lessMETless600CMPR() || SR1l_G_600lessMETlessInfCMPR() );}
bool SR1l_H_mod() {return (SR1l_H_250lessMETless450CMPR() || SR1l_H_450lessMETlessInfCMPR() );}
bool SR1l_I_mod() {return (SR1l_I_250lessMETless350_0_8() || SR1l_I_350lessMETless450_0_8() || SR1l_I_450lessMETless550_0_8() || SR1l_I_550lessMETlessInf_0_8() );}





bool SR1l_I_top_mod_inf0() {return (SR1l_I_250lessMETless350() || SR1l_I_350lessMETless450() || SR1l_I_450lessMETless550() || SR1l_I_550lessMETlessInf() &&  myEvent.topnessMod<=0 );}

bool SR1l_I_top_mod_sup0() {return (SR1l_I_250lessMETless350() || SR1l_I_350lessMETless450() || SR1l_I_450lessMETless550() || SR1l_I_550lessMETlessInf() &&  myEvent.topnessMod>0 && myEvent.topnessMod<=10 );}


bool SR1l_I_Mlb_inf175() {return (SR1l_I_250lessMETless350() || SR1l_I_350lessMETless450() || SR1l_I_450lessMETless550() || SR1l_I_550lessMETlessInf() &&  myEvent.Mlb<=175 );}

bool SR1l_I_Mlb_sup175() {return (SR1l_I_250lessMETless350() || SR1l_I_350lessMETless450() || SR1l_I_450lessMETless550() || SR1l_I_550lessMETlessInf() &&   myEvent.Mlb>175 );}


bool SR1l_I_250lessMETless350_test_mother() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350 && myEvent.lep1_mc_motherid!=24 && myEvent.lep1_mc_motherid!=(-24) && myEvent.lep1_mc_motherid!=1000024 && myEvent.lep1_mc_motherid!=(-1000024));}
bool SR1l_I_350lessMETless450_test_mother() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.lep1_mc_motherid!=24 && myEvent.lep1_mc_motherid!=(-24) && myEvent.lep1_mc_motherid!=1000024 && myEvent.lep1_mc_motherid!=(-1000024));}
bool SR1l_I_450lessMETless550_test_mother() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550 && myEvent.lep1_mc_motherid!=24 && myEvent.lep1_mc_motherid!=(-24) && myEvent.lep1_mc_motherid!=1000024 && myEvent.lep1_mc_motherid!=(-1000024));}
bool SR1l_I_550lessMETlessInf_test_mother() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550 && myEvent.lep1_mc_motherid!=24 && myEvent.lep1_mc_motherid!=(-24) && myEvent.lep1_mc_motherid!=1000024 && myEvent.lep1_mc_motherid!=(-1000024));}


bool test() {return (SR1l() &&  myEvent.pfmet>=250);}




bool SR1l_I_250lessMETless350_Pt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350  && myEvent.ak4pfjets_leadbtag_p4.Pt()>300 );}
bool SR1l_I_350lessMETless450_Pt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.ak4pfjets_leadbtag_p4.Pt()>300 );}
bool SR1l_I_450lessMETless550_Pt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550  && myEvent.ak4pfjets_leadbtag_p4.Pt()>300);}
bool SR1l_I_550lessMETlessInf_Pt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550 && myEvent.ak4pfjets_leadbtag_p4.Pt()>300 );}


bool SR1l_I3_250lessMETless350_Pt() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350  && myEvent.ak4pfjets_leadbtag_p4.Pt()>150 );}
bool SR1l_I3_350lessMETless450_Pt() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.ak4pfjets_leadbtag_p4.Pt()>150 );}
bool SR1l_I3_450lessMETless550_Pt() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550  && myEvent.ak4pfjets_leadbtag_p4.Pt()>150);}
bool SR1l_I3_550lessMETlessInf_Pt() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550 && myEvent.ak4pfjets_leadbtag_p4.Pt()>300 );}







bool SR1l_all_event() {return (SR1l() && (  ( myEvent.topnessMod<10 &&  myEvent.ngoodjets<=3 ) || ( myEvent.ngoodjets>=4)));}


bool Baseline () {return (SR1l_A() || SR1l_B() || SR1l_C() || SR1l_D() || SR1l_E() || SR1l_F() || SR1l_G() || SR1l_H() || SR1l_I_mod() ); }



bool NJets3_Modtopness () {return (myEvent.ngoodjets<=3 && SR1l());}

bool NJets3_ModtopnessA () {return (myEvent.ngoodjets<=3 && SR1l() && myEvent.Mlb<=175);}

bool NJets3_ModtopnessB () {return (myEvent.ngoodjets<=3 && SR1l() && myEvent.Mlb>175);}

bool NJets3_Mlb () {return (SR1l_A() || SR1l_B());}

bool NJets3_MET () {return (myEvent.ngoodjets<=3 && SR1l());}

bool NJets4_MET () {return (myEvent.ngoodjets>3 && SR1l() && myEvent.dphi_ak4pfjets_met>=0.8);}


/*
bool baseline(){ return myEvent.pfmet>=0 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0 && passGoodVtx() && passFilters();}
bool SR1l() { return (myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && baseline()  ); }

SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false*/

bool I_no_cut_MET() {return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false); }

bool I_no_cut_dphimet() {return (myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2  && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0 && passGoodVtx() && passFilters()&& myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto  && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 &&myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250); }



bool I_no_cut_dphilep() {return (myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2  && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0 && passGoodVtx() && passFilters() && myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto  && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150  && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.dphi_ak4pfjets_met>=0.5); }



bool I_no_cut_Ptlep() {return (myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2  && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0 && passGoodVtx() && passFilters()&&  myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto  && myEvent.ngoodjets>=5 &&  myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.lep1_dphiMET<2); }

bool I_no_cut_NJets() {return (SR1l() && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250); }

bool I_no_cut_ident() {return (SR1l() && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2  && myEvent.pfmet>=250 && myEvent.ngoodjets>=5); }


bool I_cut_tight() {return (SR1l() && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2  && myEvent.pfmet>=250 && myEvent.ngoodjets>=5 && myEvent.ak4pfjets_leadJet_identification<=0.9535 ); }


bool I_cut_med() {return (SR1l() && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2  && myEvent.pfmet>=250 && myEvent.ngoodjets>=5 && myEvent.ak4pfjets_leadJet_identification<=0.8484 ); }

bool I_cut_loose() {return (SR1l() && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2  && myEvent.pfmet>=250 && myEvent.ngoodjets>=5 && myEvent.ak4pfjets_leadJet_identification<=0.5426); }
















//cut 4 jet in I
bool SR1l_I4_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I4_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I4_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I4_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}


bool SR1l_I4_250lessMETless350_no_ptCut() { return (SR1l() && myEvent.ngoodjets>=4  && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I4_350lessMETless450_no_ptCut() { return (SR1l() && myEvent.ngoodjets>=4  && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I4_450lessMETless550_no_ptCut() { return (SR1l() && myEvent.ngoodjets>=4  && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I4_550lessMETlessInf_no_ptCut() { return (SR1l() && myEvent.ngoodjets>=4  && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}



bool SR1l_I3_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I3_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I3_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I3_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}


bool SR1l_I_equal_4_no_ptCut() {return ((SR1l_I4_250lessMETless350_no_ptCut() || SR1l_I4_350lessMETless450_no_ptCut() || SR1l_I4_450lessMETless550_no_ptCut() || SR1l_I4_550lessMETlessInf_no_ptCut() ) && myEvent.ngoodjets<5 );}


bool SR1l_I_more_4() {return (SR1l_I4_250lessMETless350() || SR1l_I4_350lessMETless450() || SR1l_I4_450lessMETless550() || SR1l_I4_550lessMETlessInf() );}

bool SR1l_I_more_3() {return (SR1l_I3_250lessMETless350() || SR1l_I3_350lessMETless450() || SR1l_I3_450lessMETless550() || SR1l_I3_550lessMETlessInf() );}



bool SR1l_I_equal4_250lessMETless350() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_equal4_350lessMETless450() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_equal4_450lessMETless550() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_equal4_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}



bool SR1l_I_equal3_250lessMETless350() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_equal3_350lessMETless450() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_equal3_450lessMETless550() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_equal3_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}


bool SR1l_I_equal_34_250lessMETless350() { return (SR1l_I_equal3_250lessMETless350() || SR1l_I_equal4_250lessMETless350());}
bool SR1l_I_equal_34_350lessMETless450() { return (SR1l_I_equal3_350lessMETless450() || SR1l_I_equal4_350lessMETless450());}
bool SR1l_I_equal_34_450lessMETless550() { return (SR1l_I_equal3_450lessMETless550() || SR1l_I_equal4_450lessMETless550());}
bool SR1l_I_equal_34_550lessMETlessInf() { return (SR1l_I_equal3_550lessMETlessInf() || SR1l_I_equal4_550lessMETlessInf());}




bool SR1l_I_equal_4() {return (SR1l_I_equal4_250lessMETless350() || SR1l_I_equal4_350lessMETless450() || SR1l_I_equal4_450lessMETless550() || SR1l_I_equal4_550lessMETlessInf() );}

bool SR1l_I_equal_3() {return (SR1l_I_equal3_250lessMETless350() || SR1l_I_equal3_350lessMETless450() || SR1l_I_equal3_450lessMETless550() || SR1l_I_equal3_550lessMETlessInf() );}


bool I_equal_34() {return ( SR1l_I_equal_4() || SR1l_I_equal_3() );}





bool SR1l_I_equal4_cut_Pt_300_250lessMETless350() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal4_cut_Pt_300_350lessMETless450() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal4_cut_Pt_300_450lessMETless550() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal4_cut_Pt_300_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550 && myEvent.ak4pfjets_pt >=300);}



bool SR1l_I_equal3_cut_Pt_300_250lessMETless350() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal3_cut_Pt_300_350lessMETless450() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal3_cut_Pt_300_450lessMETless550() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal3_cut_Pt_300_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550 && myEvent.ak4pfjets_pt >=300);}


bool SR1l_I_equal_4_cut_Pt_300() {return (SR1l_I_equal4_cut_Pt_300_250lessMETless350() || SR1l_I_equal4_cut_Pt_300_350lessMETless450() || SR1l_I_equal4_cut_Pt_300_450lessMETless550() || SR1l_I_equal4_cut_Pt_300_550lessMETlessInf() );}

bool SR1l_I_equal_3_cut_Pt_300() {return (SR1l_I_equal3_cut_Pt_300_250lessMETless350() || SR1l_I_equal3_cut_Pt_300_350lessMETless450() || SR1l_I_equal3_cut_Pt_300_450lessMETless550() || SR1l_I_equal3_cut_Pt_300_550lessMETlessInf() );}


bool I_equal_34_cut_Pt_300() {return ( SR1l_I_equal_4_cut_Pt_300() || SR1l_I_equal_3_cut_Pt_300() );}


bool I_equal_34_cut_Pt_300_250_350() {return ( SR1l_I_equal4_cut_Pt_300_250lessMETless350() ||  SR1l_I_equal3_cut_Pt_300_250lessMETless350() );}
bool I_equal_34_cut_Pt_300_350_450() {return ( SR1l_I_equal4_cut_Pt_300_350lessMETless450() || SR1l_I_equal3_cut_Pt_300_350lessMETless450() );}
bool I_equal_34_cut_Pt_300_450_550() {return ( SR1l_I_equal4_cut_Pt_300_450lessMETless550() || SR1l_I_equal3_cut_Pt_300_450lessMETless550() );}
bool I_equal_34_cut_Pt_300_550_Inf() {return ( SR1l_I_equal4_cut_Pt_300_550lessMETlessInf() || SR1l_I_equal3_cut_Pt_300_550lessMETlessInf() );}



///test for cut in modified topness 
bool I_equal_34_cut_Pt_300_modtop_Inf0() {return ( (SR1l_I_equal_4_cut_Pt_300() || SR1l_I_equal_3_cut_Pt_300()) && myEvent.topnessMod<0);}


bool I_equal_34_cut_Pt_300_250_350_modtop_Inf0() {return ( (SR1l_I_equal4_cut_Pt_300_250lessMETless350() ||  SR1l_I_equal3_cut_Pt_300_250lessMETless350()) && myEvent.topnessMod<0 );}
bool I_equal_34_cut_Pt_300_350_450_modtop_Inf0() {return ( (SR1l_I_equal4_cut_Pt_300_350lessMETless450() || SR1l_I_equal3_cut_Pt_300_350lessMETless450()) && myEvent.topnessMod<0 );}
bool I_equal_34_cut_Pt_300_450_550_modtop_Inf0() {return ( (SR1l_I_equal4_cut_Pt_300_450lessMETless550() || SR1l_I_equal3_cut_Pt_300_450lessMETless550()) && myEvent.topnessMod<0 );}
bool I_equal_34_cut_Pt_300_550_Inf_modtop_Inf0() {return ( (SR1l_I_equal4_cut_Pt_300_550lessMETlessInf() || SR1l_I_equal3_cut_Pt_300_550lessMETlessInf()) && myEvent.topnessMod<0 );}





bool I_equal_34_cut_Pt_300_modtop_Sup0() {return ( (SR1l_I_equal_4_cut_Pt_300() || SR1l_I_equal_3_cut_Pt_300()) && myEvent.topnessMod>=0);}


bool I_equal_34_cut_Pt_300_250_350_modtop_Sup0() {return ( (SR1l_I_equal4_cut_Pt_300_250lessMETless350() ||  SR1l_I_equal3_cut_Pt_300_250lessMETless350()) && myEvent.topnessMod>=0 );}
bool I_equal_34_cut_Pt_300_350_450_modtop_Sup0() {return ( (SR1l_I_equal4_cut_Pt_300_350lessMETless450() || SR1l_I_equal3_cut_Pt_300_350lessMETless450()) && myEvent.topnessMod>=0 );}
bool I_equal_34_cut_Pt_300_450_550_modtop_Sup0() {return ( (SR1l_I_equal4_cut_Pt_300_450lessMETless550() || SR1l_I_equal3_cut_Pt_300_450lessMETless550()) && myEvent.topnessMod>=0 );}
bool I_equal_34_cut_Pt_300_550_Inf_modtop_Sup0() {return ( (SR1l_I_equal4_cut_Pt_300_550lessMETlessInf() || SR1l_I_equal3_cut_Pt_300_550lessMETlessInf()) && myEvent.topnessMod>=0 );}










//-----------Loose b tagging------------------

bool SR1l_I_equal3_cut_Pt_300_250lessMETless350_looseB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=250 && myEvent.pfmet<350 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal3_cut_Pt_300_350lessMETless450_looseB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal3_cut_Pt_300_450lessMETless550_looseB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2  && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=450 && myEvent.pfmet<550 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal3_cut_Pt_300_550lessMETlessInf_looseB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=550 && myEvent.ak4pfjets_pt >=300);}



bool SR1l_I_equal4_cut_Pt_300_250lessMETless350_looseB() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=250 && myEvent.pfmet<350 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal4_cut_Pt_300_350lessMETless450_looseB() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal4_cut_Pt_300_450lessMETless550_looseB() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=450 && myEvent.pfmet<550 && myEvent.ak4pfjets_pt >=300);}
bool SR1l_I_equal4_cut_Pt_300_550lessMETlessInf_looseB() { return (SR1l() && myEvent.ngoodjets==4  && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_leadJet_identification<=0.5426 && myEvent.pfmet>=550 && myEvent.ak4pfjets_pt >=300);}


bool SR1l_I_equal_4_cut_Pt_300_looseB() {return (SR1l_I_equal4_cut_Pt_300_250lessMETless350_looseB() || SR1l_I_equal4_cut_Pt_300_350lessMETless450_looseB() || SR1l_I_equal4_cut_Pt_300_450lessMETless550_looseB() || SR1l_I_equal4_cut_Pt_300_550lessMETlessInf_looseB() );}

bool SR1l_I_equal_3_cut_Pt_300_looseB() {return (SR1l_I_equal3_cut_Pt_300_250lessMETless350_looseB() || SR1l_I_equal3_cut_Pt_300_350lessMETless450_looseB() || SR1l_I_equal3_cut_Pt_300_450lessMETless550_looseB() || SR1l_I_equal3_cut_Pt_300_550lessMETlessInf_looseB() );}


bool I_equal_34_cut_Pt_300_looseB() {return ( SR1l_I_equal_4_cut_Pt_300_looseB() || SR1l_I_equal_3_cut_Pt_300_looseB() );}


bool I_equal_34_cut_Pt_300_250_350_looseB() {return ( SR1l_I_equal4_cut_Pt_300_250lessMETless350_looseB() ||  SR1l_I_equal3_cut_Pt_300_250lessMETless350_looseB() );}
bool I_equal_34_cut_Pt_300_350_450_looseB() {return ( SR1l_I_equal4_cut_Pt_300_350lessMETless450_looseB() || SR1l_I_equal3_cut_Pt_300_350lessMETless450_looseB() );}
bool I_equal_34_cut_Pt_300_450_550_looseB() {return ( SR1l_I_equal4_cut_Pt_300_450lessMETless550_looseB() || SR1l_I_equal3_cut_Pt_300_450lessMETless550_looseB() );}
bool I_equal_34_cut_Pt_300_550_Inf_looseB() {return ( SR1l_I_equal4_cut_Pt_300_550lessMETlessInf_looseB() || SR1l_I_equal3_cut_Pt_300_550lessMETlessInf_looseB() );}






bool SR1l_I_250lessMETless350_topmod_Sup0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350 && myEvent.topnessMod>=0);}
bool SR1l_I_350lessMETless450_topmod_Sup0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.topnessMod>=0 );}
bool SR1l_I_450lessMETless550_topmod_Sup0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550 && myEvent.topnessMod>=0 );}
bool SR1l_I_550lessMETlessInf_topmod_Sup0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550 && myEvent.topnessMod>=0 );}




bool SR1l_I_250lessMETless350_topmod_Inf0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350 && myEvent.topnessMod<0);}
bool SR1l_I_350lessMETless450_topmod_Inf0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450 && myEvent.topnessMod<0 );}
bool SR1l_I_450lessMETless550_topmod_Inf0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550 && myEvent.topnessMod<0 );}
bool SR1l_I_550lessMETlessInf_topmod_Inf0() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550 && myEvent.topnessMod<0 );}


bool I_topmod_Inf0() { return (SR1l_I_250lessMETless350_topmod_Inf0() && SR1l_I_350lessMETless450_topmod_Inf0() && SR1l_I_450lessMETless550_topmod_Inf0() && SR1l_I_550lessMETlessInf_topmod_Inf0());}
bool I_topmod_Sup0() { return (SR1l_I_250lessMETless350_topmod_Sup0() && SR1l_I_350lessMETless450_topmod_Sup0() && SR1l_I_450lessMETless550_topmod_Sup0() && SR1l_I_550lessMETlessInf_topmod_Sup0());}




bool SR1l_I_250lessMETless350_nocutB() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450_nocutB() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550_nocutB() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInf_nocutB() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=550);}


bool SR1l_I4_250lessMETless350_nocutB() { return (SR1l() && myEvent.ngoodjets==4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I4_350lessMETless450_nocutB() { return (SR1l() && myEvent.ngoodjets==4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I4_450lessMETless550_nocutB() { return (SR1l() && myEvent.ngoodjets==4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I4_550lessMETlessInf_nocutB() { return (SR1l() && myEvent.ngoodjets==4 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=550);}


bool SR1l_I3_250lessMETless350_nocutB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I3_350lessMETless450_nocutB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I3_450lessMETless550_nocutB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I3_550lessMETlessInf_nocutB() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.pfmet>=550);}


bool SR1l_I34_250lessMETless350_nocutB() { return (SR1l_I3_250lessMETless350_nocutB() || SR1l_I4_250lessMETless350_nocutB() );}
bool SR1l_I34_350lessMETless450_nocutB() { return (SR1l_I3_350lessMETless450_nocutB() || SR1l_I4_350lessMETless450_nocutB() );}
bool SR1l_I34_450lessMETless550_nocutB() { return (SR1l_I3_450lessMETless550_nocutB() || SR1l_I4_450lessMETless550_nocutB() );}
bool SR1l_I34_550lessMETlessInf_nocutB() { return (SR1l_I3_550lessMETlessInf_nocutB() || SR1l_I4_550lessMETlessInf_nocutB() );}


bool SR1l_I34_nucutB() { return ( SR1l_I34_250lessMETless350_nocutB() || SR1l_I34_350lessMETless450_nocutB() || SR1l_I34_450lessMETless550_nocutB() || SR1l_I34_550lessMETlessInf_nocutB() ) ;}






bool SR1l_I_250lessMETless350_NocutPt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450_NocutPt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550_NocutPt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInf_NocutPt() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}
