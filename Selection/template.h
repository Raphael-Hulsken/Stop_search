#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;



//bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
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
    if(!myEvent.is_data)
        return true;
    else if(myEvent.filt_met && 
            myEvent.filt_badChargedCandidateFilter && 
            myEvent.filt_jetWithBadMuon && 
            myEvent.filt_pfovercalomet && 
            myEvent.filt_badMuonFilter && 
            !myEvent.filt_duplicatemuons && 
            !myEvent.filt_badmuons && 
            myEvent.filt_nobadmuons )
        return true;
    else
        return false;   
        
}

bool passGoodVtx()
{
    bool selRes = false;
    if(myEvent.nvertex>=1)
        selRes= true;
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

