#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
//bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

bool CheckTrigger(bool data, string dataset){
       if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");
        
        bool trigged = false;

        if(data)
        {
            if( dataset.find("data_single_electron")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl)
                trigged = true;
            else if( dataset.find("data_single_muon")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
                trigged = true;
            else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET100_MHT100 || myEvent.HLT_MET) && ( (fabs(myEvent.lep1_pdgid)==11 && !myEvent.HLT_SingleEl ) ||  (fabs(myEvent.lep1_pdgid)==13 && !myEvent.HLT_SingleMu) ))
                trigged = true;
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

bool tightCSVV2()
{
    for (uint32_t j =0 ; j<myEvent.ak4pfjets_CSV->size(); j++)
    {
        if( myEvent.ak4pfjets_CSV->at(j) >0.935)
            return true;
    }
    return false;
}

bool baseline(){ return myEvent.pfmet>=250 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0;}
bool SR1l() { return ( baseline() && myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && ( (myEvent.lep1_passMediumID && abs(myEvent.lep1_pdgid)==11) || (myEvent.lep1_passTightID && abs(myEvent.lep1_pdgid)==13) ) ); }
bool SR1l_AB_250lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
bool SR1l_CD_250lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<0 && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
bool SR1l_EFGH_250lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0 && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
bool SR1l_I_250lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.dphiMET<2 && myEvent.pfmet>=250);}
bool SR1l_CDEFGH_250lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
bool SR1l_NJlowTM_250lessMETlessInf() { return (SR1l() && myEvent.topnessMod<0 && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
bool SR1l_NJmidTM_250lessMETlessInf() { return (SR1l() && myEvent.topnessMod>0 && myEvent.topnessMod<10 && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
bool SR1l_NJhighTM_250lessMETlessInf() { return (SR1l() && myEvent.topnessMod>10 && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
bool SR1l_NJ_250lessMETlessInf() { return (SR1l() && myEvent.dphi_ak4pfjets_met>=0.8 && ( (myEvent.Mlb>175 && tightCSVV2()) || myEvent.Mlb<=175 ) && myEvent.pfmet>=250);}
