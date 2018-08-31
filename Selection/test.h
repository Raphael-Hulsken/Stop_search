#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

bool CheckTrigger(bool data, string dataset){
       if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");
        
        bool trigged = false;

        if(data)
        {
          if( dataset.find("data_single_electron")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl)
          {
            //cout << "electron trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << endl; 
            trigged = true;
          }
          else if( dataset.find("data_single_muon")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
          {
            //cout << "muon trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << endl; 
            trigged = true;
          }
          else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET100_MHT100 || myEvent.HLT_MET) && ( (fabs(myEvent.lep1_pdgid)==11 && !myEvent.HLT_SingleEl ) ||  (fabs(myEvent.lep1_pdgid)==13 && !myEvent.HLT_SingleMu) ))
          {
            //cout << "HT trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << " MET trigger " << myEvent.HLT_MET << endl; 
            trigged = true;
          }
          else
          {
            //cout << "no trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << " MET trigger " << myEvent.HLT_MET << endl; 
     
          }

         }
         else
         {
             trigged = true;
         }
         return trigged;

}

bool baseline(){ return myEvent.pfmet>=250 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.MT2W>=0;}
bool SR1l() { return ( baseline() && myEvent.ngoodbtags>=1 && OneLep() && ( (myEvent.lep1_passMediumID && abs(myEvent.lep1_pdgid)==11) || (myEvent.lep1_passTightID && abs(myEvent.lep1_pdgid)==13))  ); }
bool CR1l() { return ( baseline() && myEvent.ngoodbtags==0 && OneLep() ); }
bool CR2l() { return ( baseline() && myEvent.ngoodbtags>=1 && TwoLep() ); }
bool SR1l2j_MET250to350() { return (SR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l2j_MET350to450() { return (SR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l2j_MET450toInf() { return (SR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=450);}
bool SR1l3j_MET250to350() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l3j_MET350to450() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l3j_MET450to550() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l3j_MET550toInf() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=550);}
bool SR1l4jLow_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jLow_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jLow_MET450toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=450);}
bool SR1l4jHighMT2W_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jHighMT2W_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jHighMT2W_MET450to550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l4jHighMT2W_MET550to650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l4jHighMT2W_MET650toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=650);}
bool CR1l2j_MET250to350() { return (CR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l2j_MET350to450() { return (CR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l2j_MET450toInf() { return (CR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=450);}
bool CR1l3j_MET250to350() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l3j_MET350to450() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l3j_MET450to550() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l3j_MET550toInf() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=550);}
bool CR1l4jLow_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jLow_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jLow_MET450toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=450);}
bool CR1l4jHighMT2W_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jHighMT2W_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jHighMT2W_MET450to550() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l4jHighMT2W_MET550to650() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR1l4jHighMT2W_MET650toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=650);}
bool CR2l2j_MET250to350() { return (CR2l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l2j_MET350to450() { return (CR2l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l2j_MET450toInf() { return (CR2l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=450);}
bool CR2l3j_MET250to350() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l3j_MET350to450() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l3j_MET450to550() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l3j_MET550toInf() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=550);}
bool CR2l4jLow_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jLow_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jLow_MET450toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=450);}
bool CR2l4jHighMT2W_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jHighMT2W_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jHighMT2W_MET450to550() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l4jHighMT2W_MET550to650() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR2l4jHighMT2W_MET650toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=650);}
