/**************************************************************************
 **
 **   File:         Analysis_btobjvt.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_btobjvt_cxx

#include "Analysis_btobjvt.h"
#include "AnaConfig.h"
#include <TROOT.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <cstdlib>
#include <sstream>
#include "TKey.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Subtractor.hh"
#include "ANN/ANN.h"

const double PI  =3.141592653589793238463;

///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_btobjvt::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_btobjvt: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  std::string maindir(gSystem->Getenv("PROOFANADIR"));
  if(maindir==".") maindir+="/libProofAna";
  cout << maindir << " is the maindir we work from! " << endl;
  cout << gSystem->Exec(TString("ls "+maindir).Data()) << endl;

  jvt = new JetVertexTagger(0.2, maindir+"/utils/JetVertexTagger/data/JVTlikelihood_20140805.root");

  if (Debug()) cout << "Analysis_btobjvt: DEBUG Finish WorkerBegin()" << endl;  
} 

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_btobjvt::ProcessEvent()
{

  if (Debug()) cout << "Analysis_btobjvt: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();

  // -----------------------------
  // make objects truthsStable: stable truth particles for jets 
  Analysis_JetMET_Base::AddStableParticles();
  Analysis_JetMET_Base::MakeJets(fastjet::antikt_algorithm, 0.2, "truthsStable");
  
  selectTracks();
  
  MomKey newjets;
  newjets = MakeJets(fastjet::antikt_algorithm, 0.2, "tracksforjvf","tracks",2);
  Analysis_pileup::addTruthMatch(newjets,"AntiKt2Truth");
  newjets = MakeJets(fastjet::antikt_algorithm, 0.4, "tracksforjvf","tracks",2);
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");

  vector<float> trk_pt, trk_z0_wrtPV;
  for(int it=0; it< tracks(); ++it){
    trk_pt .push_back(track(it).p.Pt());
    trk_z0_wrtPV .push_back(track(it).Float("z0_wrtPV"));
    track(it).Set("JVTindex", it);
  }
  // trk to vtx assoc
  vector<vector<int> > vxp_trk_index;
  for(int iv=0; iv<vtxs(); ++iv){
    vector<int> assoc_track_indices;
    for(int it=0; it<vtx(iv).Objs("vtxTracks"); ++it){
      Particle* trk = (Particle*) vtx(iv).Obj("vtxTracks",it);
      assoc_track_indices.push_back(trk->Int("JVTindex"));
    }
    vxp_trk_index.push_back(assoc_track_indices);
  }
  // JVT
  jvt->init(trk_pt, trk_z0_wrtPV, vxp_trk_index);

  // jet JVT
  for(int iJet = 0; iJet < jets("AntiKt4LCTopo"); iJet++){
    Particle *myjet = &(jet(iJet, "AntiKt4LCTopo"));
    vector<int> assoc_trk_indices;
    for(int it=0; it<myjet->Int("nTrackAssoc"); ++it){
      Particle* trk = (Particle*) myjet->Obj("GhostAssocTrack", it);
      assoc_trk_indices.push_back(trk->Int("JVTindex"));
    }
    
    (*jvt)(myjet->p.Pt(),assoc_trk_indices);
    myjet->Set("JVT",jvt->JVT());
    myjet->Set("corrJVF",jvt->corrJVF());
    myjet->Set("RpT",jvt->RpT());
  }

  //track-jet JVT
  for(int iJet = 0; iJet < jets("AntiKt2tracks"); iJet++){
    Particle *myjet = &(jet(iJet, "AntiKt2tracks"));
    vector<int> assoc_trk_indices;
    for(int it=0; it<myjet->Objs("constituents"); ++it){
      Particle* trk = (Particle*) myjet->Obj("constituents",it);
      assoc_trk_indices.push_back(trk->Int("JVTindex"));
    }
    
    (*jvt)(myjet->p.Pt(),assoc_trk_indices);
    myjet->Set("JVT",jvt->JVT());
    myjet->Set("corrJVF",jvt->corrJVF());
    myjet->Set("RpT",jvt->RpT());
  }

  for(int iJet = 0; iJet < jets("AntiKt4tracks"); iJet++){
    Particle *myjet = &(jet(iJet, "AntiKt4tracks"));
    vector<int> assoc_trk_indices;
    for(int it=0; it<myjet->Objs("constituents"); ++it){
      Particle* trk = (Particle*) myjet->Obj("constituents",it);
      assoc_trk_indices.push_back(trk->Int("JVTindex"));
    }
    
    (*jvt)(myjet->p.Pt(),assoc_trk_indices);
    myjet->Set("JVT",jvt->JVT());
    myjet->Set("corrJVF",jvt->corrJVF());
    myjet->Set("RpT",jvt->RpT());
  }

  for(int iJet = 0; iJet < jets("AntiKt4LCTopo"); iJet++){
    Particle *myjet = &(jet(iJet, "AntiKt4LCTopo"));

    myjet->Set("BtoBJVT", -2.);
    myjet->Set("BtoBJVT_70", -2.);
    myjet->Set("BtoBJVT_75", -2.);
    myjet->Set("BtoBJVT_85", -2.);
    myjet->Set("BtoBJVT_90", -2.);
    myjet->Set("BtoBcorrJVF", -2.);
    myjet->Set("BtoBRpT", -2.);
    myjet->Set("BtoBJVT_trk4", -2.);
    myjet->Set("BtoBcorrJVF_trk4", -2.);
    myjet->Set("BtoBRpT_trk4", -2.);
    myjet->Set("BtoBJVT_trk2", -2.);
    myjet->Set("BtoBcorrJVF_trk2", -2.);
    myjet->Set("BtoBRpT_trk2", -2.);
    myjet->AddVec("jetsBtoB");
    myjet->AddVec("jetsBtoB_trk2");
    myjet->AddVec("jetsBtoB_trk4");

    if(myjet->p.Eta()<2.4) continue; 

    float mindPhi = 99.;
    float mindPhi70 = 99.;
    float mindPhi75 = 99.;
    float mindPhi85 = 99.;
    float mindPhi90 = 99.;
    Particle* backtobackjet =0;
    for(int iJ=0; iJ<jets("AntiKt4LCTopo"); ++iJ){
      Particle *otherjet = &(jet(iJ, "AntiKt4LCTopo"));
      if(myjet == otherjet)                            continue;
      if(fabs(otherjet -> p.Eta())>2.4)                continue; // only looking at jets within the tracker
      if(otherjet -> p.Pt()<10) continue; // only looking at jets with pT > 30 GeV
      float ptasymm = fabs(otherjet->p.Pt()-myjet->p.Pt())/(otherjet->p.Pt()+myjet->p.Pt());
      // if(ptasymm>0.3) continue;
      // float dPhi = fabs(fabs(otherjet->p.DeltaPhi(myjet->p))-PI);
      // if (dPhi>0.8) continue;
      // if(dPhi>mindPhi) continue; // keep track of jet that is most back to back

      float dPhi = fabs(fabs(otherjet->p.DeltaPhi(myjet->p))-PI);
      if(ptasymm<0.3){
	if(dPhi<1.5) //70% working point
	  if(dPhi<mindPhi70){
	    myjet->Set("BtoBJVT_70", otherjet->Float("JVT"));
	    mindPhi70 = dPhi;
	  }
	if(dPhi<0.8) //75% working point
	  if(dPhi<mindPhi75){
	    myjet->Set("BtoBJVT_75", otherjet->Float("JVT"));
	    mindPhi75 = dPhi;
	  }
      }//ptasymm<0.3
      if(ptasymm<0.2){
	if(dPhi<0.8) //85% working point
	  if(dPhi<mindPhi85){
	    myjet->Set("BtoBJVT_85", otherjet->Float("JVT"));
	    mindPhi85 = dPhi;
	  }
	if(dPhi<0.6) //90% working point
	  if(dPhi<mindPhi90){
	    myjet->Set("BtoBJVT_90", otherjet->Float("JVT"));
	    mindPhi90 = dPhi;
	  }
      }//ptasymm<0.3      
    }
    if(backtobackjet!=0){
      myjet->Set("BtoBJVT", backtobackjet->Float("JVT"));
      myjet->Set("BtoBcorrJVF", backtobackjet->Float("corrJVF"));
      myjet->Set("BtoBRpT", backtobackjet->Float("RpT"));
      myjet->Add("jetsBtoB",backtobackjet);
    }

    //compute BtoBJVT w/o jet
    float maxdPhi_trk2 = 0;
    Particle* backtobackjet_trk2 =0;
    for(int iJ=0; iJ<jets("AntiKt2tracks"); ++iJ){
      Particle *otherjet = &(jet(iJ, "AntiKt2tracks"));
      if(fabs(otherjet -> p.Eta())>2.4)                continue; // only looking at jets within the tracker
      // if(otherjet -> p.Pt()<20)                        continue; // only looking at jets with pT > 30 GeV
      if(fabs(otherjet->p.DeltaPhi(myjet->p))<maxdPhi_trk2) continue; // keep track of jet that is most back to back
      
      maxdPhi_trk2 = fabs(otherjet->p.DeltaPhi(myjet->p));
      backtobackjet_trk2 = otherjet;
    }
    if(backtobackjet_trk2!=0){
      myjet->Set("BtoBJVT_trk2", backtobackjet_trk2->Float("JVT"));
      myjet->Set("BtoBcorrJVF_trk2", backtobackjet_trk2->Float("corrJVF"));
      myjet->Set("BtoBRpT_trk2", backtobackjet_trk2->Float("RpT"));
      myjet->Add("jetsBtoB_trk2",backtobackjet_trk2);
    }
    float maxdPhi_trk4 = 0;
    Particle* backtobackjet_trk4 =0;
    for(int iJ=0; iJ<jets("AntiKt4tracks"); ++iJ){
      Particle *otherjet = &(jet(iJ, "AntiKt4tracks"));
      if(fabs(otherjet -> p.Eta())>2.4)                continue; // only looking at jets within the tracker
      // if(otherjet -> p.Pt()<20)                        continue; // only looking at jets with pT > 30 GeV
      if(fabs(otherjet->p.DeltaPhi(myjet->p))<maxdPhi_trk4) continue; // keep track of jet that is most back to back
      
      maxdPhi_trk4 = fabs(otherjet->p.DeltaPhi(myjet->p));
      backtobackjet_trk4 = otherjet;
    }
    if(backtobackjet_trk4!=0){
      myjet->Set("BtoBJVT_trk4", backtobackjet_trk4->Float("JVT"));
      myjet->Set("BtoBcorrJVF_trk4", backtobackjet_trk4->Float("corrJVF"));
      myjet->Set("BtoBRpT_trk4", backtobackjet_trk4->Float("RpT"));
      myjet->Set("BtoBjet_trk4", backtobackjet_trk4);
      myjet->Add("jetsBtoB_trk4",backtobackjet_trk4);
    }
    
  }
  return true;
}

/// WorkerTerminate: clean up
///=========================================
void Analysis_btobjvt::WorkerTerminate()
{


  // Nothing more

}

MomKey Analysis_btobjvt::MakeJets(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType, const MomKey extra, double minpt){

  const static MomKey SJetKey("jets");
  vector<fastjet::PseudoJet> inputConst = ObjsToPJ(constType);

  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence clustSeq(inputConst, jetDef);  

  vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(minpt));

  TString key;

  switch(algo){
    case fastjet::antikt_algorithm:
      key = "AntiKt";
      break;
    
    case fastjet::kt_algorithm:
      key = "Kt";
      break;
    
    case fastjet::cambridge_algorithm:
      key = "CamKt";
      break;
    
    default:
      cout << " this jet algorithm is not supported! quitting " << endl;
      exit(-1);
      break;
  }
  
  key+= TString::Format("%.0f",10.*jetR);

  const static MomKey LCKey("clustersLCTopo");
  const static MomKey TrackKey("tracksgood");
  const static MomKey TruthKey("truthsStable");

  if(constType==LCKey){
  	key+="LCTopo";
  } else if (constType==TrackKey){
  	key+="TrackZ";
  } else if(constType==TruthKey){
    key+="Truth";
  }
  key+=extra;

  if(Debug()) cout << "MakeJets with key " << key << endl;

  MomKey FinalKey(key);
  MomKey FFinalKey = SJetKey + FinalKey;

  AddVec(SJetKey+FinalKey);

  for(unsigned int iJet = 0 ; iJet < inclusiveJets.size() ; iJet++){
  	fastjet::PseudoJet jet = inclusiveJets[iJet];
  	vector<fastjet::PseudoJet> constituents = jet.constituents();
  	Particle* jetP = new Particle();
  	jetP->p.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.e());
  	static const MomKey ConsKey("constituents");
  	jetP->AddVec(ConsKey);
  	for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
  		const PJ_Info* info = &(constituents[iCons].user_info<PJ_Info>());
  		jetP->Add(ConsKey, info->Pointer);
  	} // end loop over cons
  	Add(FFinalKey, jetP);
  }// end loop over jets

  return FinalKey;
}

 void Analysis_btobjvt::selectTracks()
{
  const MomKey tracksforjvf("tracksforjvf");
  AddVec(tracksforjvf);
  for(int iTr = 0; iTr < tracks(); iTr++){
    if(fabs(track(iTr).Float("d0_wrtPV")) > 2.5) continue;
    if(track(iTr).p.Pt() < 0.5) continue;
    if(fabs(track(iTr).p.Eta()) > 2.5) continue;
    if(track(iTr).Int("nPixHits")+track(iTr).Int("nPixelDeadSensors") < 1) continue;
    if(track(iTr).Int("nSCTHits")+track(iTr).Int("nSCTDeadSensors") < 6) continue;
    if(track(iTr).Float("chi2")/((float)track(iTr).Int("ndof")) > 5.) continue;
    Add(tracksforjvf,&track(iTr));
  }
}
