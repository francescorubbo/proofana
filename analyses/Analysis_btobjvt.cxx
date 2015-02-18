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

  for(int iJet = 0; iJet < jets("AntiKt4LCTopo"); iJet++){
    Particle *myjet = &(jet(iJet, "AntiKt4LCTopo"));
    vector<int> assoc_trk_indices;
    for(int it=0; it<myjet->Int("nTrackAssoc"); ++it){
      Particle* trk = (Particle*) myjet->Obj("GhostAssocTrack", it);
      assoc_trk_indices.push_back(trk->Int("JVTindex"));
    }
    
    (*jvt)(myjet->p.Pt(),assoc_trk_indices);
    myjet->Set("JVT",jvt->JVT());

  }
  for(int iJet = 0; iJet < jets("AntiKt4LCTopo"); iJet++){
    Particle *myjet = &(jet(iJet, "AntiKt4LCTopo"));
    float maxdPhi = 0;
    Particle* backtobackjet =0;
    for(int iJ=0; iJ<jets("AntiKt4LCTopo"); ++iJ){
      Particle *otherjet = &(jet(iJ, "AntiKt4LCTopo"));
      if(myjet == otherjet)                            continue;
      if(fabs(otherjet -> p.Eta())>2.4)                continue; // only looking at jets within the tracker
      
      if(otherjet -> p.Pt()<20)                        continue; // only looking at jets with pT > 30 GeV
      
      if(fabs(otherjet->p.DeltaPhi(myjet->p))<maxdPhi) continue; // keep track of jet that is most back to back
      
      maxdPhi = fabs(otherjet->p.DeltaPhi(myjet->p));
      backtobackjet = otherjet;
    }
    if(backtobackjet!=0){
      myjet->Set("BtoBJVT", backtobackjet->Float("JVT"));
      myjet->Set("BtoBjet", backtobackjet);
    }
    else
      myjet->Set("BtoBJVT", -1.);
  }
  
  return true;
}

/// WorkerTerminate: clean up
///=========================================
void Analysis_btobjvt::WorkerTerminate()
{


  // Nothing more

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
