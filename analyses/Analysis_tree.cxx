/**************************************************************************
 **
 **   File:         Analysis_tree.cxx
 **
 **   Description:  See header
 **
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_tree_cxx

#include "Analysis_tree.h"
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

///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_tree::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_tree: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  OutputDir()->cd();
  fTree = new TTree("tree", "tree");
  AddBranches(fTree);

  cljvfcorr   = new vector<float>();  
  clfem   = new vector<float>();  
  clcenterlambda   = new vector<float>();  
  clpt     = new vector<float>();  
  cleta    = new vector<float>();  
  clphi    = new vector<float>();  
  clenergy = new vector<float>();  

  trkpt     = new vector<float>();  
  trketa    = new vector<float>();  
  trkphi    = new vector<float>();  
  trkispv   = new vector<bool>();  

  j0ncl     = new vector<int>();  
  j5ncl     = new vector<int>();  
  j0mass    = new vector<float>();  
  j5mass    = new vector<float>();  
  j0width   = new vector<float>();  
  j5width   = new vector<float>();  
  j0pt      = new vector<float>();  
  j0eta     = new vector<float>();  
  j0phi     = new vector<float>();  
  j0jvt     = new vector<float>();  
  j5pt      = new vector<float>();  
  j5eta     = new vector<float>();  
  j5phi     = new vector<float>();  
  tj0pt     = new vector<float>();  
  tj0eta    = new vector<float>();  
  tj5pt     = new vector<float>();  
  tj5eta    = new vector<float>();  
  
  lj0ncl     = new vector<int>();  
  lj5ncl     = new vector<int>();  
  lj0mass    = new vector<float>();  
  lj5mass    = new vector<float>();  
  lj0width   = new vector<float>();  
  lj5width   = new vector<float>();  
  lj0pt      = new vector<float>();  
  lj0eta     = new vector<float>();  
  lj0phi     = new vector<float>();  
  lj5pt      = new vector<float>();  
  lj5eta     = new vector<float>();  
  lj5phi     = new vector<float>();  
  tlj0pt     = new vector<float>();  
  tlj0eta    = new vector<float>();  
  tlj5pt     = new vector<float>();  
  tlj5eta    = new vector<float>();  

  jnoarea0ncl     = new vector<int>();  
  jnoarea5ncl     = new vector<int>();  
  jnoarea0mass    = new vector<float>();  
  jnoarea5mass    = new vector<float>();  
  jnoarea0width   = new vector<float>();  
  jnoarea5width   = new vector<float>();  
  jnoarea0pt      = new vector<float>();  
  jnoarea0eta     = new vector<float>();  
  jnoarea0phi     = new vector<float>();  
  jnoarea5pt      = new vector<float>();  
  jnoarea5eta     = new vector<float>();  
  jnoarea5phi     = new vector<float>();  
  tjnoarea0pt     = new vector<float>();  
  tjnoarea0eta    = new vector<float>();  
  tjnoarea5pt     = new vector<float>();  
  tjnoarea5eta    = new vector<float>();  

  ljnoarea0ncl     = new vector<int>();  
  ljnoarea5ncl     = new vector<int>();  
  ljnoarea0mass    = new vector<float>();  
  ljnoarea5mass    = new vector<float>();  
  ljnoarea0width   = new vector<float>();  
  ljnoarea5width   = new vector<float>();  
  ljnoarea0pt      = new vector<float>();  
  ljnoarea0eta     = new vector<float>();  
  ljnoarea0phi     = new vector<float>();  
  ljnoarea5pt      = new vector<float>();  
  ljnoarea5eta     = new vector<float>();  
  ljnoarea5phi     = new vector<float>();  
  tljnoarea0pt     = new vector<float>();  
  tljnoarea0eta    = new vector<float>();  
  tljnoarea5pt     = new vector<float>();  
  tljnoarea5eta    = new vector<float>();  

  truejetpt      = new vector<float>();  
  truejeteta     = new vector<float>();  
  truejetphi     = new vector<float>();  
  truejetenergy  = new vector<float>();  

  truelargejetpt      = new vector<float>();  
  truelargejeteta     = new vector<float>();  
  truelargejetphi     = new vector<float>();  
  truelargejetenergy  = new vector<float>();  

  if (Debug()) cout << "Analysis_tree: DEBUG Finish WorkerBegin()" << endl;
}

///=========================================
/// ProcessEvent: run the analysis
///=========================================
bool Analysis_tree::ProcessEvent()
{

  if (Debug()) cout << "Analysis_pileup: DEBUG In ProcessEvent(): RunNumber = " << RunNumber()
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();

  // -----------------------------
  // make objects truthsStable: stable truth particles for jets
  Analysis_JetMET_Base::AddStableParticles();

  FillTree("LCTopo");

  return true;
}

/// WorkerTerminate: clean up
///=========================================
void Analysis_tree::WorkerTerminate()
{

  fTree->Write();
  // Nothing more
  delete cljvfcorr;  
  delete clfem    ;  
  delete clcenterlambda    ;  
  delete clpt    ;  
  delete cleta   ;  
  delete clphi   ;  
  delete clenergy;  
  delete trkpt    ;  
  delete trketa   ;  
  delete trkphi   ;  
  delete trkispv  ;  

  delete j0ncl    ;  
  delete j5ncl    ;  
  delete j0mass    ;  
  delete j5mass    ;  
  delete j0width    ;  
  delete j5width    ;  
  delete j0pt    ;  
  delete j0eta    ;  
  delete j0phi    ;  
  delete j0jvt    ;  
  delete j5pt    ;  
  delete j5eta    ;  
  delete j5phi    ;  
  delete tj0pt    ;  
  delete tj5pt    ;  
  delete tj0eta    ;  
  delete tj5eta    ;  

  delete lj0ncl    ;  
  delete lj5ncl    ;  
  delete lj0mass    ;  
  delete lj5mass    ;  
  delete lj0width    ;  
  delete lj5width    ;  
  delete lj0pt    ;  
  delete lj0eta    ;  
  delete lj0phi    ;  
  delete lj5pt    ;  
  delete lj5eta    ;  
  delete lj5phi    ;  
  delete tlj0pt    ;  
  delete tlj5pt    ;  
  delete tlj0eta    ;  
  delete tlj5eta    ;  

  delete jnoarea0ncl    ;  
  delete jnoarea5ncl    ;  
  delete jnoarea0mass    ;  
  delete jnoarea5mass    ;  
  delete jnoarea0width    ;  
  delete jnoarea5width    ;  
  delete jnoarea0pt    ;  
  delete jnoarea0eta    ;  
  delete jnoarea0phi    ;  
  delete jnoarea5pt    ;  
  delete jnoarea5eta    ;  
  delete jnoarea5phi    ;  
  delete tjnoarea0pt    ;  
  delete tjnoarea5pt    ;  
  delete tjnoarea0eta    ;  
  delete tjnoarea5eta    ;  

  delete ljnoarea0ncl    ;  
  delete ljnoarea5ncl    ;  
  delete ljnoarea0mass    ;  
  delete ljnoarea5mass    ;  
  delete ljnoarea0width    ;  
  delete ljnoarea5width    ;  
  delete ljnoarea0pt    ;  
  delete ljnoarea0eta    ;  
  delete ljnoarea0phi    ;  
  delete ljnoarea5pt    ;  
  delete ljnoarea5eta    ;  
  delete ljnoarea5phi    ;  
  delete tljnoarea0pt    ;  
  delete tljnoarea5pt    ;  
  delete tljnoarea0eta    ;  
  delete tljnoarea5eta    ;  

  delete truejetpt    ;  
  delete truejeteta    ;  
  delete truejetphi    ;  
  delete truejetenergy    ;  

  delete truelargejetpt    ;  
  delete truelargejeteta    ;  
  delete truelargejetphi    ;  
  delete truelargejetenergy    ;  

}

void Analysis_tree::FillTree(const MomKey Key){
  if(Debug()) cout <<"Analysis_tree::FillTree Start" << endl;
  ResetBranches(fTree);
  FillEventVars(fTree);
  for(int i=0; i<clusters(Key); ++i)
    FillClVars(fTree, i, Key);
  for(int i=0; i<tracks(); ++i)
    FillTrkVars(fTree, i);

  for(int i=0; i<jets("AntiKt4jvf0"); ++i){
    FillJetVars(fTree, i, "AntiKt4jvf0",j0pt,j0eta,j0phi,tj0pt,tj0eta,j0ncl,j0mass,j0width);
    Particle  *myjet         = &(jet(i,"AntiKt4jvf0"));
    j0jvt->push_back(myjet->Float("JVT"));
  }
  
  for(int i=0; i<jets("AntiKt4jvf5"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf5",j5pt,j5eta,j5phi,tj5pt,tj5eta,j5ncl,j5mass,j5width);

  for(int i=0; i<jets("AntiKt4jvf0noarea"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf0noarea",jnoarea0pt,jnoarea0eta,jnoarea0phi,tjnoarea0pt,tjnoarea0eta,jnoarea0ncl,jnoarea0mass,jnoarea0width);
  for(int i=0; i<jets("AntiKt4jvf5noarea"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf5noarea",jnoarea5pt,jnoarea5eta,jnoarea5phi,tjnoarea5pt,tjnoarea5eta,jnoarea5ncl,jnoarea5mass,jnoarea5width);

  for(int i=0; i<jets("AntiKt10jvf0"); ++i)
    FillJetVars(fTree, i, "AntiKt10jvf0",lj0pt,lj0eta,lj0phi,tlj0pt,tlj0eta,lj0ncl,lj0mass,lj0width,"AntiKt10Truth_match");
  for(int i=0; i<jets("AntiKt10jvf5"); ++i)
    FillJetVars(fTree, i, "AntiKt10jvf5",lj5pt,lj5eta,lj5phi,tlj5pt,tlj5eta,lj5ncl,lj5mass,lj5width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt10jvf0noarea"); ++i)
    FillJetVars(fTree, i, "AntiKt10jvf0noarea",ljnoarea0pt,ljnoarea0eta,ljnoarea0phi,tljnoarea0pt,tljnoarea0eta,ljnoarea0ncl,ljnoarea0mass,ljnoarea0width,"AntiKt10Truth_match");
  for(int i=0; i<jets("AntiKt10jvf5noarea"); ++i)
    FillJetVars(fTree, i, "AntiKt10jvf5noarea",ljnoarea5pt,ljnoarea5eta,ljnoarea5phi,tljnoarea5pt,tljnoarea5eta,ljnoarea5ncl,ljnoarea5mass,ljnoarea5width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt4Truth"); ++i)
    FillTrueJetVars(fTree, i, "AntiKt4Truth");

  for(int i=0; i<jets("AntiKt10Truth"); ++i)
    FillTrueLargeJetVars(fTree, i, "AntiKt10Truth");
  
  fTree->Fill();

  if(Debug()) cout <<"Analysis_tree::FillTree End" << endl;
  return;
}

 void Analysis_tree::AddBranches(TTree *tree){
   if(Debug()) cout <<"Analysis_tree::AddBranches Start" << endl;

   // Event Info
   tree->Branch("EventNumber", &fTEventNumber,            "EventNumber/I");
   tree->Branch("RunNumber",   &fTRunNumber,              "RunNumber/I");
   tree->Branch("Weight" ,     &fTWeight,                 "Weight/F");
   tree->Branch("NPVtruth" ,   &fTNPVtruth,               "NPVtruth/F");
   tree->Branch("NPV" ,   &fTNPV,               "NPV/F");
   tree->Branch("Mu" ,         &fTMu,                     "Mu/F");
   tree->Branch("Rho" ,         &fTRho,                     "Rho/F");

   // Jet vars ------------------------------------------------------------

   tree->Branch("cljvfcorr","std::vector<float>",&cljvfcorr);
   tree->Branch("clfem","std::vector<float>",&clfem);
   tree->Branch("clcenterlambda","std::vector<float>",&clcenterlambda);
   tree->Branch("clpt","std::vector<float>",&clpt);
   tree->Branch("cleta","std::vector<float>",&cleta);
   tree->Branch("clphi","std::vector<float>",&clphi);
   tree->Branch("clenergy","std::vector<float>",&clenergy);
   tree->Branch("trkpt" ,"std::vector<float>",&trkpt);
   tree->Branch("trketa","std::vector<float>",&trketa);
   tree->Branch("trkphi","std::vector<float>",&trkphi);
   tree->Branch("trkispv","std::vector<bool>",&trkispv);

   tree->Branch("j0ncl","std::vector<int>",  &j0ncl);
   tree->Branch("j5ncl","std::vector<int>",  &j5ncl);
   tree->Branch("j0mass","std::vector<float>",  &j0mass);
   tree->Branch("j5mass","std::vector<float>",  &j5mass);
   tree->Branch("j0width","std::vector<float>",  &j0width);
   tree->Branch("j5width","std::vector<float>",  &j5width);
   tree->Branch("j0pt","std::vector<float>", &j0pt);
   tree->Branch("j0eta","std::vector<float>",&j0eta);
   tree->Branch("j0phi","std::vector<float>",&j0phi);
   tree->Branch("j0jvt","std::vector<float>",&j0jvt);
   tree->Branch("j5pt","std::vector<float>", &j5pt);
   tree->Branch("j5eta","std::vector<float>",&j5eta);
   tree->Branch("j5phi","std::vector<float>",&j5phi);
   tree->Branch("tj0pt","std::vector<float>",&tj0pt);
   tree->Branch("tj5pt","std::vector<float>",&tj5pt);
   tree->Branch("tj0eta","std::vector<float>",&tj0eta);
   tree->Branch("tj5eta","std::vector<float>",&tj5eta);

   tree->Branch("lj0ncl","std::vector<int>",  &lj0ncl);
   tree->Branch("lj5ncl","std::vector<int>",  &lj5ncl);
   tree->Branch("lj0mass","std::vector<float>",  &lj0mass);
   tree->Branch("lj5mass","std::vector<float>",  &lj5mass);
   tree->Branch("lj0width","std::vector<float>",  &lj0width);
   tree->Branch("lj5width","std::vector<float>",  &lj5width);
   tree->Branch("lj0pt","std::vector<float>", &lj0pt);
   tree->Branch("lj0eta","std::vector<float>",&lj0eta);
   tree->Branch("lj0phi","std::vector<float>",&lj0phi);
   tree->Branch("lj5pt","std::vector<float>", &lj5pt);
   tree->Branch("lj5eta","std::vector<float>",&lj5eta);
   tree->Branch("lj5phi","std::vector<float>",&lj5phi);
   tree->Branch("tlj0pt","std::vector<float>",&tlj0pt);
   tree->Branch("tlj5pt","std::vector<float>",&tlj5pt);
   tree->Branch("tlj0eta","std::vector<float>",&tlj0eta);
   tree->Branch("tlj5eta","std::vector<float>",&tlj5eta);

   tree->Branch("jnoarea0ncl","std::vector<int>",  &jnoarea0ncl);
   tree->Branch("jnoarea5ncl","std::vector<int>",  &jnoarea5ncl);
   tree->Branch("jnoarea0mass","std::vector<float>",  &jnoarea0mass);
   tree->Branch("jnoarea5mass","std::vector<float>",  &jnoarea5mass);
   tree->Branch("jnoarea0width","std::vector<float>",  &jnoarea0width);
   tree->Branch("jnoarea5width","std::vector<float>",  &jnoarea5width);
   tree->Branch("jnoarea0pt","std::vector<float>", &jnoarea0pt);
   tree->Branch("jnoarea0eta","std::vector<float>",&jnoarea0eta);
   tree->Branch("jnoarea0phi","std::vector<float>",&jnoarea0phi);
   tree->Branch("jnoarea5pt","std::vector<float>", &jnoarea5pt);
   tree->Branch("jnoarea5eta","std::vector<float>",&jnoarea5eta);
   tree->Branch("jnoarea5phi","std::vector<float>",&jnoarea5phi);
   tree->Branch("tjnoarea0pt","std::vector<float>",&tjnoarea0pt);
   tree->Branch("tjnoarea5pt","std::vector<float>",&tjnoarea5pt);
   tree->Branch("tjnoarea0eta","std::vector<float>",&tjnoarea0eta);
   tree->Branch("tjnoarea5eta","std::vector<float>",&tjnoarea5eta);

   tree->Branch("ljnoarea0ncl","std::vector<int>",  &ljnoarea0ncl);
   tree->Branch("ljnoarea5ncl","std::vector<int>",  &ljnoarea5ncl);
   tree->Branch("ljnoarea0mass","std::vector<float>",  &ljnoarea0mass);
   tree->Branch("ljnoarea5mass","std::vector<float>",  &ljnoarea5mass);
   tree->Branch("ljnoarea0width","std::vector<float>",  &ljnoarea0width);
   tree->Branch("ljnoarea5width","std::vector<float>",  &ljnoarea5width);
   tree->Branch("ljnoarea0pt","std::vector<float>", &ljnoarea0pt);
   tree->Branch("ljnoarea0eta","std::vector<float>",&ljnoarea0eta);
   tree->Branch("ljnoarea0phi","std::vector<float>",&ljnoarea0phi);
   tree->Branch("ljnoarea5pt","std::vector<float>", &ljnoarea5pt);
   tree->Branch("ljnoarea5eta","std::vector<float>",&ljnoarea5eta);
   tree->Branch("ljnoarea5phi","std::vector<float>",&ljnoarea5phi);
   tree->Branch("tljnoarea0pt","std::vector<float>",&tljnoarea0pt);
   tree->Branch("tljnoarea5pt","std::vector<float>",&tljnoarea5pt);
   tree->Branch("tljnoarea0eta","std::vector<float>",&tljnoarea0eta);
   tree->Branch("tljnoarea5eta","std::vector<float>",&tljnoarea5eta);

   tree->Branch("truejetpt","std::vector<float>",&truejetpt);
   tree->Branch("truejeteta","std::vector<float>",&truejeteta);
   tree->Branch("truejetphi","std::vector<float>",&truejetphi);
   tree->Branch("truejetenergy","std::vector<float>",&truejetenergy);

   tree->Branch("truelargejetpt","std::vector<float>",    &truelargejetpt);
   tree->Branch("truelargejeteta","std::vector<float>",   &truelargejeteta);
   tree->Branch("truelargejetphi","std::vector<float>",   &truelargejetphi);
   tree->Branch("truelargejetenergy","std::vector<float>",&truelargejetenergy);
   
   if(Debug()) cout <<"Analysis_tree::AddBranches End" << endl;
   return;
 }

void Analysis_tree::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_tree::ResetBranches Start" << endl;

  // Event Info
  fTEventNumber           = -9999;
  fTRunNumber             = -9999;
  fTWeight                = -999.99;
  fTNPVtruth              = -99;
  fTNPV                   = -99;
  fTMu                    = -99;
  fTRho                    = -99;

  // jet vars
  cljvfcorr ->clear();  
  clfem     ->clear();  
  clcenterlambda     ->clear();  
  clpt     ->clear();  
  cleta    ->clear();  
  clphi    ->clear();  
  clenergy ->clear();  
  trkpt     ->clear();  
  trketa    ->clear();  
  trkphi    ->clear();  
  trkispv    ->clear();  

  j0ncl     ->clear();  
  j5ncl     ->clear();  
  j0mass     ->clear();  
  j5mass     ->clear();  
  j0width     ->clear();  
  j5width     ->clear();  
  j0pt     ->clear();  
  j0eta     ->clear();  
  j0phi     ->clear();  
  j0jvt     ->clear();  
  j5pt     ->clear();  
  j5eta     ->clear();  
  j5phi     ->clear();  
  tj0pt     ->clear();  
  tj5pt     ->clear();  
  tj0eta     ->clear();  
  tj5eta     ->clear();  

  lj0ncl     ->clear();  
  lj5ncl     ->clear();  
  lj0mass     ->clear();  
  lj5mass     ->clear();  
  lj0width     ->clear();  
  lj5width     ->clear();  
  lj0pt     ->clear();  
  lj0eta     ->clear();  
  lj0phi     ->clear();  
  lj5pt     ->clear();  
  lj5eta     ->clear();  
  lj5phi     ->clear();  
  tlj0pt     ->clear();  
  tlj5pt     ->clear();  
  tlj0eta     ->clear();  
  tlj5eta     ->clear();  

  jnoarea0ncl     ->clear();  
  jnoarea5ncl     ->clear();  
  jnoarea0mass     ->clear();  
  jnoarea5mass     ->clear();  
  jnoarea0width     ->clear();  
  jnoarea5width     ->clear();  
  jnoarea0pt     ->clear();  
  jnoarea0eta     ->clear();  
  jnoarea0phi     ->clear();  
  jnoarea5pt     ->clear();  
  jnoarea5eta     ->clear();  
  jnoarea5phi     ->clear();  
  tjnoarea0pt     ->clear();  
  tjnoarea5pt     ->clear();  
  tjnoarea0eta     ->clear();  
  tjnoarea5eta     ->clear();  

  ljnoarea0ncl     ->clear();  
  ljnoarea5ncl     ->clear();  
  ljnoarea0mass     ->clear();  
  ljnoarea5mass     ->clear();  
  ljnoarea0width     ->clear();  
  ljnoarea5width     ->clear();  
  ljnoarea0pt     ->clear();  
  ljnoarea0eta     ->clear();  
  ljnoarea0phi     ->clear();  
  ljnoarea5pt     ->clear();  
  ljnoarea5eta     ->clear();  
  ljnoarea5phi     ->clear();  
  tljnoarea0pt     ->clear();  
  tljnoarea5pt     ->clear();  
  tljnoarea0eta     ->clear();  
  tljnoarea5eta     ->clear();  

  truejetpt     ->clear();  
  truejeteta     ->clear();  
  truejetphi     ->clear();  
  truejetenergy     ->clear();  

  truelargejetpt     ->clear();  
  truelargejeteta     ->clear();  
  truelargejetphi     ->clear();  
  truelargejetenergy     ->clear();  

  if(Debug()) cout <<"Analysis_tree::ResetBranches End" << endl;
  return;
}

void Analysis_tree::FillEventVars(TTree *tree){
  if(Debug()) cout <<"Analysis_tree::FillEventVars Begin" << endl;

  // Event Info
  fTEventNumber                 = Int("EventNumber");
  fTRunNumber                   = Int("RunNumber");
  fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
  fTNPV                         = Exists("NPV")? Int("NPV"):-1;
  fTWeight                      = DefaultWeight();
  fTMu                          = Float("averageIntPerXing");
  fTRho                          = Float("Eventshape_rhoKt4LC");

  if(Debug()) cout <<"Analysis_tree::FillEventVars End" << endl;
  return;
}

void Analysis_tree::FillClVars(TTree* tree, int jindex, const MomKey Key){
  if(Debug()) cout <<"Analysis_tree::FillClVars Begin" << endl;

  Particle  *mycl         = &(cluster(jindex, Key));

  cljvfcorr->push_back(mycl->Float("corrJVF"));
  clfem->push_back(mycl->Float("fem"));
  clcenterlambda->push_back(mycl->Float("centerlambda"));
  clpt->push_back(mycl->p.Pt());
  cleta->push_back(mycl->p.Eta());
  clphi->push_back(mycl->p.Phi());
  clenergy->push_back(mycl->p.E());

  if(Debug()) cout <<"Analysis_tree::FillClVars End" << endl;
  return;
}

void Analysis_tree::FillTrkVars(TTree* tree, int jindex){
  if(Debug()) cout <<"Analysis_tree::FillTrkVars Begin" << endl;

  Particle  *mytrk         = &(track(jindex));

  trkpt ->push_back(mytrk->p.Pt());
  trketa->push_back(mytrk->p.Eta());
  trkphi->push_back(mytrk->p.Phi());
  trkispv->push_back(mytrk->Int("origin")==0);

  if(Debug()) cout <<"Analysis_tree::FillTrkVars End" << endl;
  return;
}

void Analysis_tree::FillJetVars(TTree* tree, int jindex,
				const MomKey jetkey, 
				vector<float> *jetpt, vector<float> *jeteta,
				vector<float> *jetphi, vector<float> *tjetpt, vector<float> *tjeteta,
				vector<int> *jetncl, vector<float> *jetm, vector<float> *jetw,
				MomKey truthjetkey){
  if(Debug()) cout <<"Analysis_tree::FillJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  jetpt ->push_back(myjet->p.Pt());
  jetncl->push_back(myjet->Objs("constituents"));
  if(myjet->Bool("isHSJet")){
    Particle *tjet = (Particle*) myjet->Obj(truthjetkey);
    tjetpt ->push_back(tjet->p.Pt());
    tjeteta ->push_back(tjet->p.Eta());
  }
  else{
    tjetpt ->push_back(-1);
    tjeteta ->push_back(-99);
  }

  jeteta->push_back(myjet->p.Eta());
  jetphi->push_back(myjet->p.Phi());
  jetm->push_back(myjet->p.M());
  jetw->push_back(myjet->Float("width"));

  if(Debug()) cout <<"Analysis_tree::FillTrkVars End" << endl;
  return;
}

void Analysis_tree::FillTrueJetVars(TTree* tree, int jindex,
				    const MomKey jetkey){
  if(Debug()) cout <<"Analysis_tree::FillTrueJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  truejetpt ->push_back(myjet->p.Pt());
  truejeteta ->push_back(myjet->p.Eta());
  truejetphi ->push_back(myjet->p.Phi());
  truejetenergy ->push_back(myjet->p.E());

  if(Debug()) cout <<"Analysis_tree::FillTrkVars End" << endl;
  return;
}

void Analysis_tree::FillTrueLargeJetVars(TTree* tree, int jindex,
				    const MomKey jetkey){
  if(Debug()) cout <<"Analysis_tree::FillTrueJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  truelargejetpt ->push_back(myjet->p.Pt());
  truelargejeteta ->push_back(myjet->p.Eta());
  truelargejetphi ->push_back(myjet->p.Phi());
  truelargejetenergy ->push_back(myjet->p.E());

  if(Debug()) cout <<"Analysis_tree::FillTrkVars End" << endl;
  return;
}
