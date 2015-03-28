/**************************************************************************
 **
 **   File:         Analysis_tree_voronoi.cxx
 **
 **   Description:  See header
 **
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_tree_voronoi_cxx

#include "Analysis_tree_voronoi.h"
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
 void Analysis_tree_voronoi::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_tree_voronoi: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  OutputDir()->cd();
  fTree = new TTree("tree", "tree");
  AddBranches(fTree);

  cljvfcorr   = new vector<float>();  
  clfem   = new vector<float>();  
  clcenterlambda   = new vector<float>();  
  clpt     = new vector<float>();  
  clpt_orig     = new vector<float>();  
  clarea     = new vector<float>();  
  cleta    = new vector<float>();  
  clphi    = new vector<float>();  
  clenergy = new vector<float>();  

  j0ncl     = new vector<int>();  
  j0mass    = new vector<float>();  
  j0width   = new vector<float>();  
  j0pt      = new vector<float>();  
  j0eta     = new vector<float>();  
  j0phi     = new vector<float>();  
  j0jvt     = new vector<float>();  
  tj0pt     = new vector<float>();  
  tj0eta    = new vector<float>();  
  
  lj0ncl     = new vector<int>();  
  lj0mass    = new vector<float>();  
  lj0width   = new vector<float>();  
  lj0pt      = new vector<float>();  
  lj0eta     = new vector<float>();  
  lj0phi     = new vector<float>();  
  tlj0pt     = new vector<float>();  
  tlj0eta    = new vector<float>();  

  jnoarea0ncl     = new vector<int>();  
  jnoarea0mass    = new vector<float>();  
  jnoarea0width   = new vector<float>();  
  jnoarea0pt      = new vector<float>();  
  jnoarea0eta     = new vector<float>();  
  jnoarea0phi     = new vector<float>();  
  tjnoarea0pt     = new vector<float>();  
  tjnoarea0eta    = new vector<float>();  

  ljnoarea0ncl     = new vector<int>();  
  ljnoarea0mass    = new vector<float>();  
  ljnoarea0width   = new vector<float>();  
  ljnoarea0pt      = new vector<float>();  
  ljnoarea0eta     = new vector<float>();  
  ljnoarea0phi     = new vector<float>();  
  tljnoarea0pt     = new vector<float>();  
  tljnoarea0eta    = new vector<float>();  

  jvoro0ncl     = new vector<int>();  
  jvoro0mass    = new vector<float>();  
  jvoro0width   = new vector<float>();  
  jvoro0pt      = new vector<float>();  
  jvoro0eta     = new vector<float>();  
  jvoro0phi     = new vector<float>();  
  tjvoro0pt     = new vector<float>();  
  tjvoro0eta    = new vector<float>();  

  ljvoro0ncl     = new vector<int>();  
  ljvoro0mass    = new vector<float>();  
  ljvoro0width   = new vector<float>();  
  ljvoro0pt      = new vector<float>();  
  ljvoro0eta     = new vector<float>();  
  ljvoro0phi     = new vector<float>();  
  tljvoro0pt     = new vector<float>();  
  tljvoro0eta    = new vector<float>();  

  jvoro1ncl     = new vector<int>();  
  jvoro1mass    = new vector<float>();  
  jvoro1width   = new vector<float>();  
  jvoro1pt      = new vector<float>();  
  jvoro1eta     = new vector<float>();  
  jvoro1phi     = new vector<float>();  
  tjvoro1pt     = new vector<float>();  
  tjvoro1eta    = new vector<float>();  

  ljvoro1ncl     = new vector<int>();  
  ljvoro1mass    = new vector<float>();  
  ljvoro1width   = new vector<float>();  
  ljvoro1pt      = new vector<float>();  
  ljvoro1eta     = new vector<float>();  
  ljvoro1phi     = new vector<float>();  
  tljvoro1pt     = new vector<float>();  
  tljvoro1eta    = new vector<float>();  

  jvoro2ncl     = new vector<int>();  
  jvoro2mass    = new vector<float>();  
  jvoro2width   = new vector<float>();  
  jvoro2pt      = new vector<float>();  
  jvoro2eta     = new vector<float>();  
  jvoro2phi     = new vector<float>();  
  tjvoro2pt     = new vector<float>();  
  tjvoro2eta    = new vector<float>();  

  ljvoro2ncl     = new vector<int>();  
  ljvoro2mass    = new vector<float>();  
  ljvoro2width   = new vector<float>();  
  ljvoro2pt      = new vector<float>();  
  ljvoro2eta     = new vector<float>();  
  ljvoro2phi     = new vector<float>();  
  tljvoro2pt     = new vector<float>();  
  tljvoro2eta    = new vector<float>();  

  jvoro3ncl     = new vector<int>();  
  jvoro3mass    = new vector<float>();  
  jvoro3width   = new vector<float>();  
  jvoro3pt      = new vector<float>();  
  jvoro3eta     = new vector<float>();  
  jvoro3phi     = new vector<float>();  
  tjvoro3pt     = new vector<float>();  
  tjvoro3eta    = new vector<float>();  

  ljvoro3ncl     = new vector<int>();  
  ljvoro3mass    = new vector<float>();  
  ljvoro3width   = new vector<float>();  
  ljvoro3pt      = new vector<float>();  
  ljvoro3eta     = new vector<float>();  
  ljvoro3phi     = new vector<float>();  
  tljvoro3pt     = new vector<float>();  
  tljvoro3eta    = new vector<float>();  

  jvoro4ncl     = new vector<int>();  
  jvoro4mass    = new vector<float>();  
  jvoro4width   = new vector<float>();  
  jvoro4pt      = new vector<float>();  
  jvoro4eta     = new vector<float>();  
  jvoro4phi     = new vector<float>();  
  tjvoro4pt     = new vector<float>();  
  tjvoro4eta    = new vector<float>();  

  ljvoro4ncl     = new vector<int>();  
  ljvoro4mass    = new vector<float>();  
  ljvoro4width   = new vector<float>();  
  ljvoro4pt      = new vector<float>();  
  ljvoro4eta     = new vector<float>();  
  ljvoro4phi     = new vector<float>();  
  tljvoro4pt     = new vector<float>();  
  tljvoro4eta    = new vector<float>();  

  truejetpt      = new vector<float>();  
  truejeteta     = new vector<float>();  
  truejetphi     = new vector<float>();  
  truejetenergy  = new vector<float>();  

  truelargejetpt      = new vector<float>();  
  truelargejeteta     = new vector<float>();  
  truelargejetphi     = new vector<float>();  
  truelargejetenergy  = new vector<float>();  

  if (Debug()) cout << "Analysis_tree_voronoi: DEBUG Finish WorkerBegin()" << endl;
}

///=========================================
/// ProcessEvent: run the analysis
///=========================================
bool Analysis_tree_voronoi::ProcessEvent()
{

  if (Debug()) cout << "Analysis_pileup: DEBUG In ProcessEvent(): RunNumber = " << RunNumber()
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();

  // -----------------------------
  // make objects truthsStable: stable truth particles for jets
  Analysis_JetMET_Base::AddStableParticles();

  FillTree("Voronoi");

  return true;
}

/// WorkerTerminate: clean up
///=========================================
void Analysis_tree_voronoi::WorkerTerminate()
{

  fTree->Write();
  // Nothing more
  delete cljvfcorr;  
  delete clfem    ;  
  delete clcenterlambda    ;  
  delete clpt    ;  
  delete clpt_orig    ;  
  delete clarea    ;  
  delete cleta   ;  
  delete clphi   ;  
  delete clenergy;  

  delete j0ncl    ;  
  delete j0mass    ;  
  delete j0width    ;  
  delete j0pt    ;  
  delete j0eta    ;  
  delete j0phi    ;  
  delete j0jvt    ;  
  delete tj0pt    ;  
  delete tj0eta    ;  

  delete lj0ncl    ;  
  delete lj0mass    ;  
  delete lj0width    ;  
  delete lj0pt    ;  
  delete lj0eta    ;  
  delete lj0phi    ;  
  delete tlj0pt    ;  
  delete tlj0eta    ;  

  delete jnoarea0ncl    ;  
  delete jnoarea0mass    ;  
  delete jnoarea0width    ;  
  delete jnoarea0pt    ;  
  delete jnoarea0eta    ;  
  delete jnoarea0phi    ;  
  delete tjnoarea0pt    ;  
  delete tjnoarea0eta    ;  

  delete ljnoarea0ncl    ;  
  delete ljnoarea0mass    ;  
  delete ljnoarea0width    ;  
  delete ljnoarea0pt    ;  
  delete ljnoarea0eta    ;  
  delete ljnoarea0phi    ;  
  delete tljnoarea0pt    ;  
  delete tljnoarea0eta    ;  

  delete jvoro0ncl    ;  
  delete jvoro0mass    ;  
  delete jvoro0width    ;  
  delete jvoro0pt    ;  
  delete jvoro0eta    ;  
  delete jvoro0phi    ;  
  delete tjvoro0pt    ;  
  delete tjvoro0eta    ;  

  delete ljvoro0ncl    ;  
  delete ljvoro0mass    ;  
  delete ljvoro0width    ;  
  delete ljvoro0pt    ;  
  delete ljvoro0eta    ;  
  delete ljvoro0phi    ;  
  delete tljvoro0pt    ;  
  delete tljvoro0eta    ;  

  delete jvoro1ncl    ;  
  delete jvoro1mass    ;  
  delete jvoro1width    ;  
  delete jvoro1pt    ;  
  delete jvoro1eta    ;  
  delete jvoro1phi    ;  
  delete tjvoro1pt    ;  
  delete tjvoro1eta    ;  

  delete ljvoro1ncl    ;  
  delete ljvoro1mass    ;  
  delete ljvoro1width    ;  
  delete ljvoro1pt    ;  
  delete ljvoro1eta    ;  
  delete ljvoro1phi    ;  
  delete tljvoro1pt    ;  
  delete tljvoro1eta    ;  

  delete jvoro2ncl    ;  
  delete jvoro2mass    ;  
  delete jvoro2width    ;  
  delete jvoro2pt    ;  
  delete jvoro2eta    ;  
  delete jvoro2phi    ;  
  delete tjvoro2pt    ;  
  delete tjvoro2eta    ;  

  delete ljvoro2ncl    ;  
  delete ljvoro2mass    ;  
  delete ljvoro2width    ;  
  delete ljvoro2pt    ;  
  delete ljvoro2eta    ;  
  delete ljvoro2phi    ;  
  delete tljvoro2pt    ;  
  delete tljvoro2eta    ;  

  delete jvoro3ncl    ;  
  delete jvoro3mass    ;  
  delete jvoro3width    ;  
  delete jvoro3pt    ;  
  delete jvoro3eta    ;  
  delete jvoro3phi    ;  
  delete tjvoro3pt    ;  
  delete tjvoro3eta    ;  

  delete ljvoro3ncl    ;  
  delete ljvoro3mass    ;  
  delete ljvoro3width    ;  
  delete ljvoro3pt    ;  
  delete ljvoro3eta    ;  
  delete ljvoro3phi    ;  
  delete tljvoro3pt    ;  
  delete tljvoro3eta    ;  

  delete jvoro4ncl    ;  
  delete jvoro4mass    ;  
  delete jvoro4width    ;  
  delete jvoro4pt    ;  
  delete jvoro4eta    ;  
  delete jvoro4phi    ;  
  delete tjvoro4pt    ;  
  delete tjvoro4eta    ;  

  delete ljvoro4ncl    ;  
  delete ljvoro4mass    ;  
  delete ljvoro4width    ;  
  delete ljvoro4pt    ;  
  delete ljvoro4eta    ;  
  delete ljvoro4phi    ;  
  delete tljvoro4pt    ;  
  delete tljvoro4eta    ;  

  delete truejetpt    ;  
  delete truejeteta    ;  
  delete truejetphi    ;  
  delete truejetenergy    ;  

  delete truelargejetpt    ;  
  delete truelargejeteta    ;  
  delete truelargejetphi    ;  
  delete truelargejetenergy    ;  

}

void Analysis_tree_voronoi::FillTree(const MomKey Key){
  if(Debug()) cout <<"Analysis_tree_voronoi::FillTree Start" << endl;
  ResetBranches(fTree);
  FillEventVars(fTree);
  for(int i=0; i<clusters(Key); ++i)
    FillClVars(fTree, i, Key);

  for(int i=0; i<jets("AntiKt4VoronoiZeroSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiZeroSigma",jvoro0pt,jvoro0eta,jvoro0phi,tjvoro0pt,tjvoro0eta,jvoro0ncl,jvoro0mass,jvoro0width);
  for(int i=0; i<jets("AntiKt10VoronoiZeroSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt10VoronoiZeroSigma",ljvoro0pt,ljvoro0eta,ljvoro0phi,tljvoro0pt,tljvoro0eta,ljvoro0ncl,ljvoro0mass,ljvoro0width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt4VoronoiOneSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiOneSigma",jvoro1pt,jvoro1eta,jvoro1phi,tjvoro1pt,tjvoro1eta,jvoro1ncl,jvoro1mass,jvoro1width);
  for(int i=0; i<jets("AntiKt10VoronoiOneSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt10VoronoiOneSigma",ljvoro1pt,ljvoro1eta,ljvoro1phi,tljvoro1pt,tljvoro1eta,ljvoro1ncl,ljvoro1mass,ljvoro1width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt4VoronoiTwoSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiTwoSigma",jvoro2pt,jvoro2eta,jvoro2phi,tjvoro2pt,tjvoro2eta,jvoro2ncl,jvoro2mass,jvoro2width);
  for(int i=0; i<jets("AntiKt10VoronoiTwoSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt10VoronoiTwoSigma",ljvoro2pt,ljvoro2eta,ljvoro2phi,tljvoro2pt,tljvoro2eta,ljvoro2ncl,ljvoro2mass,ljvoro2width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt4VoronoiThreeSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiThreeSigma",jvoro3pt,jvoro3eta,jvoro3phi,tjvoro3pt,tjvoro3eta,jvoro3ncl,jvoro3mass,jvoro3width);
  for(int i=0; i<jets("AntiKt10VoronoiThreeSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt10VoronoiThreeSigma",ljvoro3pt,ljvoro3eta,ljvoro3phi,tljvoro3pt,tljvoro3eta,ljvoro3ncl,ljvoro3mass,ljvoro3width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt4VoronoiFourSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiFourSigma",jvoro4pt,jvoro4eta,jvoro4phi,tjvoro4pt,tjvoro4eta,jvoro4ncl,jvoro4mass,jvoro4width);
  for(int i=0; i<jets("AntiKt10VoronoiFourSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt10VoronoiFourSigma",ljvoro4pt,ljvoro4eta,ljvoro4phi,tljvoro4pt,tljvoro4eta,ljvoro4ncl,ljvoro4mass,ljvoro4width,"AntiKt10Truth_match");
  
  for(int i=0; i<jets("AntiKt4LCTopoj0"); ++i)
    FillJetVars(fTree, i, "AntiKt4LCTopoj0",j0pt,j0eta,j0phi,tj0pt,tj0eta,j0ncl,j0mass,j0width);
  for(int i=0; i<jets("AntiKt10LCTopolj0"); ++i)
    FillJetVars(fTree, i, "AntiKt10LCTopolj0",lj0pt,lj0eta,lj0phi,tlj0pt,tlj0eta,lj0ncl,lj0mass,lj0width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt4LCTopojnoarea0"); ++i)
    FillJetVars(fTree, i, "AntiKt4LCTopojnoarea0",jnoarea0pt,jnoarea0eta,jnoarea0phi,tjnoarea0pt,tjnoarea0eta,jnoarea0ncl,jnoarea0mass,jnoarea0width);
  for(int i=0; i<jets("AntiKt10LCTopoljnoarea0"); ++i)
    FillJetVars(fTree, i, "AntiKt10LCTopoljnoarea0",ljnoarea0pt,ljnoarea0eta,ljnoarea0phi,tljnoarea0pt,tljnoarea0eta,ljnoarea0ncl,ljnoarea0mass,ljnoarea0width,"AntiKt10Truth_match");

  for(int i=0; i<jets("AntiKt4Truth"); ++i)
    FillTrueJetVars(fTree, i, "AntiKt4Truth");
  for(int i=0; i<jets("AntiKt10Truth"); ++i)
    FillTrueLargeJetVars(fTree, i, "AntiKt10Truth");
  
  fTree->Fill();

  if(Debug()) cout <<"Analysis_tree_voronoi::FillTree End" << endl;
  return;
}

 void Analysis_tree_voronoi::AddBranches(TTree *tree){
   if(Debug()) cout <<"Analysis_tree_voronoi::AddBranches Start" << endl;

   // Event Info
   tree->Branch("EventNumber", &fTEventNumber,            "EventNumber/I");
   tree->Branch("RunNumber",   &fTRunNumber,              "RunNumber/I");
   tree->Branch("Weight" ,     &fTWeight,                 "Weight/F");
   tree->Branch("NPVtruth" ,   &fTNPVtruth,               "NPVtruth/F");
   tree->Branch("NPV" ,   &fTNPV,               "NPV/F");
   tree->Branch("Mu" ,         &fTMu,                     "Mu/F");
   tree->Branch("Rho" ,         &fTRho,                     "Rho/F");
   tree->Branch("SigmaRho" ,         &fTSigma,                     "SigmaRho/F");

   // Jet vars ------------------------------------------------------------

   tree->Branch("cljvfcorr","std::vector<float>",&cljvfcorr);
   tree->Branch("clfem","std::vector<float>",&clfem);
   tree->Branch("clcenterlambda","std::vector<float>",&clcenterlambda);
   tree->Branch("clpt","std::vector<float>",&clpt);
   tree->Branch("clpt_orig","std::vector<float>",&clpt_orig);
   tree->Branch("clarea","std::vector<float>",&clarea);
   tree->Branch("cleta","std::vector<float>",&cleta);
   tree->Branch("clphi","std::vector<float>",&clphi);
   tree->Branch("clenergy","std::vector<float>",&clenergy);

   tree->Branch("j0ncl","std::vector<int>",  &j0ncl);
   tree->Branch("j0mass","std::vector<float>",  &j0mass);
   tree->Branch("j0width","std::vector<float>",  &j0width);
   tree->Branch("j0pt","std::vector<float>", &j0pt);
   tree->Branch("j0eta","std::vector<float>",&j0eta);
   tree->Branch("j0phi","std::vector<float>",&j0phi);
   tree->Branch("j0jvt","std::vector<float>",&j0jvt);
   tree->Branch("tj0pt","std::vector<float>",&tj0pt);
   tree->Branch("tj0eta","std::vector<float>",&tj0eta);

   tree->Branch("lj0ncl","std::vector<int>",  &lj0ncl);
   tree->Branch("lj0mass","std::vector<float>",  &lj0mass);
   tree->Branch("lj0width","std::vector<float>",  &lj0width);
   tree->Branch("lj0pt","std::vector<float>", &lj0pt);
   tree->Branch("lj0eta","std::vector<float>",&lj0eta);
   tree->Branch("lj0phi","std::vector<float>",&lj0phi);
   tree->Branch("tlj0pt","std::vector<float>",&tlj0pt);
   tree->Branch("tlj0eta","std::vector<float>",&tlj0eta);

   tree->Branch("jnoarea0ncl","std::vector<int>",  &jnoarea0ncl);
   tree->Branch("jnoarea0mass","std::vector<float>",  &jnoarea0mass);
   tree->Branch("jnoarea0width","std::vector<float>",  &jnoarea0width);
   tree->Branch("jnoarea0pt","std::vector<float>", &jnoarea0pt);
   tree->Branch("jnoarea0eta","std::vector<float>",&jnoarea0eta);
   tree->Branch("jnoarea0phi","std::vector<float>",&jnoarea0phi);
   tree->Branch("tjnoarea0pt","std::vector<float>",&tjnoarea0pt);
   tree->Branch("tjnoarea0eta","std::vector<float>",&tjnoarea0eta);

   tree->Branch("ljnoarea0ncl","std::vector<int>",  &ljnoarea0ncl);
   tree->Branch("ljnoarea0mass","std::vector<float>",  &ljnoarea0mass);
   tree->Branch("ljnoarea0width","std::vector<float>",  &ljnoarea0width);
   tree->Branch("ljnoarea0pt","std::vector<float>", &ljnoarea0pt);
   tree->Branch("ljnoarea0eta","std::vector<float>",&ljnoarea0eta);
   tree->Branch("ljnoarea0phi","std::vector<float>",&ljnoarea0phi);
   tree->Branch("tljnoarea0pt","std::vector<float>",&tljnoarea0pt);
   tree->Branch("tljnoarea0eta","std::vector<float>",&tljnoarea0eta);

   tree->Branch("jvoro0ncl","std::vector<int>",     &jvoro0ncl);
   tree->Branch("jvoro0mass","std::vector<float>",  &jvoro0mass);
   tree->Branch("jvoro0width","std::vector<float>", &jvoro0width);
   tree->Branch("jvoro0pt","std::vector<float>",    &jvoro0pt);
   tree->Branch("jvoro0eta","std::vector<float>",   &jvoro0eta);
   tree->Branch("jvoro0phi","std::vector<float>",   &jvoro0phi);
   tree->Branch("tjvoro0pt","std::vector<float>",      &tjvoro0pt);
   tree->Branch("tjvoro0eta","std::vector<float>",     &tjvoro0eta);

   tree->Branch("ljvoro0ncl","std::vector<int>",     &ljvoro0ncl);
   tree->Branch("ljvoro0mass","std::vector<float>",  &ljvoro0mass);
   tree->Branch("ljvoro0width","std::vector<float>", &ljvoro0width);
   tree->Branch("ljvoro0pt","std::vector<float>",    &ljvoro0pt);
   tree->Branch("ljvoro0eta","std::vector<float>",   &ljvoro0eta);
   tree->Branch("ljvoro0phi","std::vector<float>",   &ljvoro0phi);
   tree->Branch("tljvoro0pt","std::vector<float>",      &tljvoro0pt);
   tree->Branch("tljvoro0eta","std::vector<float>",     &tljvoro0eta);

   tree->Branch("jvoro1ncl","std::vector<int>",     &jvoro1ncl);
   tree->Branch("jvoro1mass","std::vector<float>",  &jvoro1mass);
   tree->Branch("jvoro1width","std::vector<float>", &jvoro1width);
   tree->Branch("jvoro1pt","std::vector<float>",    &jvoro1pt);
   tree->Branch("jvoro1eta","std::vector<float>",   &jvoro1eta);
   tree->Branch("jvoro1phi","std::vector<float>",   &jvoro1phi);
   tree->Branch("tjvoro1pt","std::vector<float>",      &tjvoro1pt);
   tree->Branch("tjvoro1eta","std::vector<float>",     &tjvoro1eta);

   tree->Branch("ljvoro1ncl","std::vector<int>",     &ljvoro1ncl);
   tree->Branch("ljvoro1mass","std::vector<float>",  &ljvoro1mass);
   tree->Branch("ljvoro1width","std::vector<float>", &ljvoro1width);
   tree->Branch("ljvoro1pt","std::vector<float>",    &ljvoro1pt);
   tree->Branch("ljvoro1eta","std::vector<float>",   &ljvoro1eta);
   tree->Branch("ljvoro1phi","std::vector<float>",   &ljvoro1phi);
   tree->Branch("tljvoro1pt","std::vector<float>",      &tljvoro1pt);
   tree->Branch("tljvoro1eta","std::vector<float>",     &tljvoro1eta);

   tree->Branch("jvoro2ncl","std::vector<int>",     &jvoro2ncl);
   tree->Branch("jvoro2mass","std::vector<float>",  &jvoro2mass);
   tree->Branch("jvoro2width","std::vector<float>", &jvoro2width);
   tree->Branch("jvoro2pt","std::vector<float>",    &jvoro2pt);
   tree->Branch("jvoro2eta","std::vector<float>",   &jvoro2eta);
   tree->Branch("jvoro2phi","std::vector<float>",   &jvoro2phi);
   tree->Branch("tjvoro2pt","std::vector<float>",      &tjvoro2pt);
   tree->Branch("tjvoro2eta","std::vector<float>",     &tjvoro2eta);

   tree->Branch("ljvoro2ncl","std::vector<int>",     &ljvoro2ncl);
   tree->Branch("ljvoro2mass","std::vector<float>",  &ljvoro2mass);
   tree->Branch("ljvoro2width","std::vector<float>", &ljvoro2width);
   tree->Branch("ljvoro2pt","std::vector<float>",    &ljvoro2pt);
   tree->Branch("ljvoro2eta","std::vector<float>",   &ljvoro2eta);
   tree->Branch("ljvoro2phi","std::vector<float>",   &ljvoro2phi);
   tree->Branch("tljvoro2pt","std::vector<float>",      &tljvoro2pt);
   tree->Branch("tljvoro2eta","std::vector<float>",     &tljvoro2eta);

   tree->Branch("jvoro3ncl","std::vector<int>",     &jvoro3ncl);
   tree->Branch("jvoro3mass","std::vector<float>",  &jvoro3mass);
   tree->Branch("jvoro3width","std::vector<float>", &jvoro3width);
   tree->Branch("jvoro3pt","std::vector<float>",    &jvoro3pt);
   tree->Branch("jvoro3eta","std::vector<float>",   &jvoro3eta);
   tree->Branch("jvoro3phi","std::vector<float>",   &jvoro3phi);
   tree->Branch("tjvoro3pt","std::vector<float>",      &tjvoro3pt);
   tree->Branch("tjvoro3eta","std::vector<float>",     &tjvoro3eta);

   tree->Branch("ljvoro3ncl","std::vector<int>",     &ljvoro3ncl);
   tree->Branch("ljvoro3mass","std::vector<float>",  &ljvoro3mass);
   tree->Branch("ljvoro3width","std::vector<float>", &ljvoro3width);
   tree->Branch("ljvoro3pt","std::vector<float>",    &ljvoro3pt);
   tree->Branch("ljvoro3eta","std::vector<float>",   &ljvoro3eta);
   tree->Branch("ljvoro3phi","std::vector<float>",   &ljvoro3phi);
   tree->Branch("tljvoro3pt","std::vector<float>",      &tljvoro3pt);
   tree->Branch("tljvoro3eta","std::vector<float>",     &tljvoro3eta);

   tree->Branch("jvoro4ncl","std::vector<int>",     &jvoro4ncl);
   tree->Branch("jvoro4mass","std::vector<float>",  &jvoro4mass);
   tree->Branch("jvoro4width","std::vector<float>", &jvoro4width);
   tree->Branch("jvoro4pt","std::vector<float>",    &jvoro4pt);
   tree->Branch("jvoro4eta","std::vector<float>",   &jvoro4eta);
   tree->Branch("jvoro4phi","std::vector<float>",   &jvoro4phi);
   tree->Branch("tjvoro4pt","std::vector<float>",      &tjvoro4pt);
   tree->Branch("tjvoro4eta","std::vector<float>",     &tjvoro4eta);

   tree->Branch("ljvoro4ncl","std::vector<int>",     &ljvoro4ncl);
   tree->Branch("ljvoro4mass","std::vector<float>",  &ljvoro4mass);
   tree->Branch("ljvoro4width","std::vector<float>", &ljvoro4width);
   tree->Branch("ljvoro4pt","std::vector<float>",    &ljvoro4pt);
   tree->Branch("ljvoro4eta","std::vector<float>",   &ljvoro4eta);
   tree->Branch("ljvoro4phi","std::vector<float>",   &ljvoro4phi);
   tree->Branch("tljvoro4pt","std::vector<float>",      &tljvoro4pt);
   tree->Branch("tljvoro4eta","std::vector<float>",     &tljvoro4eta);

   tree->Branch("truejetpt","std::vector<float>",&truejetpt);
   tree->Branch("truejeteta","std::vector<float>",&truejeteta);
   tree->Branch("truejetphi","std::vector<float>",&truejetphi);
   tree->Branch("truejetenergy","std::vector<float>",&truejetenergy);

   tree->Branch("truelargejetpt","std::vector<float>",    &truelargejetpt);
   tree->Branch("truelargejeteta","std::vector<float>",   &truelargejeteta);
   tree->Branch("truelargejetphi","std::vector<float>",   &truelargejetphi);
   tree->Branch("truelargejetenergy","std::vector<float>",&truelargejetenergy);
   
   if(Debug()) cout <<"Analysis_tree_voronoi::AddBranches End" << endl;
   return;
 }

void Analysis_tree_voronoi::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_tree_voronoi::ResetBranches Start" << endl;

  // Event Info
  fTEventNumber           = -9999;
  fTRunNumber             = -9999;
  fTWeight                = -999.99;
  fTNPVtruth              = -99;
  fTNPV                   = -99;
  fTMu                    = -99;
  fTRho                    = -99;
  fTSigma                    = -99;

  // jet vars
  cljvfcorr ->clear();  
  clfem     ->clear();  
  clcenterlambda     ->clear();  
  clpt     ->clear();  
  clpt_orig     ->clear();  
  clarea     ->clear();  
  cleta    ->clear();  
  clphi    ->clear();  
  clenergy ->clear();  

  j0ncl     ->clear();  
  j0mass     ->clear();  
  j0width     ->clear();  
  j0pt     ->clear();  
  j0eta     ->clear();  
  j0phi     ->clear();  
  j0jvt     ->clear();  
  tj0pt     ->clear();  
  tj0eta     ->clear();  

  lj0ncl     ->clear();  
  lj0mass     ->clear();  
  lj0width     ->clear();  
  lj0pt     ->clear();  
  lj0eta     ->clear();  
  lj0phi     ->clear();  
  tlj0pt     ->clear();  
  tlj0eta     ->clear();  

  jnoarea0ncl     ->clear();  
  jnoarea0mass     ->clear();  
  jnoarea0width     ->clear();  
  jnoarea0pt     ->clear();  
  jnoarea0eta     ->clear();  
  jnoarea0phi     ->clear();  
  tjnoarea0pt     ->clear();  
  tjnoarea0eta     ->clear();  

  ljnoarea0ncl     ->clear();  
  ljnoarea0mass     ->clear();  
  ljnoarea0width     ->clear();  
  ljnoarea0pt     ->clear();  
  ljnoarea0eta     ->clear();  
  ljnoarea0phi     ->clear();  
  tljnoarea0pt     ->clear();  
  tljnoarea0eta     ->clear();  

  jvoro0ncl     ->clear();  
  jvoro0mass     ->clear();  
  jvoro0width     ->clear();  
  jvoro0pt     ->clear();  
  jvoro0eta     ->clear();  
  jvoro0phi     ->clear();  
  tjvoro0pt     ->clear();  
  tjvoro0eta     ->clear();  

  ljvoro0ncl     ->clear();  
  ljvoro0mass     ->clear();  
  ljvoro0width     ->clear();  
  ljvoro0pt     ->clear();  
  ljvoro0eta     ->clear();  
  ljvoro0phi     ->clear();  
  tljvoro0pt     ->clear();  
  tljvoro0eta     ->clear();  

  jvoro1ncl     ->clear();  
  jvoro1mass     ->clear();  
  jvoro1width     ->clear();  
  jvoro1pt     ->clear();  
  jvoro1eta     ->clear();  
  jvoro1phi     ->clear();  
  tjvoro1pt     ->clear();  
  tjvoro1eta     ->clear();  

  ljvoro1ncl     ->clear();  
  ljvoro1mass     ->clear();  
  ljvoro1width     ->clear();  
  ljvoro1pt     ->clear();  
  ljvoro1eta     ->clear();  
  ljvoro1phi     ->clear();  
  tljvoro1pt     ->clear();  
  tljvoro1eta     ->clear();  

  jvoro2ncl     ->clear();  
  jvoro2mass     ->clear();  
  jvoro2width     ->clear();  
  jvoro2pt     ->clear();  
  jvoro2eta     ->clear();  
  jvoro2phi     ->clear();  
  tjvoro2pt     ->clear();  
  tjvoro2eta     ->clear();  

  ljvoro2ncl     ->clear();  
  ljvoro2mass     ->clear();  
  ljvoro2width     ->clear();  
  ljvoro2pt     ->clear();  
  ljvoro2eta     ->clear();  
  ljvoro2phi     ->clear();  
  tljvoro2pt     ->clear();  
  tljvoro2eta     ->clear();  

  jvoro3ncl     ->clear();  
  jvoro3mass     ->clear();  
  jvoro3width     ->clear();  
  jvoro3pt     ->clear();  
  jvoro3eta     ->clear();  
  jvoro3phi     ->clear();  
  tjvoro3pt     ->clear();  
  tjvoro3eta     ->clear();  

  ljvoro3ncl     ->clear();  
  ljvoro3mass     ->clear();  
  ljvoro3width     ->clear();  
  ljvoro3pt     ->clear();  
  ljvoro3eta     ->clear();  
  ljvoro3phi     ->clear();  
  tljvoro3pt     ->clear();  
  tljvoro3eta     ->clear();  

  jvoro4ncl     ->clear();  
  jvoro4mass     ->clear();  
  jvoro4width     ->clear();  
  jvoro4pt     ->clear();  
  jvoro4eta     ->clear();  
  jvoro4phi     ->clear();  
  tjvoro4pt     ->clear();  
  tjvoro4eta     ->clear();  

  ljvoro4ncl     ->clear();  
  ljvoro4mass     ->clear();  
  ljvoro4width     ->clear();  
  ljvoro4pt     ->clear();  
  ljvoro4eta     ->clear();  
  ljvoro4phi     ->clear();  
  tljvoro4pt     ->clear();  
  tljvoro4eta     ->clear();  

  truejetpt     ->clear();  
  truejeteta     ->clear();  
  truejetphi     ->clear();  
  truejetenergy     ->clear();  

  truelargejetpt     ->clear();  
  truelargejeteta     ->clear();  
  truelargejetphi     ->clear();  
  truelargejetenergy     ->clear();  

  if(Debug()) cout <<"Analysis_tree_voronoi::ResetBranches End" << endl;
  return;
}

void Analysis_tree_voronoi::FillEventVars(TTree *tree){
  if(Debug()) cout <<"Analysis_tree_voronoi::FillEventVars Begin" << endl;

  // Event Info
  fTEventNumber                 = Int("EventNumber");
  fTRunNumber                   = Int("RunNumber");
  fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
  fTNPV                         = Exists("NPV")? Int("NPV"):-1;
  fTWeight                      = DefaultWeight();
  fTMu                          = Float("averageIntPerXing");
  fTRho                          = Float("rho_voronoi");
  fTSigma                          = Float("sigma_voronoi");

  if(Debug()) cout <<"Analysis_tree_voronoi::FillEventVars End" << endl;
  return;
}

void Analysis_tree_voronoi::FillClVars(TTree* tree, int jindex, const MomKey Key){
  if(Debug()) cout <<"Analysis_tree_voronoi::FillClVars Begin" << endl;

  Particle  *mycl         = &(cluster(jindex, Key));

  clpt->push_back(mycl->p.Pt());
  clarea->push_back(mycl->Float("area"));
  float orig_pt=mycl->p.Pt()+mycl->Float("area")*Float("rho_voronoi");
  clpt_orig->push_back(orig_pt);
  cleta->push_back(mycl->p.Eta());
  clphi->push_back(mycl->p.Phi());
  clenergy->push_back(mycl->p.E());

  if(Debug()) cout <<"Analysis_tree_voronoi::FillClVars End" << endl;
  return;
}


void Analysis_tree_voronoi::FillJetVars(TTree* tree, int jindex,
				const MomKey jetkey, 
				vector<float> *jetpt, vector<float> *jeteta,
				vector<float> *jetphi, vector<float> *tjetpt, vector<float> *tjeteta,
				vector<int> *jetncl, vector<float> *jetm, vector<float> *jetw,
				MomKey truthjetkey){
  if(Debug()) cout <<"Analysis_tree_voronoi::FillJetVars Begin" << endl;

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

  if(Debug()) cout <<"Analysis_tree_voronoi::FillJetVars End" << endl;
  return;
}

void Analysis_tree_voronoi::FillTrueJetVars(TTree* tree, int jindex,
				    const MomKey jetkey){
  if(Debug()) cout <<"Analysis_tree_voronoi::FillTrueJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  truejetpt ->push_back(myjet->p.Pt());
  truejeteta ->push_back(myjet->p.Eta());
  truejetphi ->push_back(myjet->p.Phi());
  truejetenergy ->push_back(myjet->p.E());

  if(Debug()) cout <<"Analysis_tree_voronoi::FillTrkVars End" << endl;
  return;
}

void Analysis_tree_voronoi::FillTrueLargeJetVars(TTree* tree, int jindex,
				    const MomKey jetkey){
  if(Debug()) cout <<"Analysis_tree_voronoi::FillTrueJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  truelargejetpt ->push_back(myjet->p.Pt());
  truelargejeteta ->push_back(myjet->p.Eta());
  truelargejetphi ->push_back(myjet->p.Phi());
  truelargejetenergy ->push_back(myjet->p.E());

  if(Debug()) cout <<"Analysis_tree_voronoi::FillTrkVars End" << endl;
  return;
}
