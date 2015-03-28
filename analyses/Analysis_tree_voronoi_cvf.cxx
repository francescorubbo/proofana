/**************************************************************************
 **
 **   File:         Analysis_tree_voronoi_cvf.cxx
 **
 **   Description:  See header
 **
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_tree_voronoi_cvf_cxx

#include "Analysis_tree_voronoi_cvf.h"
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
 void Analysis_tree_voronoi_cvf::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_tree_voronoi_cvf: DEBUG In WorkerBegin() " << endl;

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
  
  j0cvfncl     = new vector<int>();  
  j0cvfmass    = new vector<float>();  
  j0cvfwidth   = new vector<float>();  
  j0cvfpt      = new vector<float>();  
  j0cvfeta     = new vector<float>();  
  j0cvfphi     = new vector<float>();  
  j0cvfjvt     = new vector<float>();  
  tj0cvfpt     = new vector<float>();  
  tj0cvfeta    = new vector<float>();  
  
  jnoarea0ncl     = new vector<int>();  
  jnoarea0mass    = new vector<float>();  
  jnoarea0width   = new vector<float>();  
  jnoarea0pt      = new vector<float>();  
  jnoarea0eta     = new vector<float>();  
  jnoarea0phi     = new vector<float>();  
  tjnoarea0pt     = new vector<float>();  
  tjnoarea0eta    = new vector<float>();  

  jvoro0ncl     = new vector<int>();  
  jvoro0mass    = new vector<float>();  
  jvoro0width   = new vector<float>();  
  jvoro0pt      = new vector<float>();  
  jvoro0eta     = new vector<float>();  
  jvoro0phi     = new vector<float>();  
  tjvoro0pt     = new vector<float>();  
  tjvoro0eta    = new vector<float>();  

  jvorocvfxsncl     = new vector<int>();  
  jvorocvfxsmass    = new vector<float>();  
  jvorocvfxswidth   = new vector<float>();  
  jvorocvfxspt      = new vector<float>();  
  jvorocvfxseta     = new vector<float>();  
  jvorocvfxsphi     = new vector<float>();  
  tjvorocvfxspt     = new vector<float>();  
  tjvorocvfxseta    = new vector<float>();  

  jvorocvf5sncl     = new vector<int>();  
  jvorocvf5smass    = new vector<float>();  
  jvorocvf5swidth   = new vector<float>();  
  jvorocvf5spt      = new vector<float>();  
  jvorocvf5seta     = new vector<float>();  
  jvorocvf5sphi     = new vector<float>();  
  tjvorocvf5spt     = new vector<float>();  
  tjvorocvf5seta    = new vector<float>();  

  jvorosncl     = new vector<int>();  
  jvorosmass    = new vector<float>();  
  jvoroswidth   = new vector<float>();  
  jvorospt      = new vector<float>();  
  jvoroseta     = new vector<float>();  
  jvorosphi     = new vector<float>();  
  tjvorospt     = new vector<float>();  
  tjvoroseta    = new vector<float>();  

  jvoro1ncl     = new vector<int>();  
  jvoro1mass    = new vector<float>();  
  jvoro1width   = new vector<float>();  
  jvoro1pt      = new vector<float>();  
  jvoro1eta     = new vector<float>();  
  jvoro1phi     = new vector<float>();  
  tjvoro1pt     = new vector<float>();  
  tjvoro1eta    = new vector<float>();  

  jvoro0cvf5ncl     = new vector<int>();  
  jvoro0cvf5mass    = new vector<float>();  
  jvoro0cvf5width   = new vector<float>();  
  jvoro0cvf5pt      = new vector<float>();  
  jvoro0cvf5eta     = new vector<float>();  
  jvoro0cvf5phi     = new vector<float>();  
  tjvoro0cvf5pt     = new vector<float>();  
  tjvoro0cvf5eta    = new vector<float>();  

  jvoro1cvf5ncl     = new vector<int>();  
  jvoro1cvf5mass    = new vector<float>();  
  jvoro1cvf5width   = new vector<float>();  
  jvoro1cvf5pt      = new vector<float>();  
  jvoro1cvf5eta     = new vector<float>();  
  jvoro1cvf5phi     = new vector<float>();  
  tjvoro1cvf5pt     = new vector<float>();  
  tjvoro1cvf5eta    = new vector<float>();  

  jvoro0cvfxncl     = new vector<int>();  
  jvoro0cvfxmass    = new vector<float>();  
  jvoro0cvfxwidth   = new vector<float>();  
  jvoro0cvfxpt      = new vector<float>();  
  jvoro0cvfxeta     = new vector<float>();  
  jvoro0cvfxphi     = new vector<float>();  
  tjvoro0cvfxpt     = new vector<float>();  
  tjvoro0cvfxeta    = new vector<float>();  

  jvoro1cvfxncl     = new vector<int>();  
  jvoro1cvfxmass    = new vector<float>();  
  jvoro1cvfxwidth   = new vector<float>();  
  jvoro1cvfxpt      = new vector<float>();  
  jvoro1cvfxeta     = new vector<float>();  
  jvoro1cvfxphi     = new vector<float>();  
  tjvoro1cvfxpt     = new vector<float>();  
  tjvoro1cvfxeta    = new vector<float>();  

  jvoro10ncl     = new vector<int>();  
  jvoro10mass    = new vector<float>();  
  jvoro10width   = new vector<float>();  
  jvoro10pt      = new vector<float>();  
  jvoro10eta     = new vector<float>();  
  jvoro10phi     = new vector<float>();  
  tjvoro10pt     = new vector<float>();  
  tjvoro10eta    = new vector<float>();  


  truejetpt      = new vector<float>();  
  truejeteta     = new vector<float>();  
  truejetphi     = new vector<float>();  
  truejetenergy  = new vector<float>();  

  if (Debug()) cout << "Analysis_tree_voronoi_cvf: DEBUG Finish WorkerBegin()" << endl;
}

///=========================================
/// ProcessEvent: run the analysis
///=========================================
bool Analysis_tree_voronoi_cvf::ProcessEvent()
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
void Analysis_tree_voronoi_cvf::WorkerTerminate()
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

  delete j0cvfncl    ;  
  delete j0cvfmass    ;  
  delete j0cvfwidth    ;  
  delete j0cvfpt    ;  
  delete j0cvfeta    ;  
  delete j0cvfphi    ;  
  delete j0cvfjvt    ;  
  delete tj0cvfpt    ;  
  delete tj0cvfeta    ;  

  delete jnoarea0ncl    ;  
  delete jnoarea0mass    ;  
  delete jnoarea0width    ;  
  delete jnoarea0pt    ;  
  delete jnoarea0eta    ;  
  delete jnoarea0phi    ;  
  delete tjnoarea0pt    ;  
  delete tjnoarea0eta    ;  

  delete jvoro0ncl    ;  
  delete jvoro0mass    ;  
  delete jvoro0width    ;  
  delete jvoro0pt    ;  
  delete jvoro0eta    ;  
  delete jvoro0phi    ;  
  delete tjvoro0pt    ;  
  delete tjvoro0eta    ;  

  delete jvorosncl    ;  
  delete jvorosmass    ;  
  delete jvoroswidth    ;  
  delete jvorospt    ;  
  delete jvoroseta    ;  
  delete jvorosphi    ;  
  delete tjvorospt    ;  
  delete tjvoroseta    ;  

  delete jvorocvf5sncl    ;  
  delete jvorocvf5smass    ;  
  delete jvorocvf5swidth    ;  
  delete jvorocvf5spt    ;  
  delete jvorocvf5seta    ;  
  delete jvorocvf5sphi    ;  
  delete tjvorocvf5spt    ;  
  delete tjvorocvf5seta    ;  

  delete jvorocvfxsncl    ;  
  delete jvorocvfxsmass    ;  
  delete jvorocvfxswidth    ;  
  delete jvorocvfxspt    ;  
  delete jvorocvfxseta    ;  
  delete jvorocvfxsphi    ;  
  delete tjvorocvfxspt    ;  
  delete tjvorocvfxseta    ;  

  delete jvoro1ncl    ;  
  delete jvoro1mass    ;  
  delete jvoro1width    ;  
  delete jvoro1pt    ;  
  delete jvoro1eta    ;  
  delete jvoro1phi    ;  
  delete tjvoro1pt    ;  
  delete tjvoro1eta    ;  

  delete jvoro0cvf5ncl    ;  
  delete jvoro0cvf5mass    ;  
  delete jvoro0cvf5width    ;  
  delete jvoro0cvf5pt    ;  
  delete jvoro0cvf5eta    ;  
  delete jvoro0cvf5phi    ;  
  delete tjvoro0cvf5pt    ;  
  delete tjvoro0cvf5eta    ;  

  delete jvoro1cvf5ncl    ;  
  delete jvoro1cvf5mass    ;  
  delete jvoro1cvf5width    ;  
  delete jvoro1cvf5pt    ;  
  delete jvoro1cvf5eta    ;  
  delete jvoro1cvf5phi    ;  
  delete tjvoro1cvf5pt    ;  
  delete tjvoro1cvf5eta    ;  

  delete jvoro0cvfxncl    ;  
  delete jvoro0cvfxmass    ;  
  delete jvoro0cvfxwidth    ;  
  delete jvoro0cvfxpt    ;  
  delete jvoro0cvfxeta    ;  
  delete jvoro0cvfxphi    ;  
  delete tjvoro0cvfxpt    ;  
  delete tjvoro0cvfxeta    ;  

  delete jvoro1cvfxncl    ;  
  delete jvoro1cvfxmass    ;  
  delete jvoro1cvfxwidth    ;  
  delete jvoro1cvfxpt    ;  
  delete jvoro1cvfxeta    ;  
  delete jvoro1cvfxphi    ;  
  delete tjvoro1cvfxpt    ;  
  delete tjvoro1cvfxeta    ;  


  delete jvoro10ncl    ;  
  delete jvoro10mass    ;  
  delete jvoro10width    ;  
  delete jvoro10pt    ;  
  delete jvoro10eta    ;  
  delete jvoro10phi    ;  
  delete tjvoro10pt    ;  
  delete tjvoro10eta    ;  

  delete truejetpt    ;  
  delete truejeteta    ;  
  delete truejetphi    ;  
  delete truejetenergy    ;  

}

void Analysis_tree_voronoi_cvf::FillTree(const MomKey Key){
  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillTree Start" << endl;
  ResetBranches(fTree);
  FillEventVars(fTree);
  for(int i=0; i<clusters(Key); ++i)
    FillClVars(fTree, i, Key);

  for(int i=0; i<jets("AntiKt4VoronoiZeroSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiZeroSigma",jvoro0pt,jvoro0eta,jvoro0phi,tjvoro0pt,tjvoro0eta,jvoro0ncl,jvoro0mass,jvoro0width);
  for(int i=0; i<jets("AntiKt4VoronoiOneSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiOneSigma",jvoro1pt,jvoro1eta,jvoro1phi,tjvoro1pt,tjvoro1eta,jvoro1ncl,jvoro1mass,jvoro1width);
  for(int i=0; i<jets("AntiKt4VoronoiTenSigma"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiTenSigma",jvoro10pt,jvoro10eta,jvoro10phi,tjvoro10pt,tjvoro10eta,jvoro10ncl,jvoro10mass,jvoro10width);

  for(int i=0; i<jets("AntiKt4LCTopoj0"); ++i)
    FillJetVars(fTree, i, "AntiKt4LCTopoj0",j0pt,j0eta,j0phi,tj0pt,tj0eta,j0ncl,j0mass,j0width);
  for(int i=0; i<jets("AntiKt4j0_CVF"); ++i)
    FillJetVars(fTree, i, "AntiKt4j0_CVF",j0cvfpt,j0cvfeta,j0cvfphi,tj0cvfpt,tj0cvfeta,j0cvfncl,j0cvfmass,j0cvfwidth);
  for(int i=0; i<jets("AntiKt4LCTopojnoarea0"); ++i)
    FillJetVars(fTree, i, "AntiKt4LCTopojnoarea0",jnoarea0pt,jnoarea0eta,jnoarea0phi,tjnoarea0pt,tjnoarea0eta,jnoarea0ncl,jnoarea0mass,jnoarea0width);

  for(int i=0; i<jets("AntiKt4VoronoiZeroSigma_CVF5"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiZeroSigma_CVF5",jvoro0cvf5pt,jvoro0cvf5eta,jvoro0cvf5phi,tjvoro0cvf5pt,tjvoro0cvf5eta,jvoro0cvf5ncl,jvoro0cvf5mass,jvoro0cvf5width);
  for(int i=0; i<jets("AntiKt4VoronoiOneSigma_CVF5"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiOneSigma_CVF5",jvoro1cvf5pt,jvoro1cvf5eta,jvoro1cvf5phi,tjvoro1cvf5pt,tjvoro1cvf5eta,jvoro1cvf5ncl,jvoro1cvf5mass,jvoro1cvf5width);
  for(int i=0; i<jets("AntiKt4VoronoiZeroSigma_CVFx"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiZeroSigma_CVFx",jvoro0cvfxpt,jvoro0cvfxeta,jvoro0cvfxphi,tjvoro0cvfxpt,tjvoro0cvfxeta,jvoro0cvfxncl,jvoro0cvfxmass,jvoro0cvfxwidth);
  for(int i=0; i<jets("AntiKt4VoronoiOneSigma_CVFx"); ++i)
    FillJetVars(fTree, i, "AntiKt4VoronoiOneSigma_CVFx",jvoro1cvfxpt,jvoro1cvfxeta,jvoro1cvfxphi,tjvoro1cvfxpt,tjvoro1cvfxeta,jvoro1cvfxncl,jvoro1cvfxmass,jvoro1cvfxwidth);

  for(int i=0; i<jets("AntiKt4Voronoi_SpreadPT"); ++i)
    FillJetVars(fTree, i, "AntiKt4Voronoi_SpreadPT",jvorospt,jvoroseta,jvorosphi,tjvorospt,tjvoroseta,jvorosncl,jvorosmass,jvoroswidth);
  for(int i=0; i<jets("AntiKt4Voronoi_CVF5_SpreadPT"); ++i)
    FillJetVars(fTree, i, "AntiKt4Voronoi_CVF5_SpreadPT",jvorocvf5spt,jvorocvf5seta,jvorocvf5sphi,tjvorocvf5spt,tjvorocvf5seta,jvorocvf5sncl,jvorocvf5smass,jvorocvf5swidth);
  for(int i=0; i<jets("AntiKt4Voronoi_CVFx_SpreadPT"); ++i)
    FillJetVars(fTree, i, "AntiKt4Voronoi_CVFx_SpreadPT",jvorocvfxspt,jvorocvfxseta,jvorocvfxsphi,tjvorocvfxspt,tjvorocvfxseta,jvorocvfxsncl,jvorocvfxsmass,jvorocvfxswidth);

  for(int i=0; i<jets("AntiKt4Truth"); ++i)
    FillTrueJetVars(fTree, i, "AntiKt4Truth");
  
  fTree->Fill();

  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillTree End" << endl;
  return;
}

 void Analysis_tree_voronoi_cvf::AddBranches(TTree *tree){
   if(Debug()) cout <<"Analysis_tree_voronoi_cvf::AddBranches Start" << endl;

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

   tree->Branch("j0cvfncl","std::vector<int>",  &j0cvfncl);
   tree->Branch("j0cvfmass","std::vector<float>",  &j0cvfmass);
   tree->Branch("j0cvfwidth","std::vector<float>",  &j0cvfwidth);
   tree->Branch("j0cvfpt","std::vector<float>", &j0cvfpt);
   tree->Branch("j0cvfeta","std::vector<float>",&j0cvfeta);
   tree->Branch("j0cvfphi","std::vector<float>",&j0cvfphi);
   tree->Branch("j0cvfjvt","std::vector<float>",&j0cvfjvt);
   tree->Branch("tj0cvfpt","std::vector<float>",&tj0cvfpt);
   tree->Branch("tj0cvfeta","std::vector<float>",&tj0cvfeta);

   tree->Branch("jnoarea0ncl","std::vector<int>",  &jnoarea0ncl);
   tree->Branch("jnoarea0mass","std::vector<float>",  &jnoarea0mass);
   tree->Branch("jnoarea0width","std::vector<float>",  &jnoarea0width);
   tree->Branch("jnoarea0pt","std::vector<float>", &jnoarea0pt);
   tree->Branch("jnoarea0eta","std::vector<float>",&jnoarea0eta);
   tree->Branch("jnoarea0phi","std::vector<float>",&jnoarea0phi);
   tree->Branch("tjnoarea0pt","std::vector<float>",&tjnoarea0pt);
   tree->Branch("tjnoarea0eta","std::vector<float>",&tjnoarea0eta);

   tree->Branch("jvoro0ncl","std::vector<int>",     &jvoro0ncl);
   tree->Branch("jvoro0mass","std::vector<float>",  &jvoro0mass);
   tree->Branch("jvoro0width","std::vector<float>", &jvoro0width);
   tree->Branch("jvoro0pt","std::vector<float>",    &jvoro0pt);
   tree->Branch("jvoro0eta","std::vector<float>",   &jvoro0eta);
   tree->Branch("jvoro0phi","std::vector<float>",   &jvoro0phi);
   tree->Branch("tjvoro0pt","std::vector<float>",      &tjvoro0pt);
   tree->Branch("tjvoro0eta","std::vector<float>",     &tjvoro0eta);

   tree->Branch("jvoro1ncl","std::vector<int>",     &jvoro1ncl);
   tree->Branch("jvoro1mass","std::vector<float>",  &jvoro1mass);
   tree->Branch("jvoro1width","std::vector<float>", &jvoro1width);
   tree->Branch("jvoro1pt","std::vector<float>",    &jvoro1pt);
   tree->Branch("jvoro1eta","std::vector<float>",   &jvoro1eta);
   tree->Branch("jvoro1phi","std::vector<float>",   &jvoro1phi);
   tree->Branch("tjvoro1pt","std::vector<float>",      &tjvoro1pt);
   tree->Branch("tjvoro1eta","std::vector<float>",     &tjvoro1eta);

   tree->Branch("jvoro0cvf5ncl","std::vector<int>",     &jvoro0cvf5ncl);
   tree->Branch("jvoro0cvf5mass","std::vector<float>",  &jvoro0cvf5mass);
   tree->Branch("jvoro0cvf5width","std::vector<float>", &jvoro0cvf5width);
   tree->Branch("jvoro0cvf5pt","std::vector<float>",    &jvoro0cvf5pt);
   tree->Branch("jvoro0cvf5eta","std::vector<float>",   &jvoro0cvf5eta);
   tree->Branch("jvoro0cvf5phi","std::vector<float>",   &jvoro0cvf5phi);
   tree->Branch("tjvoro0cvf5pt","std::vector<float>",      &tjvoro0cvf5pt);
   tree->Branch("tjvoro0cvf5eta","std::vector<float>",     &tjvoro0cvf5eta);

   tree->Branch("jvoro1cvf5ncl","std::vector<int>",     &jvoro1cvf5ncl);
   tree->Branch("jvoro1cvf5mass","std::vector<float>",  &jvoro1cvf5mass);
   tree->Branch("jvoro1cvf5width","std::vector<float>", &jvoro1cvf5width);
   tree->Branch("jvoro1cvf5pt","std::vector<float>",    &jvoro1cvf5pt);
   tree->Branch("jvoro1cvf5eta","std::vector<float>",   &jvoro1cvf5eta);
   tree->Branch("jvoro1cvf5phi","std::vector<float>",   &jvoro1cvf5phi);
   tree->Branch("tjvoro1cvf5pt","std::vector<float>",      &tjvoro1cvf5pt);
   tree->Branch("tjvoro1cvf5eta","std::vector<float>",     &tjvoro1cvf5eta);

   tree->Branch("jvoro0cvfxncl","std::vector<int>",     &jvoro0cvfxncl);
   tree->Branch("jvoro0cvfxmass","std::vector<float>",  &jvoro0cvfxmass);
   tree->Branch("jvoro0cvfxwidth","std::vector<float>", &jvoro0cvfxwidth);
   tree->Branch("jvoro0cvfxpt","std::vector<float>",    &jvoro0cvfxpt);
   tree->Branch("jvoro0cvfxeta","std::vector<float>",   &jvoro0cvfxeta);
   tree->Branch("jvoro0cvfxphi","std::vector<float>",   &jvoro0cvfxphi);
   tree->Branch("tjvoro0cvfxpt","std::vector<float>",      &tjvoro0cvfxpt);
   tree->Branch("tjvoro0cvfxeta","std::vector<float>",     &tjvoro0cvfxeta);

   tree->Branch("jvoro1cvfxncl","std::vector<int>",     &jvoro1cvfxncl);
   tree->Branch("jvoro1cvfxmass","std::vector<float>",  &jvoro1cvfxmass);
   tree->Branch("jvoro1cvfxwidth","std::vector<float>", &jvoro1cvfxwidth);
   tree->Branch("jvoro1cvfxpt","std::vector<float>",    &jvoro1cvfxpt);
   tree->Branch("jvoro1cvfxeta","std::vector<float>",   &jvoro1cvfxeta);
   tree->Branch("jvoro1cvfxphi","std::vector<float>",   &jvoro1cvfxphi);
   tree->Branch("tjvoro1cvfxpt","std::vector<float>",      &tjvoro1cvfxpt);
   tree->Branch("tjvoro1cvfxeta","std::vector<float>",     &tjvoro1cvfxeta);

   tree->Branch("jvorosncl","std::vector<int>",     &jvorosncl);
   tree->Branch("jvorosmass","std::vector<float>",  &jvorosmass);
   tree->Branch("jvoroswidth","std::vector<float>", &jvoroswidth);
   tree->Branch("jvorospt","std::vector<float>",    &jvorospt);
   tree->Branch("jvoroseta","std::vector<float>",   &jvoroseta);
   tree->Branch("jvorosphi","std::vector<float>",   &jvorosphi);
   tree->Branch("tjvorospt","std::vector<float>",      &tjvorospt);
   tree->Branch("tjvoroseta","std::vector<float>",     &tjvoroseta);

   tree->Branch("jvorocvfxsncl","std::vector<int>",     &jvorocvfxsncl);
   tree->Branch("jvorocvfxsmass","std::vector<float>",  &jvorocvfxsmass);
   tree->Branch("jvorocvfxswidth","std::vector<float>", &jvorocvfxswidth);
   tree->Branch("jvorocvfxspt","std::vector<float>",    &jvorocvfxspt);
   tree->Branch("jvorocvfxseta","std::vector<float>",   &jvorocvfxseta);
   tree->Branch("jvorocvfxsphi","std::vector<float>",   &jvorocvfxsphi);
   tree->Branch("tjvorocvfxspt","std::vector<float>",      &tjvorocvfxspt);
   tree->Branch("tjvorocvfxseta","std::vector<float>",     &tjvorocvfxseta);

   tree->Branch("jvorocvf5sncl","std::vector<int>",     &jvorocvf5sncl);
   tree->Branch("jvorocvf5smass","std::vector<float>",  &jvorocvf5smass);
   tree->Branch("jvorocvf5swidth","std::vector<float>", &jvorocvf5swidth);
   tree->Branch("jvorocvf5spt","std::vector<float>",    &jvorocvf5spt);
   tree->Branch("jvorocvf5seta","std::vector<float>",   &jvorocvf5seta);
   tree->Branch("jvorocvf5sphi","std::vector<float>",   &jvorocvf5sphi);
   tree->Branch("tjvorocvf5spt","std::vector<float>",      &tjvorocvf5spt);
   tree->Branch("tjvorocvf5seta","std::vector<float>",     &tjvorocvf5seta);

   tree->Branch("jvoro10ncl","std::vector<int>",     &jvoro10ncl);
   tree->Branch("jvoro10mass","std::vector<float>",  &jvoro10mass);
   tree->Branch("jvoro10width","std::vector<float>", &jvoro10width);
   tree->Branch("jvoro10pt","std::vector<float>",    &jvoro10pt);
   tree->Branch("jvoro10eta","std::vector<float>",   &jvoro10eta);
   tree->Branch("jvoro10phi","std::vector<float>",   &jvoro10phi);
   tree->Branch("tjvoro10pt","std::vector<float>",      &tjvoro10pt);
   tree->Branch("tjvoro10eta","std::vector<float>",     &tjvoro10eta);

   tree->Branch("truejetpt","std::vector<float>",&truejetpt);
   tree->Branch("truejeteta","std::vector<float>",&truejeteta);
   tree->Branch("truejetphi","std::vector<float>",&truejetphi);
   tree->Branch("truejetenergy","std::vector<float>",&truejetenergy);
   
   if(Debug()) cout <<"Analysis_tree_voronoi_cvf::AddBranches End" << endl;
   return;
 }

void Analysis_tree_voronoi_cvf::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::ResetBranches Start" << endl;

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

  j0cvfncl     ->clear();  
  j0cvfmass     ->clear();  
  j0cvfwidth     ->clear();  
  j0cvfpt     ->clear();  
  j0cvfeta     ->clear();  
  j0cvfphi     ->clear();  
  j0cvfjvt     ->clear();  
  tj0cvfpt     ->clear();  
  tj0cvfeta     ->clear();  

  jnoarea0ncl     ->clear();  
  jnoarea0mass     ->clear();  
  jnoarea0width     ->clear();  
  jnoarea0pt     ->clear();  
  jnoarea0eta     ->clear();  
  jnoarea0phi     ->clear();  
  tjnoarea0pt     ->clear();  
  tjnoarea0eta     ->clear();  

  jvoro0ncl     ->clear();  
  jvoro0mass     ->clear();  
  jvoro0width     ->clear();  
  jvoro0pt     ->clear();  
  jvoro0eta     ->clear();  
  jvoro0phi     ->clear();  
  tjvoro0pt     ->clear();  
  tjvoro0eta     ->clear();  

  jvorosncl     ->clear();  
  jvorosmass     ->clear();  
  jvoroswidth     ->clear();  
  jvorospt     ->clear();  
  jvoroseta     ->clear();  
  jvorosphi     ->clear();  
  tjvorospt     ->clear();  
  tjvoroseta     ->clear();  

  jvorocvf5sncl     ->clear();  
  jvorocvf5smass     ->clear();  
  jvorocvf5swidth     ->clear();  
  jvorocvf5spt     ->clear();  
  jvorocvf5seta     ->clear();  
  jvorocvf5sphi     ->clear();  
  tjvorocvf5spt     ->clear();  
  tjvorocvf5seta     ->clear();  

  jvorocvfxsncl     ->clear();  
  jvorocvfxsmass     ->clear();  
  jvorocvfxswidth     ->clear();  
  jvorocvfxspt     ->clear();  
  jvorocvfxseta     ->clear();  
  jvorocvfxsphi     ->clear();  
  tjvorocvfxspt     ->clear();  
  tjvorocvfxseta     ->clear();  

  jvoro1ncl     ->clear();  
  jvoro1mass     ->clear();  
  jvoro1width     ->clear();  
  jvoro1pt     ->clear();  
  jvoro1eta     ->clear();  
  jvoro1phi     ->clear();  
  tjvoro1pt     ->clear();  
  tjvoro1eta     ->clear();  

  jvoro0cvf5ncl     ->clear();  
  jvoro0cvf5mass     ->clear();  
  jvoro0cvf5width     ->clear();  
  jvoro0cvf5pt     ->clear();  
  jvoro0cvf5eta     ->clear();  
  jvoro0cvf5phi     ->clear();  
  tjvoro0cvf5pt     ->clear();  
  tjvoro0cvf5eta     ->clear();  

  jvoro1cvf5ncl     ->clear();  
  jvoro1cvf5mass     ->clear();  
  jvoro1cvf5width     ->clear();  
  jvoro1cvf5pt     ->clear();  
  jvoro1cvf5eta     ->clear();  
  jvoro1cvf5phi     ->clear();  
  tjvoro1cvf5pt     ->clear();  
  tjvoro1cvf5eta     ->clear();  

  jvoro0cvfxncl     ->clear();  
  jvoro0cvfxmass     ->clear();  
  jvoro0cvfxwidth     ->clear();  
  jvoro0cvfxpt     ->clear();  
  jvoro0cvfxeta     ->clear();  
  jvoro0cvfxphi     ->clear();  
  tjvoro0cvfxpt     ->clear();  
  tjvoro0cvfxeta     ->clear();  

  jvoro1cvfxncl     ->clear();  
  jvoro1cvfxmass     ->clear();  
  jvoro1cvfxwidth     ->clear();  
  jvoro1cvfxpt     ->clear();  
  jvoro1cvfxeta     ->clear();  
  jvoro1cvfxphi     ->clear();  
  tjvoro1cvfxpt     ->clear();  
  tjvoro1cvfxeta     ->clear();  

  jvoro10ncl     ->clear();  
  jvoro10mass     ->clear();  
  jvoro10width     ->clear();  
  jvoro10pt     ->clear();  
  jvoro10eta     ->clear();  
  jvoro10phi     ->clear();  
  tjvoro10pt     ->clear();  
  tjvoro10eta     ->clear();  

  truejetpt     ->clear();  
  truejeteta     ->clear();  
  truejetphi     ->clear();  
  truejetenergy     ->clear();  

  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::ResetBranches End" << endl;
  return;
}

void Analysis_tree_voronoi_cvf::FillEventVars(TTree *tree){
  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillEventVars Begin" << endl;

  // Event Info
  fTEventNumber                 = Int("EventNumber");
  fTRunNumber                   = Int("RunNumber");
  fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
  fTNPV                         = Exists("NPV")? Int("NPV"):-1;
  fTWeight                      = DefaultWeight();
  fTMu                          = Float("averageIntPerXing");
  fTRho                          = Float("rho_voronoi");
  fTSigma                          = Float("sigma_voronoi");

  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillEventVars End" << endl;
  return;
}

void Analysis_tree_voronoi_cvf::FillClVars(TTree* tree, int jindex, const MomKey Key){
  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillClVars Begin" << endl;

  Particle  *mycl         = &(cluster(jindex, Key));

  clpt->push_back(mycl->p.Pt());
  clarea->push_back(mycl->Float("area"));
  float orig_pt=mycl->p.Pt()+mycl->Float("area")*Float("rho_voronoi");
  clpt_orig->push_back(orig_pt);
  cleta->push_back(mycl->p.Eta());
  clphi->push_back(mycl->p.Phi());
  clenergy->push_back(mycl->p.E());

  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillClVars End" << endl;
  return;
}


void Analysis_tree_voronoi_cvf::FillJetVars(TTree* tree, int jindex,
				const MomKey jetkey, 
				vector<float> *jetpt, vector<float> *jeteta,
				vector<float> *jetphi, vector<float> *tjetpt, vector<float> *tjeteta,
				vector<int> *jetncl, vector<float> *jetm, vector<float> *jetw,
				MomKey truthjetkey){
  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillJetVars Begin" << endl;

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

  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillJetVars End" << endl;
  return;
}

void Analysis_tree_voronoi_cvf::FillTrueJetVars(TTree* tree, int jindex,
				    const MomKey jetkey){
  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillTrueJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  truejetpt ->push_back(myjet->p.Pt());
  truejeteta ->push_back(myjet->p.Eta());
  truejetphi ->push_back(myjet->p.Phi());
  truejetenergy ->push_back(myjet->p.E());

  if(Debug()) cout <<"Analysis_tree_voronoi_cvf::FillTrkVars End" << endl;
  return;
}

