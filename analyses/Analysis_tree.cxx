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
  j1ncl     = new vector<int>();  
  j2ncl     = new vector<int>();  
  j3ncl     = new vector<int>();  
  j4ncl     = new vector<int>();  
  j5ncl     = new vector<int>();  
  j0pt     = new vector<float>();  
  j0eta     = new vector<float>();  
  j0phi     = new vector<float>();  
  j1pt     = new vector<float>();  
  j1eta     = new vector<float>();  
  j1phi     = new vector<float>();  
  j2pt     = new vector<float>();  
  j2eta     = new vector<float>();  
  j2phi     = new vector<float>();  
  j3pt     = new vector<float>();  
  j3eta     = new vector<float>();  
  j3phi     = new vector<float>();  
  j4pt     = new vector<float>();  
  j4eta     = new vector<float>();  
  j4phi     = new vector<float>();  
  j5pt     = new vector<float>();  
  j5eta     = new vector<float>();  
  j5phi     = new vector<float>();  
  tj0pt     = new vector<float>();  
  tj1pt     = new vector<float>();  
  tj2pt     = new vector<float>();  
  tj3pt     = new vector<float>();  
  tj4pt     = new vector<float>();  
  tj5pt     = new vector<float>();  

  truejetpt     = new vector<float>();  
  truejeteta     = new vector<float>();  
  truejetphi     = new vector<float>();  
  truejetenergy     = new vector<float>();  

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
  delete j1ncl    ;  
  delete j2ncl    ;  
  delete j3ncl    ;  
  delete j4ncl    ;  
  delete j5ncl    ;  
  delete j0pt    ;  
  delete j0eta    ;  
  delete j0phi    ;  
  delete j1pt    ;  
  delete j1eta    ;  
  delete j1phi    ;  
  delete j2pt    ;  
  delete j2eta    ;  
  delete j2phi    ;  
  delete j3pt    ;  
  delete j3eta    ;  
  delete j3phi    ;  
  delete j4pt    ;  
  delete j4eta    ;  
  delete j4phi    ;  
  delete j5pt    ;  
  delete j5eta    ;  
  delete j5phi    ;  
  delete tj0pt    ;  
  delete tj1pt    ;  
  delete tj2pt    ;  
  delete tj3pt    ;  
  delete tj4pt    ;  
  delete tj5pt    ;  

  delete truejetpt    ;  
  delete truejeteta    ;  
  delete truejetphi    ;  
  delete truejetenergy    ;  

}

void Analysis_tree::FillTree(const MomKey Key){
  if(Debug()) cout <<"Analysis_tree::FillTree Start" << endl;
  ResetBranches(fTree);
  FillEventVars(fTree);
  for(int i=0; i<clusters(Key); ++i)
    FillClVars(fTree, i, Key);
  for(int i=0; i<tracks(); ++i)
    FillTrkVars(fTree, i);
  for(int i=0; i<jets("AntiKt4jvf0"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf0",j0pt,j0eta,j0phi,tj0pt,j0ncl);
  for(int i=0; i<jets("AntiKt4jvf1"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf1",j1pt,j1eta,j1phi,tj1pt,j1ncl);
  for(int i=0; i<jets("AntiKt4jvf2"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf2",j2pt,j2eta,j2phi,tj2pt,j2ncl);
  for(int i=0; i<jets("AntiKt4jvf3"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf3",j3pt,j3eta,j3phi,tj3pt,j3ncl);
  for(int i=0; i<jets("AntiKt4jvf4"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf4",j4pt,j4eta,j4phi,tj4pt,j4ncl);
  for(int i=0; i<jets("AntiKt4jvf5"); ++i)
    FillJetVars(fTree, i, "AntiKt4jvf5",j5pt,j5eta,j5phi,tj5pt,j5ncl);
  for(int i=0; i<jets("AntiKt4Truth"); ++i)
    FillTrueJetVars(fTree, i, "AntiKt4Truth");
  
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

   tree->Branch("j0ncl","std::vector<int>",&j0ncl);
   tree->Branch("j1ncl","std::vector<int>",&j1ncl);
   tree->Branch("j2ncl","std::vector<int>",&j2ncl);
   tree->Branch("j3ncl","std::vector<int>",&j3ncl);
   tree->Branch("j4ncl","std::vector<int>",&j4ncl);
   tree->Branch("j5ncl","std::vector<int>",&j5ncl);
   tree->Branch("j0pt","std::vector<float>",&j0pt);
   tree->Branch("j0eta","std::vector<float>",&j0eta);
   tree->Branch("j0phi","std::vector<float>",&j0phi);
   tree->Branch("j1pt","std::vector<float>",&j1pt);
   tree->Branch("j1eta","std::vector<float>",&j1eta);
   tree->Branch("j1phi","std::vector<float>",&j1phi);
   tree->Branch("j2pt","std::vector<float>",&j2pt);
   tree->Branch("j2eta","std::vector<float>",&j2eta);
   tree->Branch("j2phi","std::vector<float>",&j2phi);
   tree->Branch("j3pt","std::vector<float>",&j3pt);
   tree->Branch("j3eta","std::vector<float>",&j3eta);
   tree->Branch("j3phi","std::vector<float>",&j3phi);
   tree->Branch("j4pt","std::vector<float>",&j4pt);
   tree->Branch("j4eta","std::vector<float>",&j4eta);
   tree->Branch("j4phi","std::vector<float>",&j4phi);
   tree->Branch("j5pt","std::vector<float>",&j5pt);
   tree->Branch("j5eta","std::vector<float>",&j5eta);
   tree->Branch("j5phi","std::vector<float>",&j5phi);
   tree->Branch("tj0pt","std::vector<float>",&tj0pt);
   tree->Branch("tj1pt","std::vector<float>",&tj1pt);
   tree->Branch("tj2pt","std::vector<float>",&tj2pt);
   tree->Branch("tj3pt","std::vector<float>",&tj3pt);
   tree->Branch("tj4pt","std::vector<float>",&tj4pt);
   tree->Branch("tj5pt","std::vector<float>",&tj5pt);

   tree->Branch("truejetpt","std::vector<float>",&truejetpt);
   tree->Branch("truejeteta","std::vector<float>",&truejeteta);
   tree->Branch("truejetphi","std::vector<float>",&truejetphi);
   tree->Branch("truejetenergy","std::vector<float>",&truejetenergy);
   
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
  j1ncl     ->clear();  
  j2ncl     ->clear();  
  j3ncl     ->clear();  
  j4ncl     ->clear();  
  j5ncl     ->clear();  
  j0pt     ->clear();  
  j0eta     ->clear();  
  j0phi     ->clear();  
  j1pt     ->clear();  
  j1eta     ->clear();  
  j1phi     ->clear();  
  j2pt     ->clear();  
  j2eta     ->clear();  
  j2phi     ->clear();  
  j3pt     ->clear();  
  j3eta     ->clear();  
  j3phi     ->clear();  
  j4pt     ->clear();  
  j4eta     ->clear();  
  j4phi     ->clear();  
  j5pt     ->clear();  
  j5eta     ->clear();  
  j5phi     ->clear();  
  tj0pt     ->clear();  
  tj1pt     ->clear();  
  tj2pt     ->clear();  
  tj3pt     ->clear();  
  tj4pt     ->clear();  
  tj5pt     ->clear();  

  truejetpt     ->clear();  
  truejeteta     ->clear();  
  truejetphi     ->clear();  
  truejetenergy     ->clear();  

  if(Debug()) cout <<"Analysis_tree::ResetBranches End" << endl;
  return;
}

void Analysis_tree::FillEventVars(TTree *tree){
  if(Debug()) cout <<"Analysis_tree::FillEventVars Begin" << endl;

  // Event Info
  fTEventNumber                 = Int("EventNumber");
  fTRunNumber                   = Int("RunNumber");
  fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
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
				vector<float> *jetphi, vector<float> *tjetpt,
				vector<int> *jetncl){
  if(Debug()) cout <<"Analysis_tree::FillJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  float area = myjet->Float("area")/1000;
  jetpt ->push_back(myjet->p.Pt()-fTRho*area);
  jetncl->push_back(myjet->Objs("constituents"));
  if(myjet->Bool("isHSJet")){
    Particle *tjet = (Particle*) myjet->Obj("AntiKt4Truth_match");
    tjetpt ->push_back(tjet->p.Pt());}
  else
    tjetpt ->push_back(-1);

  jeteta->push_back(myjet->p.Eta());
  jetphi->push_back(myjet->p.Phi());

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
