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

  fTree = new TTree("tree", "tree");
  AddBranches(fTree);

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

}

void Analysis_tree::FillTree(const MomKey Key){
  if(Debug()) cout <<"Analysis_tree::FillTree Start" << endl;
  ResetBranches(fTree);
  FillEventVars(fTree);
  for(int i=0; i<clusters(Key); ++i)
    FillClVars(fTree, i, Key);
  for(int i=0; i<tracks(); ++i)
    FillTrkVars(fTree, i);
  
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

   // Jet vars ------------------------------------------------------------

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

  // jet vars
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

  if(Debug()) cout <<"Analysis_tree::FillEventVars End" << endl;
  return;
}

void Analysis_tree::FillClVars(TTree* tree, int jindex, const MomKey Key){
  if(Debug()) cout <<"Analysis_tree::FillClVars Begin" << endl;

  Particle  *mycl         = &(cluster(jindex, Key));

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
