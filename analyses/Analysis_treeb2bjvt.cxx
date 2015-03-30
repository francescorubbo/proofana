/**************************************************************************
 **
 **   File:         Analysis_treeb2bjvt.cxx
 **
 **   Description:  See header
 **
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_treeb2bjvt_cxx

#include "Analysis_treeb2bjvt.h"
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
 void Analysis_treeb2bjvt::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_treeb2bjvt: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  OutputDir()->cd();
  fTree = new TTree("tree", "tree");
  AddBranches(fTree);

  jpt     = new vector<float>();  
  jeta     = new vector<float>();  
  jphi     = new vector<float>();
  jjvt     = new vector<float>();
  jcorrjvf = new vector<float>();
  jrpt     = new vector<float>();
  jb2bpt  = new vector<float>();
  jb2bjvt  = new vector<float>();
  jb2bjvt_70  = new vector<float>();
  jb2bjvt_75  = new vector<float>();
  jb2bjvt_85  = new vector<float>();
  jb2bjvt_90  = new vector<float>();
  jb2bcorrjvf  = new vector<float>();
  jb2brpt  = new vector<float>();
  jb2bpt_trk2  = new vector<float>();
  jb2bjvt_trk2  = new vector<float>();
  jb2bcorrjvf_trk2  = new vector<float>();
  jb2brpt_trk2  = new vector<float>();
  jb2bpt_trk4  = new vector<float>();
  jb2bjvt_trk4  = new vector<float>();
  jb2bcorrjvf_trk4  = new vector<float>();
  jb2brpt_trk4  = new vector<float>();
  tjpt     = new vector<float>();  
  jishs     = new vector<bool>();  
  jispu     = new vector<bool>();  

  trk2jpt     = new vector<float>();  
  trk2jeta     = new vector<float>();  
  trk2jphi     = new vector<float>();
  trk2jjvt     = new vector<float>();
  trk2jcorrjvf = new vector<float>();
  trk2jrpt     = new vector<float>();
  ttrk2jpt     = new vector<float>();  
  trk2jishs     = new vector<bool>();  
  trk2jispu     = new vector<bool>();  

  trk4jpt     = new vector<float>();  
  trk4jeta     = new vector<float>();  
  trk4jphi     = new vector<float>();
  trk4jjvt     = new vector<float>();
  trk4jcorrjvf = new vector<float>();
  trk4jrpt     = new vector<float>();
  ttrk4jpt     = new vector<float>();  
  trk4jishs     = new vector<bool>();  
  trk4jispu     = new vector<bool>();  

  truejetpt     = new vector<float>();  
  truejeteta     = new vector<float>();  
  truejetphi     = new vector<float>();  
  truejetenergy     = new vector<float>();  

  if (Debug()) cout << "Analysis_treeb2bjvt: DEBUG Finish WorkerBegin()" << endl;
}

///=========================================
/// ProcessEvent: run the analysis
///=========================================
bool Analysis_treeb2bjvt::ProcessEvent()
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
void Analysis_treeb2bjvt::WorkerTerminate()
{

  fTree->Write();
  // Nothing more
  delete jpt    ;  
  delete jeta    ;  
  delete jphi    ;  
  delete jjvt    ;  
  delete jrpt    ;  
  delete jcorrjvf    ;  
  delete jb2bpt ;  
  delete jb2bjvt ;  
  delete jb2bjvt_70 ;  
  delete jb2bjvt_75 ;  
  delete jb2bjvt_85 ;  
  delete jb2bjvt_90 ;  
  delete jb2bcorrjvf ;  
  delete jb2brpt ;  
  delete jb2bpt_trk2 ;  
  delete jb2bjvt_trk2 ;  
  delete jb2bcorrjvf_trk2 ;  
  delete jb2brpt_trk2 ;  
  delete jb2bpt_trk4 ;  
  delete jb2bjvt_trk4 ;  
  delete jb2bcorrjvf_trk4 ;  
  delete jb2brpt_trk4 ;  
  delete tjpt    ;  
  delete jishs    ;  
  delete jispu   ;  

  delete trk2jpt    ;  
  delete trk2jeta    ;  
  delete trk2jphi    ;  
  delete trk2jjvt    ;  
  delete trk2jrpt    ;  
  delete trk2jcorrjvf    ;  
  delete ttrk2jpt    ;  
  delete trk2jishs    ;  
  delete trk2jispu   ;  

  delete trk4jpt    ;  
  delete trk4jeta    ;  
  delete trk4jphi    ;  
  delete trk4jjvt    ;  
  delete trk4jrpt    ;  
  delete trk4jcorrjvf    ;  
  delete ttrk4jpt    ;  
  delete trk4jishs    ;  
  delete trk4jispu   ;  

  delete truejetpt    ;  
  delete truejeteta    ;  
  delete truejetphi    ;  
  delete truejetenergy    ;  

}

void Analysis_treeb2bjvt::FillTree(const MomKey Key){
  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTree Start" << endl;
  ResetBranches(fTree);
  FillEventVars(fTree);
  for(int i=0; i<jets("AntiKt4LCTopo"); ++i)
    FillJetVars(fTree, i, "AntiKt4LCTopo");
  for(int i=0; i<jets("AntiKt2tracks"); ++i)
    FillTrk2JetVars(fTree, i, "AntiKt2tracks");
  for(int i=0; i<jets("AntiKt4tracks"); ++i)
    FillTrk4JetVars(fTree, i, "AntiKt4tracks");
  for(int i=0; i<jets("AntiKt4Truth"); ++i)
    FillTrueJetVars(fTree, i, "AntiKt4Truth");
  
  fTree->Fill();

  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTree End" << endl;
  return;
}

 void Analysis_treeb2bjvt::AddBranches(TTree *tree){
   if(Debug()) cout <<"Analysis_treeb2bjvt::AddBranches Start" << endl;

   // Event Info
   tree->Branch("EventNumber", &fTEventNumber,            "EventNumber/I");
   tree->Branch("RunNumber",   &fTRunNumber,              "RunNumber/I");
   tree->Branch("Weight" ,     &fTWeight,                 "Weight/F");
   tree->Branch("NPVtruth" ,   &fTNPVtruth,               "NPVtruth/F");
   tree->Branch("NPV" ,   &fTNPV,               "NPV/F");
   tree->Branch("Mu" ,         &fTMu,                     "Mu/F");
   tree->Branch("VtxDzTruth",  &fTVtxDzTruth ,            "VtxDzTruth/F");

   // Jet vars ------------------------------------------------------------

   tree->Branch("jpt","std::vector<float>",&jpt);
   tree->Branch("jeta","std::vector<float>",&jeta);
   tree->Branch("jphi","std::vector<float>",&jphi);
   tree->Branch("jjvt","std::vector<float>",&jjvt);
   tree->Branch("jcorrjvf","std::vector<float>",&jcorrjvf);
   tree->Branch("jrpt","std::vector<float>",&jrpt);
   tree->Branch("jb2bpt","std::vector<float>",&jb2bpt);
   tree->Branch("jb2bjvt","std::vector<float>",&jb2bjvt);
   tree->Branch("jb2bjvt_70","std::vector<float>",&jb2bjvt_70);
   tree->Branch("jb2bjvt_75","std::vector<float>",&jb2bjvt_75);
   tree->Branch("jb2bjvt_85","std::vector<float>",&jb2bjvt_85);
   tree->Branch("jb2bjvt_90","std::vector<float>",&jb2bjvt_90);
   tree->Branch("jb2bcorrjvf","std::vector<float>",&jb2bcorrjvf);
   tree->Branch("jb2brpt","std::vector<float>",&jb2brpt);
   tree->Branch("jb2bpt_trk2","std::vector<float>",&jb2bpt_trk2);
   tree->Branch("jb2bjvt_trk2","std::vector<float>",&jb2bjvt_trk2);
   tree->Branch("jb2bcorrjvf_trk2","std::vector<float>",&jb2bcorrjvf_trk2);
   tree->Branch("jb2brpt_trk2","std::vector<float>",&jb2brpt_trk2);
   tree->Branch("jb2bpt_trk4","std::vector<float>",&jb2bpt_trk4);
   tree->Branch("jb2bjvt_trk4","std::vector<float>",&jb2bjvt_trk4);
   tree->Branch("jb2bcorrjvf_trk4","std::vector<float>",&jb2bcorrjvf_trk4);
   tree->Branch("jb2brpt_trk4","std::vector<float>",&jb2brpt_trk4);
   tree->Branch("tjpt","std::vector<float>",&tjpt);
   tree->Branch("jispu","std::vector<bool>",&jispu);
   tree->Branch("jishs","std::vector<bool>",&jishs);

   tree->Branch("trk2jpt","std::vector<float>", &trk2jpt);
   tree->Branch("trk2jeta","std::vector<float>",&trk2jeta);
   tree->Branch("trk2jphi","std::vector<float>",&trk2jphi);
   tree->Branch("trk2jjvt","std::vector<float>",&trk2jjvt);
   tree->Branch("trk2jcorrjvf","std::vector<float>",&trk2jcorrjvf);
   tree->Branch("trk2jrpt","std::vector<float>",&trk2jrpt);
   tree->Branch("ttrk2jpt","std::vector<float>",&ttrk2jpt);
   tree->Branch("trk2jispu","std::vector<bool>",&trk2jispu);
   tree->Branch("trk2jishs","std::vector<bool>",&trk2jishs);

   tree->Branch("trk4jpt","std::vector<float>", &trk4jpt);
   tree->Branch("trk4jeta","std::vector<float>",&trk4jeta);
   tree->Branch("trk4jphi","std::vector<float>",&trk4jphi);
   tree->Branch("trk4jjvt","std::vector<float>",&trk4jjvt);
   tree->Branch("trk4jcorrjvf","std::vector<float>",&trk4jcorrjvf);
   tree->Branch("trk4jrpt","std::vector<float>",&trk4jrpt);
   tree->Branch("ttrk4jpt","std::vector<float>",&ttrk4jpt);
   tree->Branch("trk4jispu","std::vector<bool>",&trk4jispu);
   tree->Branch("trk4jishs","std::vector<bool>",&trk4jishs);

   tree->Branch("truejetpt","std::vector<float>",&truejetpt);
   tree->Branch("truejeteta","std::vector<float>",&truejeteta);
   tree->Branch("truejetphi","std::vector<float>",&truejetphi);
   tree->Branch("truejetenergy","std::vector<float>",&truejetenergy);
   
   if(Debug()) cout <<"Analysis_treeb2bjvt::AddBranches End" << endl;
   return;
 }

void Analysis_treeb2bjvt::ResetBranches(TTree *tree){
  if(Debug()) cout <<"Analysis_treeb2bjvt::ResetBranches Start" << endl;

  // Event Info
  fTEventNumber           = -9999;
  fTRunNumber             = -9999;
  fTWeight                = -999.99;
  fTNPVtruth              = -99;
  fTNPV                   = -99;
  fTMu                    = -99;
  fTVtxDzTruth            = -999.99;

  // jet vars
  jpt     ->clear();  
  jeta     ->clear();  
  jphi     ->clear();  
  jjvt     ->clear();  
  jcorrjvf     ->clear();  
  jrpt     ->clear();  
  jb2bpt  ->clear();  
  jb2bjvt  ->clear();  
  jb2bjvt_70  ->clear();  
  jb2bjvt_75  ->clear();  
  jb2bjvt_85  ->clear();  
  jb2bjvt_90  ->clear();  
  jb2bcorrjvf  ->clear();  
  jb2brpt  ->clear();  
  jb2bpt_trk2  ->clear();  
  jb2bjvt_trk2  ->clear();  
  jb2bcorrjvf_trk2  ->clear();  
  jb2brpt_trk2  ->clear();  
  jb2bpt_trk4  ->clear();  
  jb2bjvt_trk4  ->clear();  
  jb2bcorrjvf_trk4  ->clear();  
  jb2brpt_trk4  ->clear();  
  tjpt     ->clear();  
  jispu    ->clear();  
  jishs    ->clear();  

  trk2jpt     ->clear();  
  trk2jeta     ->clear();  
  trk2jphi     ->clear();  
  trk2jjvt     ->clear();  
  trk2jcorrjvf     ->clear();  
  trk2jrpt     ->clear();  
  ttrk2jpt     ->clear();  
  trk2jispu    ->clear();  
  trk2jishs    ->clear();  

  trk4jpt     ->clear();  
  trk4jeta     ->clear();  
  trk4jphi     ->clear();  
  trk4jjvt     ->clear();  
  trk4jcorrjvf     ->clear();  
  trk4jrpt     ->clear();  
  ttrk4jpt     ->clear();  
  trk4jispu    ->clear();  
  trk4jishs    ->clear();  

  truejetpt     ->clear();  
  truejeteta     ->clear();  
  truejetphi     ->clear();  
  truejetenergy     ->clear();  

  if(Debug()) cout <<"Analysis_treeb2bjvt::ResetBranches End" << endl;
  return;
}

void Analysis_treeb2bjvt::FillEventVars(TTree *tree){
  if(Debug()) cout <<"Analysis_treeb2bjvt::FillEventVars Begin" << endl;

  // Event Info
  fTEventNumber                 = Int("EventNumber");
  fTRunNumber                   = Int("RunNumber");
  fTNPVtruth                    = Exists("NPVTruth")? Int("NPVTruth"):-1;
  fTNPV                         = Exists("NPV")? Int("NPV"):-1;
  fTWeight                      = DefaultWeight();
  fTMu                          = Float("averageIntPerXing");
  fTVtxDzTruth = vtxs()>0 ? fabs(vtx(0).x.z() - vtx(0, "Truth").x.z()): -1;

  if(Debug()) cout <<"Analysis_treeb2bjvt::FillEventVars End" << endl;
  return;
}

void Analysis_treeb2bjvt::FillJetVars(TTree* tree, int jindex,
				const MomKey jetkey){
  if(Debug()) cout <<"Analysis_treeb2bjvt::FillJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  jpt ->push_back(myjet->p.Pt());
  if(myjet->Bool("isHSJet")){
    Particle *tjet = (Particle*) myjet->Obj("AntiKt4Truth_match");
    tjpt ->push_back(tjet->p.Pt());}
  else
    tjpt ->push_back(-1);
  jispu->push_back(myjet->Bool("isPUJet"));
  jishs->push_back(myjet->Bool("isHSJet"));

  Particle *b2bjet = NULL;
  if(myjet->Objs("jetsBtoB")==1) b2bjet = (Particle*) myjet->Obj("jetsBtoB",0);
  Particle *b2bjet_trk2 = NULL;
  if(myjet->Objs("jetsBtoB_trk2")==1) b2bjet_trk2 = (Particle*) myjet->Obj("jetsBtoB_trk2",0);
  Particle *b2bjet_trk4 = NULL;
  if(myjet->Objs("jetsBtoB_trk4")==1) b2bjet_trk4 = (Particle*) myjet->Obj("jetsBtoB_trk4",0);
  
  jeta->push_back(myjet->p.Eta());
  jphi->push_back(myjet->p.Phi());
  jjvt->push_back(myjet->Float("JVT"));
  jcorrjvf->push_back(myjet->Float("corrJVF"));
  jrpt->push_back(myjet->Float("RpT"));
  if(myjet->Objs("jetsBtoB")==1) jb2bpt->push_back(b2bjet->p.Pt());
  else jb2bpt->push_back(-99);
  jb2bjvt->push_back(myjet->Float("BtoBJVT"));
  jb2bjvt_70->push_back(myjet->Float("BtoBJVT_70"));
  jb2bjvt_75->push_back(myjet->Float("BtoBJVT_75"));
  jb2bjvt_85->push_back(myjet->Float("BtoBJVT_85"));
  jb2bjvt_90->push_back(myjet->Float("BtoBJVT_90"));
  jb2bcorrjvf->push_back(myjet->Float("BtoBcorrJVF"));
  jb2brpt->push_back(myjet->Float("BtoBRpT"));
  if(myjet->Objs("jetsBtoB_trk2")==1) jb2bpt_trk2->push_back(b2bjet_trk2->p.Pt());
  else jb2bpt_trk2->push_back(-99);
  jb2bjvt_trk2->push_back(myjet->Float("BtoBJVT_trk2"));
  jb2bcorrjvf_trk2->push_back(myjet->Float("BtoBcorrJVF_trk2"));
  jb2brpt_trk2->push_back(myjet->Float("BtoBRpT_trk2"));
  if(myjet->Objs("jetsBtoB_trk4")==1) jb2bpt_trk4->push_back(b2bjet_trk4->p.Pt());
  else jb2bpt_trk4->push_back(-99);
  jb2bjvt_trk4->push_back(myjet->Float("BtoBJVT_trk4"));
  jb2bcorrjvf_trk4->push_back(myjet->Float("BtoBcorrJVF_trk4"));
  jb2brpt_trk4->push_back(myjet->Float("BtoBRpT_trk4"));

  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTrkVars End" << endl;
  return;
}

void Analysis_treeb2bjvt::FillTrk2JetVars(TTree* tree, int jindex,
				const MomKey jetkey){
  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTrk2JetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  trk2jpt ->push_back(myjet->p.Pt());
  if(myjet->Bool("isHSJet")){
    Particle *ttrk2jet = (Particle*) myjet->Obj("AntiKt2Truth_match");
    ttrk2jpt ->push_back(ttrk2jet->p.Pt());}
  else
    ttrk2jpt ->push_back(-1);
  trk2jispu->push_back(myjet->Bool("isPUJet"));
  trk2jishs->push_back(myjet->Bool("isHSJet"));
  
  trk2jeta->push_back(myjet->p.Eta());
  trk2jphi->push_back(myjet->p.Phi());
  trk2jjvt->push_back(myjet->Float("JVT"));
  trk2jcorrjvf->push_back(myjet->Float("corrJVF"));
  trk2jrpt->push_back(myjet->Float("RpT"));

  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTrk2JetVars End" << endl;
  return;
}

void Analysis_treeb2bjvt::FillTrk4JetVars(TTree* tree, int jindex,
				const MomKey jetkey){
  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTrk4JetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  trk4jpt ->push_back(myjet->p.Pt());
  if(myjet->Bool("isHSJet")){
    Particle *ttrk4jet = (Particle*) myjet->Obj("AntiKt4Truth_match");
    ttrk4jpt ->push_back(ttrk4jet->p.Pt());}
  else
    ttrk4jpt ->push_back(-1);
  trk4jispu->push_back(myjet->Bool("isPUJet"));
  trk4jishs->push_back(myjet->Bool("isHSJet"));
  
  trk4jeta->push_back(myjet->p.Eta());
  trk4jphi->push_back(myjet->p.Phi());
  trk4jjvt->push_back(myjet->Float("JVT"));
  trk4jcorrjvf->push_back(myjet->Float("corrJVF"));
  trk4jrpt->push_back(myjet->Float("RpT"));

  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTrk4JetVars End" << endl;
  return;
}

void Analysis_treeb2bjvt::FillTrueJetVars(TTree* tree, int jindex,
				    const MomKey jetkey){
  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTrueJetVars Begin" << endl;

  Particle  *myjet         = &(jet(jindex,jetkey));

  truejetpt ->push_back(myjet->p.Pt());
  truejeteta ->push_back(myjet->p.Eta());
  truejetphi ->push_back(myjet->p.Phi());
  truejetenergy ->push_back(myjet->p.E());

  if(Debug()) cout <<"Analysis_treeb2bjvt::FillTrkVars End" << endl;
  return;
}
