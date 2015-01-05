/**************************************************************************
 **
 **   File:         Analysis_pileup.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_pileup_cxx

#include "Analysis_pileup.h"
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
 void Analysis_pileup::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_pileup: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  std::string maindir(gSystem->Getenv("PROOFANADIR"));
  if(maindir==".") maindir+="/libProofAna";
  cout << maindir << " is the maindir we work from! " << endl;
  cout << gSystem->Exec(TString("ls "+maindir).Data()) << endl;

  jvt = new JetVertexTagger(0.2, maindir+"/utils/JetVertexTagger/data/JVTlikelihood_20140805.root");

  if (Debug()) cout << "Analysis_pileup: DEBUG Finish WorkerBegin()" << endl;  
} 

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_pileup::ProcessEvent()
{

  if (Debug()) cout << "Analysis_pileup: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
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
  
   for(int it=0; it< clusters("LCTopo"); ++it){
     Particle  *mycluster = &(cluster(it, "LCTopo"));
     associateTrackstoCluster(mycluster);
     vector<int> assoc_trk_indices;
     for(int it=0; it<mycluster->Objs("assoctrks"); ++it){
       Particle* trk = (Particle*) mycluster->Obj("assoctrks", it);
       assoc_trk_indices.push_back(trk->Int("JVTindex"));
     }
     (*jvt)(mycluster->p.Pt(),assoc_trk_indices);
     mycluster->Set("corrJVF",jvt->corrJVF());
   }

   MomKey newjets;
   selectClusters(0.0,"jvf0");
   newjets = MakeJets(fastjet::antikt_algorithm,0.4,"clustersjvf0","jvf0");
   addTruthMatch(newjets,"AntiKt4Truth");
   selectClusters(0.1,"jvf1");
   newjets = MakeJets(fastjet::antikt_algorithm,0.4,"clustersjvf1","jvf1");
   addTruthMatch(newjets,"AntiKt4Truth");
   selectClusters(0.2,"jvf2");
   newjets = MakeJets(fastjet::antikt_algorithm,0.4,"clustersjvf2","jvf2");
   addTruthMatch(newjets,"AntiKt4Truth");
   selectClusters(0.3,"jvf3");
   newjets = MakeJets(fastjet::antikt_algorithm,0.4,"clustersjvf3","jvf3");
   addTruthMatch(newjets,"AntiKt4Truth");
   selectClusters(0.4,"jvf4");
   newjets = MakeJets(fastjet::antikt_algorithm,0.4,"clustersjvf4","jvf4");
   addTruthMatch(newjets,"AntiKt4Truth");
   selectClusters(0.5,"jvf5");
   newjets = MakeJets(fastjet::antikt_algorithm,0.4,"clustersjvf5","jvf5");
   addTruthMatch(newjets,"AntiKt4Truth");
   
   return true;
}

/// WorkerTerminate: clean up
///=========================================
void Analysis_pileup::WorkerTerminate()
{


  // Nothing more

}

void Analysis_pileup::associateTrackstoCluster(Particle *thecluster)
{
  ANNpoint qpoint = annAllocPt(2,0);
  qpoint[0] = thecluster->p.Eta();
  qpoint[1] = thecluster->p.Phi();
  ANNpoint qpointbis = annAllocPt(2,0);
  qpointbis[0] = thecluster->p.Eta();
  qpointbis[1] = thecluster->p.Phi() + 2*PI;
  
  selectTracks();
  int ntrks = tracks("forjvf");
  ANNpointArray points = annAllocPts(2*ntrks,2);
   for(int it=0; it< ntrks; ++it){
     points[it][0] = track(it,"forjvf").p.Eta();
     points[it][1] = track(it,"forjvf").p.Phi();
   }  
   for(int it=0; it< ntrks; ++it){
     points[ntrks+it][0] = track(it,"forjvf").p.Eta();
     points[ntrks+it][1] = track(it,"forjvf").p.Phi()+2*PI;
   } 

   ANNdist radius = 0.3*0.3;
   ANNkd_tree* kdTree = new ANNkd_tree(points,2*ntrks,2);
   ANNidxArray nnIdx = new ANNidx[2*ntrks];
   ANNidxArray nnIdxbis = new ANNidx[2*ntrks];
   int ntrksfound = kdTree->annkFRSearch(qpoint,radius,2*ntrks,nnIdx,NULL,0.0);
   int ntrksfoundbis = kdTree->annkFRSearch(qpointbis,radius,2*ntrks,nnIdxbis,NULL,0.0);

   thecluster->AddVec("assoctrks");
   for(int it=0;it<ntrksfound;++it){
     int correctIdx = nnIdx[it];
     if(points[correctIdx][1]>2*PI || points[correctIdx][1]<0.) continue;
     if(correctIdx>ntrks-1) correctIdx -= ntrks;
     thecluster->Add("assoctrks",&(track(correctIdx,"forjvf")));
   }
   for(int it=0;it<ntrksfoundbis;++it){
     int correctIdx = nnIdxbis[it];
     if(points[correctIdx][1]>2*PI || points[correctIdx][1]<0.) continue;
     if(correctIdx>ntrks-1) correctIdx -= ntrks;
     thecluster->Add("assoctrks",&(track(correctIdx,"forjvf")));
   }
  annDeallocPt(qpoint);
  annDeallocPt(qpointbis);
  annDeallocPts(points);
  delete [] nnIdx;
  delete [] nnIdxbis;
  delete kdTree;
  annClose();
}

 void Analysis_pileup::selectTracks()
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

void Analysis_pileup::selectClusters(float jvfcut,string suffix)
{
  const MomKey clustersjvf("clusters"+suffix);
  AddVec(clustersjvf);
  for(int iCl = 0; iCl < clusters("LCTopo"); iCl++){
    float cljvf = cluster(iCl,"LCTopo").Float("corrJVF");
    if(cljvf>-1 && cljvf < jvfcut) continue;
    Add(clustersjvf,&cluster(iCl,"LCTopo"));
  }
}

void Analysis_pileup::addTruthMatch(const MomKey JetType, const MomKey TruthJetType){
  for(int iJet = 0; iJet < jets(JetType); iJet++) {
    Particle* thejet = &(jet(iJet, JetType));
    
    // defaults
    thejet->Set("isPUJet", false);
    thejet->Set("isHSJet", false);
    
    float mindR= 999.99; 
    float maxPt=-999.99; 
    int minDRindex =-1;
    int maxPtIndex =-1;
    for(int iTrueJ=0; iTrueJ < jets(TruthJetType); ++iTrueJ){
      Particle* trueJ = &(jet(iTrueJ, TruthJetType));

      float dR = (thejet->p).DeltaR(trueJ->p);
      if(dR<mindR){ mindR = dR; minDRindex = iTrueJ;}
      if(dR < Float("MAXJETTRUTHMATCHDR") && maxPt < trueJ->p.Pt())   
	{ maxPt = trueJ->p.Pt(); maxPtIndex = iTrueJ;} 
    }//true jets

    if(maxPtIndex != -1){
      thejet->Add(TruthJetType+"_match", &(jet(maxPtIndex, TruthJetType)));
      thejet->Set("isHSJet", true);  
      thejet->Set("isPUJet", false); 
      if (!(jet(maxPtIndex, TruthJetType).Exists(JetType+"_match"))) 
	jet(maxPtIndex, TruthJetType).Add(JetType+"_match", thejet, true); 
    } else if (mindR> (Float("MINJETPUDR"))){
      thejet->Set("isPUJet", true);
      thejet->Set("isHSJet", false);  
    } else {
      thejet->Set("isPUJet", false);
      thejet->Set("isHSJet", false);  
    }//alternatives to maxPtIndex != -1
  }//jet loop
}
