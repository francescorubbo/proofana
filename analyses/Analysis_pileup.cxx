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

  // jvt = new JetVertexTagger(0.2, maindir+"/utils/JetVertexTagger/data/JVTlikelihood_20140805.root");

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
  
  // vector<float> trk_pt, trk_z0_wrtPV;
  // for(int it=0; it< tracks(); ++it){
  //   trk_pt                    .push_back(track(it).p.Pt());
  //   trk_z0_wrtPV              .push_back(track(it).Float("z0_wrtPV"));
  //   track(it).Set("JVTindex", it);
  // }

  // // trk to vtx assoc                                                                                                                                               
  // vector<vector<int> > vxp_trk_index;
  // for(int iv=0; iv<vtxs(); ++iv){
  //   vector<int> assoc_track_indices;
  //   for(int it=0; it<vtx(iv).Objs("vtxTracks"); ++it){
  //     Particle* trk = (Particle*) vtx(iv).Obj("vtxTracks",it);
  //     assoc_track_indices.push_back(trk->Int("JVTindex"));
  //   }
  //   vxp_trk_index.push_back(assoc_track_indices);
  // }

  // // JVT                                                                                                                                                            
  // jvt->init(trk_pt, trk_z0_wrtPV, vxp_trk_index);

  // Particle  *myjet         = &(jet(0, "AntiKt4LCTopo"));
  // // myjet->Show();

  // vector<int> assoc_trk_indices;
  // for(int it=0; it<myjet->Int("nTrackAssoc"); ++it){
  //   Particle* trk = (Particle*) myjet->Obj("GhostAssocTrack", it);
  //   assoc_trk_indices.push_back(trk->Int("JVTindex"));
  // }
  // bool pass = (*jvt)(myjet->p.Pt(), assoc_trk_indices);

  // float corrJVF = jvt->corrJVF();
  // float RpT = jvt->RpT();
  // float JVT = jvt->JVT();

  for(int it=0; it< clusters("LCTopo"); ++it){
    Particle  *mycluster = &(cluster(it, "LCTopo"));
    associateTrackstoCluster(mycluster);
    computeJVF(mycluster);
  }

  MomKey newjets = MakeJets(fastjet::antikt_algorithm,0.4,"clustersLCTopo");
  
  
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
    ANNpoint point = annAllocPt(2,0);
    point[0] = track(it,"forjvf").p.Eta();
    point[1] = track(it,"forjvf").p.Phi();
    points[it] = point;
  }  
  for(int it=0; it< ntrks; ++it){
    ANNpoint point = annAllocPt(2,0);
    point[0] = track(it,"forjvf").p.Eta();
    point[1] = track(it,"forjvf").p.Phi()+2*PI;
    points[ntrks+it] = point;
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
  delete [] nnIdx;
  delete [] nnIdxbis;
  delete kdTree;
  annClose();
}

 void Analysis_pileup::computeJVF(Particle *thecluster)
{
  int nassoctrks = thecluster->Objs("assoctrks");
  if(nassoctrks<1){
    thecluster->Set("corrJVF",-1);
    thecluster->Set("JVF",-1);
    return;
  }
  double ptpu = 0.,ptpv = 0.;
  for(int it=0; it< nassoctrks; ++it){
    Particle* trk = (Particle*) thecluster->Obj("assoctrks",it);
    (trk->Int("origin")==0 ? ptpv : ptpu) += trk->p.Pt();
  }
  float kfac = 0.01;
  thecluster->Set("JVF",ptpv/(ptpv+ptpu));
  thecluster->Set("corrJVF",ptpv/(ptpv+ptpu/(tracks("forjvf")*kfac)));
}

 void Analysis_pileup::selectTracks()
{
  const static MomKey tracksforjvf("tracksforjvf");
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

