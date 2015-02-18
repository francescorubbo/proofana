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
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
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

  Analysis_JetMET_Base::MakeJets(fastjet::antikt_algorithm, 1.0, "truthsStable");

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
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,0.4,"clustersjvf0","jvf0");
   addTruthMatch(newjets,"AntiKt4Truth");
   AddGhostMatch("AntiKt4jvf0","tracks","clustersjvf0",fastjet::antikt_algorithm,0.4);
   for(int iJet = 0; iJet < jets("AntiKt4jvf0"); iJet++){
     Particle *myjet = &(jet(iJet, "AntiKt4jvf0"));
     vector<int> assoc_trk_indices;
     for(int it=0; it<myjet->Objs("tracksGhost"); ++it){
       Particle* trk = (Particle*) myjet->Obj("tracksGhost", it);
       assoc_trk_indices.push_back(trk->Int("JVTindex"));
     } 
     (*jvt)(myjet->p.Pt(),assoc_trk_indices);
     myjet->Set("JVT",jvt->JVT());
   }
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,1.0,"clustersjvf0","jvf0");
   addTruthMatch(newjets,"AntiKt10Truth");
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,0.4,"clustersjvf0","jvf0noarea",false);
   addTruthMatch(newjets,"AntiKt4Truth");
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,1.0,"clustersjvf0","jvf0noarea",false);
   addTruthMatch(newjets,"AntiKt10Truth");

   selectClusters(0.5,"jvf5");
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,0.4,"clustersjvf5","jvf5");
   addTruthMatch(newjets,"AntiKt4Truth");
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,1.0,"clustersjvf5","jvf5");
   addTruthMatch(newjets,"AntiKt10Truth");
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,0.4,"clustersjvf5","jvf5noarea",false);
   addTruthMatch(newjets,"AntiKt4Truth");
   newjets = MakeJetsWArea(fastjet::antikt_algorithm,1.0,"clustersjvf5","jvf5noarea",false);
   addTruthMatch(newjets,"AntiKt10Truth");

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
     // points[it][0] = track(it,"forjvf").p.Eta();
     // points[it][1] = track(it,"forjvf").p.Phi();
     points[it][0] = track(it,"forjvf").Float("eta_atCalo");
     points[it][1] = track(it,"forjvf").Float("phi_atCalo");
   }  
   for(int it=0; it< ntrks; ++it){
     // points[ntrks+it][0] = track(it,"forjvf").p.Eta();
     // points[ntrks+it][1] = track(it,"forjvf").p.Phi()+2*PI;
     points[ntrks+it][0] = track(it,"forjvf").Float("eta_atCalo");
     points[ntrks+it][1] = track(it,"forjvf").Float("phi_atCalo")+2*PI;
   } 

   ANNdist radius = 0.1*0.1;
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
    float clpt = cluster(iCl,"LCTopo").p.Pt();
    // if(clpt<2.)
      if(cljvf>-1 && cljvf < jvfcut) continue;
    Add(clustersjvf,&cluster(iCl,"LCTopo"));
  }
}

void Analysis_pileup::addTruthMatch(const MomKey JetType, const MomKey TruthJetType){

  vector<int> matched;
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
      if(std::find(matched.begin(), matched.end(), iTrueJ) != matched.end())
	continue;
      Particle* trueJ = &(jet(iTrueJ, TruthJetType));

      float dR = (thejet->p).DeltaR(trueJ->p);
      if(dR<mindR){ mindR = dR; minDRindex = iTrueJ;}
      if(dR < Float("MAXJETTRUTHMATCHDR") && maxPt < trueJ->p.Pt())   
	{ maxPt = trueJ->p.Pt(); maxPtIndex = iTrueJ;} 
    }//true jets

    if(maxPtIndex != -1){
      matched.push_back(maxPtIndex);
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

MomKey Analysis_pileup::MakeJetsWArea(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType, const MomKey extra,bool doareasub){

  const static MomKey SJetKey("jets");
  vector<fastjet::PseudoJet> inputConst = ObjsToPJ(constType);

  fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);
  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::JetDefinition jetDefkt(fastjet::kt_algorithm,0.4);
  fastjet::AreaDefinition active_area = fastjet::AreaDefinition(fastjet::active_area);
  fastjet::AreaDefinition voronoi_area = fastjet::AreaDefinition(fastjet::voronoi_area,fastjet::VoronoiAreaSpec(0.9));
  fastjet::ClusterSequenceArea clustSeq(inputConst, jetDef, active_area);

  fastjet::JetMedianBackgroundEstimator bge(jselector,jetDefkt,voronoi_area);
  fastjet::Subtractor subtractor(&bge);
  bge.set_particles(inputConst);
  double rho = bge.rho();
  MomKey rhokey("rho_");
  rhokey+=extra;
  Set(rhokey,rho);

  vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(10.));
  vector<fastjet::PseudoJet> subtractedJets = inclusiveJets;
  if(doareasub)
    subtractedJets = subtractor(inclusiveJets);

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

  for(unsigned int iJet = 0 ; iJet < subtractedJets.size() ; iJet++){
  	fastjet::PseudoJet jet = subtractedJets[iJet];
	if(jet.perp2()<5.*5.) continue;
  	vector<fastjet::PseudoJet> constituents = jet.constituents();
  	Particle* jetP = new Particle();
  	jetP->p.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.e());
	jetP->Set("area", jet.area());
  	static const MomKey ConsKey("constituents");
  	jetP->AddVec(ConsKey);
	double numwidth = 0.;
	double denwidth = 0.;
  	for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
  		const PJ_Info* info = &(constituents[iCons].user_info<PJ_Info>());
  		jetP->Add(ConsKey, info->Pointer);
		double ptcons = constituents[iCons].pt();
		double drcons = constituents[iCons].delta_R(jet);
		numwidth += drcons*ptcons;
		denwidth += ptcons;
  	} // end loop over cons
	jetP->Set("width",numwidth/denwidth);
	jetP->Set("constscale_pt",jet.pt());
  	Add(FFinalKey, jetP);
  }// end loop over jets

  return FinalKey;
}
