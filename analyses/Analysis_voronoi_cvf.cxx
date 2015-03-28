/**************************************************************************
 **
 **   File:         Analysis_voronoi.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_voronoi_cvf_cxx

#include "Analysis_voronoi_cvf.h"
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
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

const double PI  =3.141592653589793238463;
///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_voronoi_cvf::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_voronoi_cvf: DEBUG In WorkerBegin() " << endl;

  Analysis_pileup::WorkerBegin();

  std::string maindir(gSystem->Getenv("PROOFANADIR"));
  if(maindir==".") maindir+="/libProofAna";
  cout << maindir << " is the maindir we work from! " << endl;
  cout << gSystem->Exec(TString("ls "+maindir).Data()) << endl;

  if (Debug()) cout << "Analysis_voronoi_cvf: DEBUG Finish WorkerBegin()" << endl;  
} 

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_voronoi_cvf::ProcessEvent()
{

  if (Debug()) cout << "Analysis_voronoi_cvf: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();

  Analysis_JetMET_Base::AddStableParticles();
  Analysis_JetMET_Base::MakeJets(fastjet::antikt_algorithm, 1.0, "truthsStable");

SetCVF("LCTopo");
selectLCTopoCVFClusters();
 // cout << "CVF LCTopo jets: " << clusters("LCTopoZeroSigma_CVFcutx") << endl; 

  MakeVoronoiClusters(fastjet::kt_algorithm , 0.4, "clustersLCTopo");
SetCVF("Voronoi");

  selectClusters("Voronoi",0);
  selectClusters("Voronoi",1);
  selectClusters("Voronoi",10);

  selectClusters("Voronoi",0,true,5);
  selectClusters("Voronoi",0,true,-1);

  selectClusters("Voronoi",1,true,5);
  selectClusters("Voronoi",1,true,-1);


selectVoronoiCVFClusters(5);
selectVoronoiCVFClusters(-1);
float spreadr=0.4;
MakeSpreadVoronoiClusters(spreadr,"Voronoi");
MakeSpreadVoronoiClusters(spreadr,"Voronoi_CVFcut5");
MakeSpreadVoronoiClusters(spreadr,"Voronoi_CVFcutx");
//MakeSpreadVoronoiClusters(spreadr,"Voronoi_CVFcut");

//Show();


  MomKey newjets;
  
  //Standard cuts:
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersLCTopo_CVFcut","j0_CVF",true);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm,0.4,"clustersLCTopo","j0",true);
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm,0.4,"clustersLCTopo","jnoarea0",false);
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");

  //Voronoi nsigma:
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiZeroSigma","VoronoiZeroSigma",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiOneSigma","VoronoiOneSigma",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiTenSigma","VoronoiTenSigma",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");

  //Voronoi with CVF:
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiZeroSigma_CVFcut5","VoronoiZeroSigma_CVF5",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiZeroSigma_CVFcutx","VoronoiZeroSigma_CVFx",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiOneSigma_CVFcut5","VoronoiOneSigma_CVF5",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiOneSigma_CVFcutx","VoronoiOneSigma_CVFx",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");

  //Spreading with CVF:
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoi_SpreadPT","Voronoi_SpreadPT",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoi_CVFcut5_SpreadPT","Voronoi_CVF5_SpreadPT",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoi_CVFcutx_SpreadPT","Voronoi_CVFx_SpreadPT",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");

 // Show();
 
  return true;
}

/// WorkerTerminate: clean up
///=========================================
void Analysis_voronoi_cvf::WorkerTerminate()
{

}

void Analysis_voronoi_cvf::SetCVF(MomKey key){
//assign indices to tracks
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

 for(int it=0; it< clusters(key); ++it){
     Particle  *mycluster = &(cluster(it, key));
     associateTrackstoCluster(mycluster);
     vector<int> assoc_trk_indices;
     for(int it2=0; it2<mycluster->Objs("assoctrks"); ++it2){
       Particle* trk = (Particle*) mycluster->Obj("assoctrks", it2);
       assoc_trk_indices.push_back(trk->Int("JVTindex"));
     }
     (*jvt)(mycluster->p.Pt(),assoc_trk_indices);
     mycluster->Set("CVF",jvt->corrJVF());
   }
}

void Analysis_voronoi_cvf::MakeSpreadVoronoiClusters(float spreadr,MomKey key){

vector<float> spreadPT(clusters(key));
for(int iCl = 0; iCl < clusters(key); iCl++){
	spreadPT[iCl]=cluster(iCl,key).Float("correctedPT");
}

  ANNpointArray points = annAllocPts(2*clusters(key),2);
  for(int i=0; i<clusters(key); i++){
  //set up points to look through: phi in [-pi,3pi]
points[i][0] = cluster(i,key).p.Eta();
points[i][1] = cluster(i,key).p.Phi();
points[clusters(key)+i][0] = cluster(i,key).p.Eta();
points[clusters(key)+i][1] = cluster(i,key).p.Phi()+2*PI;
} 

for(int i = 0; i < clusters(key); i++){
	if(!(spreadPT[i]<0)) continue;
	//find closest positive PT cluster:
/*	int iclosest=0;
	float drmin=100;
if(i==3)	cout << "dr: " << endl;
	for(int j = 0; j < clusters(key); j++){
		
		if(!(spreadPT[j]>0)) continue;
		float dr=cluster(i,key).p.DeltaR(cluster(j,key).p);
		if(dr<drmin) {drmin=dr; iclosest=j;}
	}*/
  ANNpoint qpoint = annAllocPt(2,0);
  qpoint[0] = cluster(i,key).p.Eta();
  qpoint[1] = cluster(i,key).p.Phi();
  qpoint[1]+=2*PI*(qpoint[1]<0); //point you're looking at has phi in [0,2pi]
  
ANNdist radius = spreadr*spreadr;
ANNkd_tree* kdTree = new ANNkd_tree(points,2*clusters(key),2);
ANNidxArray nnIdx = new ANNidx[2*clusters(key)];
ANNdistArray dists = new ANNdist[2*clusters(key)];
int nclfound = kdTree->annkFRSearch(qpoint,radius,2*clusters(key),nnIdx,dists,0.0);

for(int j=0; j<nclfound; j++){
int realid=nnIdx[j]%clusters(key);
if(realid==i) continue; //don't spread to itself!
if(!(spreadPT[realid]>0)) continue; //only spread to positive PT cells
//cout << dists[j] << endl; //algorithm lists neighbors in order of distance; dist=square(dr)
//cout << "Before: " << realid << ";" << dists[j] << ";" << spreadPT[i] << ";" << spreadPT[realid] << endl;
	if(fabs(spreadPT[i])>spreadPT[realid]) { spreadPT[i]+=spreadPT[realid];spreadPT[realid]=0;}
	else { spreadPT[realid]+=spreadPT[i];spreadPT[i]=0;}
//cout << "After: " << realid << ";" << dists[j] << ";" << spreadPT[i] << ";" << spreadPT[realid] << endl;
}
//cout << i << ";" << cluster(i,key).Float("correctedPT") << ";" << spreadPT[i]<< endl;
  annDeallocPt(qpoint);
  delete [] nnIdx;
  delete [] dists;
  delete kdTree;
}
  annDeallocPts(points);
  annClose();
/*float totalcorrpt=0, totalspreadpt=0;
for(int i=0; i<clusters(key); i++){ totalcorrpt+=cluster(i,key).Float("correctedPT"); totalspreadpt+=spreadPT[i];}
cout << totalcorrpt << ";" << totalspreadpt << endl; //should be the same*/


MomKey newkey="clusters"+key;
//newkey+="clustersVoronoi";
//newkey+=(int)10*spreadr;
newkey+="_SpreadPT";
//const MomKey clusterskey(newkey);
AddVec(newkey);

for(int iCl = 0; iCl < clusters(key); iCl++){
if(!(spreadPT[iCl]>0)) continue;
  		Particle * constituentP = new Particle();
		constituentP->p.SetPtEtaPhiM(spreadPT[iCl],
					     cluster(iCl,key).p.Eta(),
					     cluster(iCl,key).p.Phi(),
					     cluster(iCl,key).p.M());
  		Add(newkey,constituentP);
}

}


void Analysis_voronoi_cvf::selectClusters(MomKey key,int nsigma,bool doCVF,float threshold)
{
MomKey clusterskey="clusters"+key;
//key+="clustersVoronoi";
//key+=clusters;
if(nsigma==0)  clusterskey+="ZeroSigma";
if(nsigma==1)  clusterskey+="OneSigma";
if(nsigma==2)  clusterskey+="TwoSigma";
if(nsigma==3)  clusterskey+="ThreeSigma";
if(nsigma==4)  clusterskey+="FourSigma";
if(nsigma==10)  clusterskey+="TenSigma";

if(doCVF) clusterskey+="_CVFcut";
if(doCVF){
if(threshold>0) clusterskey+=TString::Format("%i",(int)threshold);
else clusterskey+="x";
}
//const MomKey clusterskey(key);
  AddVec(clusterskey);
  float sigma=Float("sigma_voronoi");
  for(int iCl = 0; iCl < clusters(key); iCl++){
    float area = cluster(iCl,key).Float("area");
    float clpt = cluster(iCl,key).Float("correctedPT");
    float cvf = cluster(iCl,key).Float("CVF");
    //if no CVF include only clusters greater than sigma*nsigma
    if(!doCVF) if(clpt<sqrt(area)*sigma*(float)nsigma) continue;
    //if CVF:
    //for CVF=0, only include if pt>threshold. If threshold>100 don't include any.
    //for CVF=1, include all > 0
    //for CVF=-1, include only greater than sigma*nsigma
    if(doCVF){
	if(fabs(cvf)<0.5){
		if(threshold>0){ if(clpt<threshold) continue;}
		else continue;
	}
	if(cvf>0.5) if(clpt<0) continue;
	if(cvf<-0.5) if(clpt<sqrt(area)*sigma*(float)nsigma) continue;
    }
  		Particle * constituentP = new Particle();
		constituentP->p.SetPtEtaPhiM(clpt,
					     cluster(iCl,key).p.Eta(),
					     cluster(iCl,key).p.Phi(),
					     cluster(iCl,key).p.M());
    Add(clusterskey,constituentP);
  }
  
}


void Analysis_voronoi_cvf::selectLCTopoCVFClusters()
{
MomKey key="LCTopo";
MomKey clusterskey="clusters"+key;
clusterskey+="_CVFcut";
  AddVec(clusterskey);
  for(int iCl = 0; iCl < clusters(key); iCl++){
    float cvf = cluster(iCl,key).Float("CVF");
        if(fabs(cvf)<0.5) continue;
    
  		Particle * constituentP = new Particle();
		constituentP->p.SetPtEtaPhiM(cluster(iCl,key).p.Pt(),
					     cluster(iCl,key).p.Eta(),
					     cluster(iCl,key).p.Phi(),
					     cluster(iCl,key).p.M());
		constituentP->Set("CVF",cvf);
    Add(clusterskey,constituentP);
  }
 
}

void Analysis_voronoi_cvf::selectVoronoiCVFClusters(float threshold){
MomKey key="Voronoi";
MomKey clusterskey="clusters"+key;
clusterskey+="_CVFcut";
if(threshold>0) clusterskey+=TString::Format("%i",(int)threshold);
else clusterskey+="x";
  AddVec(clusterskey);
  for(int iCl = 0; iCl < clusters(key); iCl++){
    float area = cluster(iCl,key).Float("area"); 
    float clpt = cluster(iCl,key).Float("correctedPT");
    float cvf = cluster(iCl,key).Float("CVF");
        if(fabs(cvf)<0.5){
		if(threshold>0){if(clpt<threshold) continue;}
		else continue;
	}
    
  		Particle * constituentP = new Particle();
		constituentP->p.SetPtEtaPhiM(cluster(iCl,key).p.Pt(),
					     cluster(iCl,key).p.Eta(),
					     cluster(iCl,key).p.Phi(),
					     cluster(iCl,key).p.M());
		constituentP->Set("correctedPT",clpt);
		constituentP->Set("area",area);
		constituentP->Set("CVF",cvf);
    Add(clusterskey,constituentP);
  }
 
}

MomKey Analysis_voronoi_cvf::MakeVoronoiClusters(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType){

  const static MomKey SJetKey("clustersVoronoi");
  vector<fastjet::PseudoJet> inputConst = ObjsToPJ(constType); //constType="clustersLCTopo"
  fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);

  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::AreaDefinition area_def(fastjet::voronoi_area, fastjet::VoronoiAreaSpec(0.9));

  fastjet::ClusterSequenceArea clustSeq(inputConst, jetDef, area_def);  
  fastjet::JetMedianBackgroundEstimator bge(jselector,jetDef,area_def);
  bge.set_particles(inputConst);

  vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(0));

  MomKey FinalKey = SJetKey; //"clustersVoronoi"

  AddVec(SJetKey); //AddVec "clustersVoronoi"

//  float rho = Float("Eventshape_rhoKt4LC")*0.001;
  Set("rho_voronoi",bge.rho());
  Set("sigma_voronoi",bge.sigma());

  for(unsigned int iJet = 0 ; iJet < inclusiveJets.size() ; iJet++){
  	fastjet::PseudoJet jet = inclusiveJets[iJet];
  	vector<fastjet::PseudoJet> constituents = jet.constituents();
  	static const MomKey ConsKey("constituents");
  	for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
  		Particle * constituentP = new Particle();
		float area = constituents[iCons].area();
		float correctedPt = constituents[iCons].pt()-Float("rho_voronoi")*area;
//		if(correctedPt<0.) continue;
		constituentP->p.SetPtEtaPhiM(constituents[iCons].pt(),
					     constituents[iCons].eta(),
					     constituents[iCons].phi(),
					     constituents[iCons].m());
		constituentP->Set("area",area);
		constituentP->Set("correctedPT",correctedPt);
		
  		Add(FinalKey, constituentP);
  	} // end loop over cons
  }// end loop over jets

  return FinalKey;//"clustersVoronoi"
}
