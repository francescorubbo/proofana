/**************************************************************************
 **
 **   File:         Analysis_voronoi.cxx
 **
 **   Description:  See header
 **                 
 **   Authors:      F. Rubbo
 **
 **************************************************************************/

#define Analysis_voronoi_cxx

#include "Analysis_voronoi.h"
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

///=========================================
/// WorkerBegin: setup binning, etc
///=========================================
 void Analysis_voronoi::WorkerBegin()
 {
  if (Debug()) cout << "Analysis_voronoi: DEBUG In WorkerBegin() " << endl;

  Analysis_JetMET_Base::WorkerBegin();

  std::string maindir(gSystem->Getenv("PROOFANADIR"));
  if(maindir==".") maindir+="/libProofAna";
  cout << maindir << " is the maindir we work from! " << endl;
  cout << gSystem->Exec(TString("ls "+maindir).Data()) << endl;

  if (Debug()) cout << "Analysis_voronoi: DEBUG Finish WorkerBegin()" << endl;  
} 

///=========================================
/// ProcessEvent: run the analysis
///=========================================  
bool Analysis_voronoi::ProcessEvent()
{

  if (Debug()) cout << "Analysis_voronoi: DEBUG In ProcessEvent(): RunNumber = " << RunNumber() 
		    << "; Channel = " << ChannelNumber() << "; EventNumber = " << EventNumber() << endl;

  OutputDir()->cd();

  MakeVoronoiClusters(fastjet::kt_algorithm , 0.4, "clustersLCTopo");
  selectClusters();
  MomKey newjets;
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 0.4, "clustersVoronoiOneSigma","VoronoiOneSigma",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt4Truth");
  newjets = Analysis_pileup::MakeJetsWArea(fastjet::antikt_algorithm, 1.0, "clustersVoronoiOneSigma","VoronoiOneSigma",false);  
  Analysis_pileup::addTruthMatch(newjets,"AntiKt10Truth");
  
  return true;
}

/// WorkerTerminate: clean up
///=========================================
void Analysis_voronoi::WorkerTerminate()
{

}

void Analysis_voronoi::selectClusters()
{
  const MomKey clusterskey("clustersVoronoiOneSigma");
  AddVec(clusterskey);
  float sigma = Float("sigma_voronoi");
  for(int iCl = 0; iCl < clusters("Voronoi"); iCl++){
    float area = cluster(iCl,"Voronoi").Float("area");
    float clpt = cluster(iCl,"Voronoi").p.Pt();
    if(clpt<area*sqrt(area)*sigma) continue;
    Add(clusterskey,&cluster(iCl,"Voronoi"));
  }
}

MomKey Analysis_voronoi::MakeVoronoiClusters(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType){

  const static MomKey SJetKey("clustersVoronoi");
  vector<fastjet::PseudoJet> inputConst = ObjsToPJ(constType);
  fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);

  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::AreaDefinition area_def(fastjet::voronoi_area, fastjet::VoronoiAreaSpec(0.9));

  fastjet::ClusterSequenceArea clustSeq(inputConst, jetDef, area_def);  
  fastjet::JetMedianBackgroundEstimator bge(jselector,jetDef,area_def);
  bge.set_particles(inputConst);

  vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(0));

  MomKey FinalKey = SJetKey;

  AddVec(SJetKey);

  float rho = Float("Eventshape_rhoKt4LC")*0.001;
  Set("rho_voronoi",bge.rho());
  Set("sigma_voronoi",bge.sigma());

  for(unsigned int iJet = 0 ; iJet < inclusiveJets.size() ; iJet++){
  	fastjet::PseudoJet jet = inclusiveJets[iJet];
  	vector<fastjet::PseudoJet> constituents = jet.constituents();
  	static const MomKey ConsKey("constituents");
  	for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
  		Particle * constituentP = new Particle();
		float area = constituents[iCons].area();
		float correctedPt = constituents[iCons].pt()-rho*area;
		if(correctedPt<0.) continue;
		constituentP->p.SetPtEtaPhiM(correctedPt,
					     constituents[iCons].eta(),
					     constituents[iCons].phi(),
					     constituents[iCons].m());
		constituentP->Set("area",area);
		
  		Add(FinalKey, constituentP);
  	} // end loop over cons
  }// end loop over jets

  return FinalKey;
}
