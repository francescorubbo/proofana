/**************************************************************************
 **
 **   File:         Analysis_treeb2bjvt.h
 **
 **   Description:  Tree filler
 **
 **
 **   Authors:      F. Rubbo
 **
 **   Created:      11/25/2014
 **
 **************************************************************************/

#ifndef Analysis_treeb2bjvt_h
#define Analysis_treeb2bjvt_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

using std::cout;
using std::endl;


class Analysis_treeb2bjvt : public Analysis_JetMET_Base {

 public :

  Analysis_treeb2bjvt(TTree* /*tree*/ = 0) {
    fDetail = false;
  }

  virtual ~Analysis_treeb2bjvt() { }

  ClassDef(Analysis_treeb2bjvt, 0);

  Bool_t  fDetail;

  virtual bool    ProcessEvent();
  virtual void    WorkerBegin();
  virtual void    WorkerTerminate();

  TTree *fTree;
  void AddBranches(TTree *tree);
  void FillTree(const MomKey JetKey);
  void FillEventVars(TTree *tree);
  void FillClVars(TTree *tree, int jindex, const MomKey JetKey);
  void FillJetVars(TTree *tree, int jindex, const MomKey JetKey);
  void FillTrueJetVars(TTree *tree, int jindex, const MomKey JetKey);
  void FillTrkVars(TTree *tree, int jindex);
  void ResetBranches(TTree *tree);

  // per Event variables
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTMu;
  float fTNPVtruth;

  // 
  vector<float>   *jpt;
  vector<float>   *jeta;
  vector<float>   *jphi;
  vector<float>   *jb2bjvt;
  vector<float>   *tjpt;

  vector<float>   *truejetpt;
  vector<float>   *truejeteta;
  vector<float>   *truejetphi;
  vector<float>   *truejetenergy;

  private :

};

#endif

