/**************************************************************************
 **
 **   File:         Analysis_tree.h
 **
 **   Description:  Tree filler
 **
 **
 **   Authors:      F. Rubbo
 **
 **   Created:      11/25/2014
 **
 **************************************************************************/

#ifndef Analysis_tree_h
#define Analysis_tree_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

using std::cout;
using std::endl;


class Analysis_tree : public Analysis_JetMET_Base {

 public :

  Analysis_tree(TTree* /*tree*/ = 0) {
    fDetail = false;
  }

  virtual ~Analysis_tree() { }

  ClassDef(Analysis_tree, 0);

  Bool_t  fDetail;

  virtual bool    ProcessEvent();
  virtual void    WorkerBegin();
  virtual void    WorkerTerminate();

  TTree *fTree;
  void AddBranches(TTree *tree);
  void FillTree(const MomKey JetKey);
  void FillEventVars(TTree *tree);
  void FillClVars(TTree *tree, int jindex, const MomKey JetKey);
  void FillJetVars(TTree *tree, int jindex, const MomKey JetKey,vector<float> *jetpt,vector<float> *tjetpt,vector<int> *jetncl);
  void FillTrkVars(TTree *tree, int jindex);
  void ResetBranches(TTree *tree);

  // per Event variables
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTMu;
  float fTNPVtruth;

  // 
  vector<float>   *clcenterlambda;
  vector<float>   *cljvfcorr;
  vector<float>   *cljvf;
  vector<float>   *clfem;
  vector<float>   *clpt;
  vector<float>   *cleta;
  vector<float>   *clphi;
  vector<float>   *clenergy;
  vector<float>   *trkpt;
  vector<float>   *trketa;
  vector<float>   *trkphi;
  vector<bool>    *trkispv;

  vector<int>   *j0ncl;
  vector<int>   *j1ncl;
  vector<int>   *j2ncl;
  vector<int>   *j3ncl;
  vector<int>   *j4ncl;
  vector<int>   *j5ncl;
  vector<float>   *j0pt;
  vector<float>   *j1pt;
  vector<float>   *j2pt;
  vector<float>   *j3pt;
  vector<float>   *j4pt;
  vector<float>   *j5pt;
  vector<float>   *tj0pt;
  vector<float>   *tj1pt;
  vector<float>   *tj2pt;
  vector<float>   *tj3pt;
  vector<float>   *tj4pt;
  vector<float>   *tj5pt;

  private :

};

#endif

