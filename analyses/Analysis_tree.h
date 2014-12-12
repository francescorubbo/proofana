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
#include "JetVertexTagger/JetVertexTagger.h"

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
  void FillTrkVars(TTree *tree, int jindex);
  void ResetBranches(TTree *tree);

  // per Event variables (but filled on jet-by-jet basis)
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTMu;
  float fTNPVtruth;

  // jet by jet
  vector<float>   *clcenterlambda;
  vector<float>   *clfem;
  vector<float>   *clpt;
  vector<float>   *cleta;
  vector<float>   *clphi;
  vector<float>   *clenergy;
  vector<float>   *trkpt;
  vector<float>   *trketa;
  vector<float>   *trkphi;
  vector<bool>   *trkispv;

  private :

};

#endif

