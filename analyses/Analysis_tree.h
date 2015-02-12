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
  void FillJetVars(TTree *tree, int jindex, const MomKey JetKey,
		   vector<float> *jetpt,vector<float> *jeteta,
		   vector<float> *jetphi,vector<float> *tjetpt,vector<float> *tjeteta,
		   vector<int> *jetncl,vector<float> *jetm,vector<float> *jetw,MomKey truthjetkey = "AntiKt4Truth_match");
  void FillTrueJetVars(TTree *tree, int jindex, const MomKey JetKey);
  void FillTrueLargeJetVars(TTree *tree, int jindex, const MomKey JetKey);
  void FillTrkVars(TTree *tree, int jindex);
  void ResetBranches(TTree *tree);

  // per Event variables
  int   fTEventNumber;
  int   fTRunNumber;
  float fTWeight;
  float fTMu;
  float fTRho;
  float fTNPVtruth;
  float fTNPV;

  // 
  vector<float>   *clcenterlambda;
  vector<float>   *cljvfcorr;
  vector<float>   *clfem;
  vector<float>   *clpt;
  vector<float>   *cleta;
  vector<float>   *clphi;
  vector<float>   *clenergy;
  vector<float>   *trkpt;
  vector<float>   *trketa;
  vector<float>   *trkphi;
  vector<bool>    *trkispv;

  vector<int>     *j0ncl;
  vector<int>     *j5ncl;
  vector<float>     *j0mass;
  vector<float>     *j5mass;
  vector<float>     *j0width;
  vector<float>     *j5width;
  vector<float>   *j0pt;
  vector<float>   *j0eta;
  vector<float>   *j0phi;
  vector<float>   *j0jvt;
  vector<float>   *j5pt;
  vector<float>   *j5eta;
  vector<float>   *j5phi;
  vector<float>   *tj0pt;
  vector<float>   *tj5pt;
  vector<float>   *tj0eta;
  vector<float>   *tj5eta;

  vector<int>     *lj0ncl;
  vector<int>     *lj5ncl;
  vector<float>     *lj0mass;
  vector<float>     *lj5mass;
  vector<float>     *lj0width;
  vector<float>     *lj5width;
  vector<float>   *lj0pt;
  vector<float>   *lj0eta;
  vector<float>   *lj0phi;
  vector<float>   *lj5pt;
  vector<float>   *lj5eta;
  vector<float>   *lj5phi;
  vector<float>   *tlj0pt;
  vector<float>   *tlj5pt;
  vector<float>   *tlj0eta;
  vector<float>   *tlj5eta;

  vector<int>     *jnoarea0ncl;
  vector<int>     *jnoarea5ncl;
  vector<float>     *jnoarea0mass;
  vector<float>     *jnoarea5mass;
  vector<float>     *jnoarea0width;
  vector<float>     *jnoarea5width;
  vector<float>   *jnoarea0pt;
  vector<float>   *jnoarea0eta;
  vector<float>   *jnoarea0phi;
  vector<float>   *jnoarea5pt;
  vector<float>   *jnoarea5eta;
  vector<float>   *jnoarea5phi;
  vector<float>   *tjnoarea0pt;
  vector<float>   *tjnoarea5pt;
  vector<float>   *tjnoarea0eta;
  vector<float>   *tjnoarea5eta;

  vector<int>     *ljnoarea0ncl;
  vector<int>     *ljnoarea5ncl;
  vector<float>     *ljnoarea0mass;
  vector<float>     *ljnoarea5mass;
  vector<float>     *ljnoarea0width;
  vector<float>     *ljnoarea5width;
  vector<float>   *ljnoarea0pt;
  vector<float>   *ljnoarea0eta;
  vector<float>   *ljnoarea0phi;
  vector<float>   *ljnoarea5pt;
  vector<float>   *ljnoarea5eta;
  vector<float>   *ljnoarea5phi;
  vector<float>   *tljnoarea0pt;
  vector<float>   *tljnoarea5pt;
  vector<float>   *tljnoarea0eta;
  vector<float>   *tljnoarea5eta;

  vector<float>   *truejetpt;
  vector<float>   *truejeteta;
  vector<float>   *truejetphi;
  vector<float>   *truejetenergy;

  vector<float>   *truelargejetpt;
  vector<float>   *truelargejeteta;
  vector<float>   *truelargejetphi;
  vector<float>   *truelargejetenergy;

  private :

};

#endif

