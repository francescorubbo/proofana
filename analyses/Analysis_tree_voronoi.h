/**************************************************************************
 **
 **   File:         Analysis_tree_voronoi.h
 **
 **   Description:  Tree filler
 **
 **
 **   Authors:      F. Rubbo
 **
 **   Created:      11/25/2014
 **
 **************************************************************************/

#ifndef Analysis_tree_voronoi_h
#define Analysis_tree_voronoi_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

using std::cout;
using std::endl;


class Analysis_tree_voronoi : public Analysis_JetMET_Base {

 public :

  Analysis_tree_voronoi(TTree* /*tree*/ = 0) {
    fDetail = false;
  }

  virtual ~Analysis_tree_voronoi() { }

  ClassDef(Analysis_tree_voronoi, 0);

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
  float fTSigma;
  float fTNPVtruth;
  float fTNPV;

  // 
  vector<float>   *clcenterlambda;
  vector<float>   *cljvfcorr;
  vector<float>   *clfem;
  vector<float>   *clpt;
  vector<float>   *clpt_orig;
  vector<float>   *clarea;
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

  vector<int>     *jvoro0ncl;
  vector<float>   *jvoro0mass;
  vector<float>   *jvoro0width;
  vector<float>   *jvoro0pt;
  vector<float>   *jvoro0eta;
  vector<float>   *jvoro0phi;
  vector<float>   *tjvoro0pt;
  vector<float>   *tjvoro0eta;

  vector<int>     *ljvoro0ncl;
  vector<float>   *ljvoro0mass;
  vector<float>   *ljvoro0width;
  vector<float>   *ljvoro0pt;
  vector<float>   *ljvoro0eta;
  vector<float>   *ljvoro0phi;
  vector<float>   *tljvoro0pt;
  vector<float>   *tljvoro0eta;

  vector<int>     *jvoro1ncl;
  vector<float>   *jvoro1mass;
  vector<float>   *jvoro1width;
  vector<float>   *jvoro1pt;
  vector<float>   *jvoro1eta;
  vector<float>   *jvoro1phi;
  vector<float>   *tjvoro1pt;
  vector<float>   *tjvoro1eta;

  vector<int>     *ljvoro1ncl;
  vector<float>   *ljvoro1mass;
  vector<float>   *ljvoro1width;
  vector<float>   *ljvoro1pt;
  vector<float>   *ljvoro1eta;
  vector<float>   *ljvoro1phi;
  vector<float>   *tljvoro1pt;
  vector<float>   *tljvoro1eta;

  vector<int>     *jvoro2ncl;
  vector<float>   *jvoro2mass;
  vector<float>   *jvoro2width;
  vector<float>   *jvoro2pt;
  vector<float>   *jvoro2eta;
  vector<float>   *jvoro2phi;
  vector<float>   *tjvoro2pt;
  vector<float>   *tjvoro2eta;

  vector<int>     *ljvoro2ncl;
  vector<float>   *ljvoro2mass;
  vector<float>   *ljvoro2width;
  vector<float>   *ljvoro2pt;
  vector<float>   *ljvoro2eta;
  vector<float>   *ljvoro2phi;
  vector<float>   *tljvoro2pt;
  vector<float>   *tljvoro2eta;

  vector<int>     *jvoro3ncl;
  vector<float>   *jvoro3mass;
  vector<float>   *jvoro3width;
  vector<float>   *jvoro3pt;
  vector<float>   *jvoro3eta;
  vector<float>   *jvoro3phi;
  vector<float>   *tjvoro3pt;
  vector<float>   *tjvoro3eta;

  vector<int>     *ljvoro3ncl;
  vector<float>   *ljvoro3mass;
  vector<float>   *ljvoro3width;
  vector<float>   *ljvoro3pt;
  vector<float>   *ljvoro3eta;
  vector<float>   *ljvoro3phi;
  vector<float>   *tljvoro3pt;
  vector<float>   *tljvoro3eta;

  vector<int>     *jvoro4ncl;
  vector<float>   *jvoro4mass;
  vector<float>   *jvoro4width;
  vector<float>   *jvoro4pt;
  vector<float>   *jvoro4eta;
  vector<float>   *jvoro4phi;
  vector<float>   *tjvoro4pt;
  vector<float>   *tjvoro4eta;

  vector<int>     *ljvoro4ncl;
  vector<float>   *ljvoro4mass;
  vector<float>   *ljvoro4width;
  vector<float>   *ljvoro4pt;
  vector<float>   *ljvoro4eta;
  vector<float>   *ljvoro4phi;
  vector<float>   *tljvoro4pt;
  vector<float>   *tljvoro4eta;

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

