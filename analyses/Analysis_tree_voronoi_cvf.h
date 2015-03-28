/**************************************************************************
 **
 **   File:         Analysis_tree_voronoi_cvf.h
 **
 **   Description:  Tree filler
 **
 **
 **   Authors:      F. Rubbo
 **
 **   Created:      11/25/2014
 **
 **************************************************************************/

#ifndef Analysis_tree_voronoi_cvf_h
#define Analysis_tree_voronoi_cvf_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"

using std::cout;
using std::endl;


class Analysis_tree_voronoi_cvf : public Analysis_JetMET_Base {

 public :

  Analysis_tree_voronoi_cvf(TTree* /*tree*/ = 0) {
    fDetail = false;
  }

  virtual ~Analysis_tree_voronoi_cvf() { }

  ClassDef(Analysis_tree_voronoi_cvf, 0);

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
  vector<float>     *j0mass;
  vector<float>     *j0width;
  vector<float>   *j0pt;
  vector<float>   *j0eta;
  vector<float>   *j0phi;
  vector<float>   *j0jvt;
  vector<float>   *tj0pt;
  vector<float>   *tj0eta;

  vector<int>     *j0cvfncl;
  vector<float>     *j0cvfmass;
  vector<float>     *j0cvfwidth;
  vector<float>   *j0cvfpt;
  vector<float>   *j0cvfeta;
  vector<float>   *j0cvfphi;
  vector<float>   *j0cvfjvt;
  vector<float>   *tj0cvfpt;
  vector<float>   *tj0cvfeta;

  vector<int>     *jnoarea0ncl;
  vector<float>     *jnoarea0mass;
  vector<float>     *jnoarea0width;
  vector<float>   *jnoarea0pt;
  vector<float>   *jnoarea0eta;
  vector<float>   *jnoarea0phi;
  vector<float>   *tjnoarea0pt;
  vector<float>   *tjnoarea0eta;

  vector<int>     *jvoro0ncl;
  vector<float>   *jvoro0mass;
  vector<float>   *jvoro0width;
  vector<float>   *jvoro0pt;
  vector<float>   *jvoro0eta;
  vector<float>   *jvoro0phi;
  vector<float>   *tjvoro0pt;
  vector<float>   *tjvoro0eta;

  vector<int>     *jvoro1ncl;
  vector<float>   *jvoro1mass;
  vector<float>   *jvoro1width;
  vector<float>   *jvoro1pt;
  vector<float>   *jvoro1eta;
  vector<float>   *jvoro1phi;
  vector<float>   *tjvoro1pt;
  vector<float>   *tjvoro1eta;

  vector<int>     *jvoro0cvf5ncl;
  vector<float>   *jvoro0cvf5mass;
  vector<float>   *jvoro0cvf5width;
  vector<float>   *jvoro0cvf5pt;
  vector<float>   *jvoro0cvf5eta;
  vector<float>   *jvoro0cvf5phi;
  vector<float>   *tjvoro0cvf5pt;
  vector<float>   *tjvoro0cvf5eta;

  vector<int>     *jvoro1cvf5ncl;
  vector<float>   *jvoro1cvf5mass;
  vector<float>   *jvoro1cvf5width;
  vector<float>   *jvoro1cvf5pt;
  vector<float>   *jvoro1cvf5eta;
  vector<float>   *jvoro1cvf5phi;
  vector<float>   *tjvoro1cvf5pt;
  vector<float>   *tjvoro1cvf5eta;

  vector<int>     *jvoro0cvfxncl;
  vector<float>   *jvoro0cvfxmass;
  vector<float>   *jvoro0cvfxwidth;
  vector<float>   *jvoro0cvfxpt;
  vector<float>   *jvoro0cvfxeta;
  vector<float>   *jvoro0cvfxphi;
  vector<float>   *tjvoro0cvfxpt;
  vector<float>   *tjvoro0cvfxeta;

  vector<int>     *jvoro1cvfxncl;
  vector<float>   *jvoro1cvfxmass;
  vector<float>   *jvoro1cvfxwidth;
  vector<float>   *jvoro1cvfxpt;
  vector<float>   *jvoro1cvfxeta;
  vector<float>   *jvoro1cvfxphi;
  vector<float>   *tjvoro1cvfxpt;
  vector<float>   *tjvoro1cvfxeta;

  vector<int>     *jvoro10ncl;
  vector<float>   *jvoro10mass;
  vector<float>   *jvoro10width;
  vector<float>   *jvoro10pt;
  vector<float>   *jvoro10eta;
  vector<float>   *jvoro10phi;
  vector<float>   *tjvoro10pt;
  vector<float>   *tjvoro10eta;

  vector<int>     *jvorosncl;
  vector<float>   *jvorosmass;
  vector<float>   *jvoroswidth;
  vector<float>   *jvorospt;
  vector<float>   *jvoroseta;
  vector<float>   *jvorosphi;
  vector<float>   *tjvorospt;
  vector<float>   *tjvoroseta;

  vector<int>     *jvorocvf5sncl;
  vector<float>   *jvorocvf5smass;
  vector<float>   *jvorocvf5swidth;
  vector<float>   *jvorocvf5spt;
  vector<float>   *jvorocvf5seta;
  vector<float>   *jvorocvf5sphi;
  vector<float>   *tjvorocvf5spt;
  vector<float>   *tjvorocvf5seta;

  vector<int>     *jvorocvfxsncl;
  vector<float>   *jvorocvfxsmass;
  vector<float>   *jvorocvfxswidth;
  vector<float>   *jvorocvfxspt;
  vector<float>   *jvorocvfxseta;
  vector<float>   *jvorocvfxsphi;
  vector<float>   *tjvorocvfxspt;
  vector<float>   *tjvorocvfxseta;

  vector<float>   *truejetpt;
  vector<float>   *truejeteta;
  vector<float>   *truejetphi;
  vector<float>   *truejetenergy;

  private :

};

#endif

