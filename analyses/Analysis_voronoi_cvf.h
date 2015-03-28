/**************************************************************************
 **
 **   File:         Analysis_voronoi_cvf.h
 **
 **   Description:  PileUp Analysis
 **                 
 ** 
 **   Authors:      F. Rubbo
 **
 **   Created:      11/21/2014
 **
 **************************************************************************/

#ifndef Analysis_voronoi_cvf_h
#define Analysis_voronoi_cvf_h

#include "Analysis_pileup.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"
#include "ANN/ANN.h"
 
using std::cout;
using std::endl;

typedef ANNcoord* ANNpoint;
typedef ANNpoint* ANNpointArray;
typedef ANNdist* ANNdistArray;
typedef ANNidx* ANNidxArray;

class Analysis_voronoi_cvf : public Analysis_pileup {

 public :
  
  Analysis_voronoi_cvf(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_voronoi_cvf() { }
  
  ClassDef(Analysis_voronoi_cvf, 0);
  
  Bool_t  fDetail;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();

  MomKey MakeVoronoiClusters(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType);
  void selectClusters(MomKey clusterskey,int nsigma=0,bool doCVF=false,float threshold=-1);
  void selectLCTopoCVFClusters();
  void selectVoronoiCVFClusters(float threshold);
  void MakeSpreadVoronoiClusters(float spreadr, MomKey key);
  void SetCVF(MomKey key);
  void SpreadPT();

  private :			  

};

#endif

