/**************************************************************************
 **
 **   File:         Analysis_btobjvt.h
 **
 **   Description:  PileUp Analysis
 **                 
 ** 
 **   Authors:      F. Rubbo
 **
 **   Created:      11/21/2014
 **
 **************************************************************************/

#ifndef Analysis_btobjvt_h
#define Analysis_btobjvt_h

#include "Analysis_JetMET_Base.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"
#include "JetVertexTagger/JetVertexTagger.h" 
#include "ANN/ANN.h"
 
using std::cout;
using std::endl;

typedef ANNcoord* ANNpoint;
typedef ANNpoint* ANNpointArray;
typedef ANNdist* ANNdistArray;
typedef ANNidx* ANNidxArray;

class Analysis_btobjvt : public Analysis_JetMET_Base {

 public :
  
  Analysis_btobjvt(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_btobjvt() { }
  
  ClassDef(Analysis_btobjvt, 0);
  
  Bool_t  fDetail;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();

  void associateTrackstoCluster(Particle *thecluster);
  void selectTracks();
  void selectClusters(float jvfcut,string suffix);
  void addTruthMatch(const MomKey JetType, const MomKey TruthJetType);

  private :			  

  JetVertexTagger* jvt; 


};

#endif

