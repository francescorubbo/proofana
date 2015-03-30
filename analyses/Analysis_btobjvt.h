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
#include "Analysis_pileup.h"
 
using std::cout;
using std::endl;


class Analysis_btobjvt : public Analysis_pileup {

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

  MomKey  MakeJets(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType, const MomKey extra = "",double minpt=10);
  void selectTracks();

  private :			  

  JetVertexTagger* jvt; 


};

#endif

