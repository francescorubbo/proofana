/**************************************************************************
 **
 **   File:         Analysis_voronoi.h
 **
 **   Description:  PileUp Analysis
 **                 
 ** 
 **   Authors:      F. Rubbo
 **
 **   Created:      11/21/2014
 **
 **************************************************************************/

#ifndef Analysis_voronoi_h
#define Analysis_voronoi_h

#include "Analysis_pileup.h"
#include <TTree.h>
#include "Particle.h"
#include "TMVA/Reader.h"
 
using std::cout;
using std::endl;

class Analysis_voronoi : public Analysis_pileup {

 public :
  
  Analysis_voronoi(TTree* /*tree*/ = 0) { 
    fDetail = false; 
  }
  
  virtual ~Analysis_voronoi() { }
  
  ClassDef(Analysis_voronoi, 0);
  
  Bool_t  fDetail;
  
  virtual bool    ProcessEvent();
  virtual void    WorkerBegin(); 
  virtual void    WorkerTerminate();

  MomKey MakeVoronoiClusters(const fastjet::JetAlgorithm algo, const double jetR, const MomKey constType);
  void selectClusters();

  private :			  

};

#endif

