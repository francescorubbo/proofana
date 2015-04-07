#include "../scripts/runLocal.C"
#include "../scripts/runProof.C"
#include "../scripts/helperFunc.C"
#include "../scripts/helperJetMETCommon.C"
#include "../scripts/loadLibraries.C"
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

void runVoronoi_cvf(TString mode       = "cluster",         // local, lite, or cluster
TString identifier = "pileup",                      // tag 
	       // TString dataset   = "mc12_14TeV_Pythia8_J2_ITK_140_140_COMMON.jetmet2012pileupcustom",  // dataset name
	       TString dataset   = "PythJ1and2mc12aJETMET_short.jetmet2012pileupcustom",  // dataset name
	       // TString dataset   = "PythiaNoPU_COMMON.jetmet2012pileupcustom",  // dataset name
	        // TString dataset   = "PythiaPU40_COMMON.jetmet2012pileupcustom",  // dataset name
	       // TString dataset   = "PythiaPU80_COMMON.jetmet2012pileupcustom", 
	       // dataset name
TString username   = "acukierm",                               // username (e.g. swiatlow, fizisist)
bool mcweights     = false,                                 // use mc weights?
bool debug         = false,                                // debug mode
bool doPRWGen      = false,                                 // PRWgen
Long64_t nentries  = 100                           // nevents
    ) 
{ 
    
    TString date(currentDateTime().c_str());
    identifier = date+"_"+identifier;

    ///----------------------------------------------------------------
    /// Load libraries , set the config file, treenam, and cluster info
    ///----------------------------------------------------------------

    cout << "trying to load libraries" << endl; 
    loadLibraries();

    cout << " Libraries loaded " << endl;

    // SetConfig 
    TString configfile("../config/pileupstudies.config");

    
    // Best to leave alone  
    TString pathLite("");
    TString pathCluster("root://atlprf01.slac.stanford.edu:2094//atlas/output/");
    pathCluster.Append(username);
    pathCluster.Append("/");

    // Determine eventbuilder from dataset name
    TString eventbuilder(dataset);
    eventbuilder.Remove(0,eventbuilder.Last('.')+1); 
    //Change if defaulting to wrong TTree, otherwise leave
    TString treename;
    if(dataset.Contains("SMWZ") || dataset.Contains("COMMON")){
        treename = "physics";
    } else{
        treename = "qcd";
    }
   
 
    ///----------------------------------------------------------------
    /// Filename paths, URLs for PROOF running
    ///----------------------------------------------------------------
    bool runCluster(false);
    TString url(mode);
    TString path("");
    if(mode.CompareTo("lite")==0) {
        url = "lite://";
        path = pathLite;
    }
    else if(mode.CompareTo("cluster")==0) {
        url = TString(username+"@atlprf01.slac.stanford.edu");
        path = pathCluster;
        runCluster = true;
    }
    
    // Make an options file, edit as needed
    TFile* options = new TFile("options.root","RECREATE");

    ///----------------------------------------------------------------
    /// Overall Configuration
    ///----------------------------------------------------------------
    bool doBasic          = true;
    bool doJetCalibrations= true;
    bool doJetTriggers    = false;
    bool doPRW            = false;
    bool doParentChild    = false;
    bool doTrack          = true;
    bool doVertex         = true;
    bool doLCCluster      = true; 
    bool doEMCluster      = false;
    bool doEMJets         = false;
    bool doLCJets         = true;
    bool doJet4           = true;
    bool doJet6           = false; 
    bool doVectorJVFLinks = true;
    bool doTruth          = true;
    bool doTruthJets      = true;
    bool doOOTtruthJet4   = false;
    bool doTruthLinks     = true;
    bool doPhotons        = false;
    bool doElectrons      = false;
    bool doConstitLinks   = true;
    bool doMuons          = false;
    int  counterMax       = -1;
    TString prwTypes      = "EF_j280_a4tchad";

    bool doCOMMON         = false;
    bool doITKFixes       = false;
    if(dataset.Contains("COMMON")){
        doCOMMON= true;
    }
    if(dataset.Contains("ITK")){
        doITKFixes= false;
    }

    float maxJetTruthMatchdR = 0.3;
    float minJetPUdR         = 0.6;
    bool  requireTrueHSvertex= false;
    float LUMI               =1;
    bool  doSMWZfixes         = false;
    bool  doCOMMONfixes       = false;
    bool  addLeptonZInfo      = false;
    bool  doQCDSelection      = true;
    bool  doInTimeTruthJet4   = false;
    bool  addTracksToTree     = false;
    bool  doTrkFromVxt        = true;
    bool  RecoverBtracks      = true;



    ///----------------------------------------------------------------
    /// Nominal Configuration
    /// set up the actual analyses
    ///----------------------------------------------------------------
    Config* QCDSelection = new Config("QCDSelection",configfile);
    QCDSelection->Set("ANALYSIS","QCDCommonSelection");
    QCDSelection->Set("DEBUG",debug);
    
    Config* voronoi = new Config("voronoi_cvf",configfile);
    voronoi->Set("ANALYSIS","voronoi_cvf");
    voronoi->Set("DEBUG",debug);

    Config* tree0 = new Config("tree0",configfile);
    tree0->Set("ANALYSIS","tree_voronoi_cvf");
    tree0->Set("DEBUG",debug);

    
    Config* chain = new Config("chain",configfile);
    chain->AddVec("ANALYSIS");
    chain->Add("ANALYSIS",QCDSelection);
    chain->Add("ANALYSIS",voronoi);
    chain->Add("ANALYSIS",tree0);


    // set up configurations, this overwrites configs from configfile
    chain->Set("DOJETCALIBRATIONS",doJetCalibrations);
    chain->Set("DOJETTRIGGERS"   , doJetTriggers);
    chain->Set("DOJET4"          , doJet4          );
    chain->Set("DOJET6"          , doJet6          );
    chain->Set("DOVECTORJVFLINKS", doVectorJVFLinks);
    chain->Set("DOLCJETS"        , doLCJets        );
    chain->Set("DOEMJETS"        , doEMJets        );
    chain->Set("COUNTERMAX"      , counterMax      );
    chain->Set("DEBUG"           , debug           );
    chain->Set("MCWEIGHTS"       , mcweights       );
    chain->Set("PILE"            , doPRW           );
    chain->Set("DOBASIC"         , doBasic         );
    chain->Set("DOTRUTHLINKS"    , doTruthLinks    );
    chain->Set("DOCONSTITLINKS"  , doConstitLinks  );
    chain->Set("DOPARENTCHILD"   , doParentChild   );
    chain->Set("DOTRACK"         , doTrack         );
    chain->Set("DOLCCLUSTER"     , doLCCluster     );
    chain->Set("DOEMCLUSTER"     , doEMCluster     );
    chain->Set("DOTRUTH"         , doTruth         );
    chain->Set("DOTRUTHJETS"     , doTruthJets     );
    chain->Set("DOOOTTRUTHJET4"  , doOOTtruthJet4  );
    chain->Set("DOINTIMETRUTHJET4", doInTimeTruthJet4);
    chain->Set("DOVTX"           , doVertex        );
    chain->Set("DOPHOTON"        , doPhotons       );
    chain->Set("DOELECTRONS"     , doElectrons     );
    chain->Set("DOMUONS"         , doMuons         );
    chain->Set("JETTYPES"        , ""              );
    chain->Set("BTAGS"           , ""              ); 
    chain->Set("PRWTYPES"        , prwTypes        );
    chain->Set("ADDTRACKSTOTREE" , addTracksToTree );
    chain->Set("ADDLEPTONZINFOTOTREE", addLeptonZInfo);
    chain->Set("LUMI",             LUMI);
    chain->Set("requireTrueHSvertex", requireTrueHSvertex);
    chain->Set("DOSMWZfixes"       , doSMWZfixes          );
    chain->Set("DOCOMMON"          , doCOMMON          );
    chain->Set("DOCOMMONfixes"     , doCOMMONfixes          );
    chain->Set("doQCDSelection"    , doQCDSelection); 
    chain->Set("DOTRKFROMVTX",       doTrkFromVxt);
    chain->Set("runCluster",       runCluster);
    chain->Set("RecoverBtracks",   RecoverBtracks   );
    chain->Set("DOITKFIXES",   doITKFixes);
  

    // set cuts and selections
    chain->Set("MAXJETTRUTHMATCHDR"       , maxJetTruthMatchdR);
    chain->Set("MINJETPUDR"               , minJetPUdR);
    chain->Write();


    
    if(doPRWGen){
    WritePRWConfigure(options, "MC12a");
    }
    WriteJetCalibrationObjects(options);
    WriteGRLObject("data12_8TeV.periodAB_HEAD_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
    WriteDijetPRWO(options, "EF_j15_a4tchad");
    WriteDijetPRWO(options, "EF_j25_a4tchad");
    WriteDijetPRWO(options, "EF_j35_a4tchad");
    WriteDijetPRWO(options, "EF_j45_a4tchad");
    WriteDijetPRWO(options, "EF_j45_a4tchad_L2FS_L1J15");
    WriteDijetPRWO(options, "EF_j55_a4tchad");
    WriteDijetPRWO(options, "EF_j80_a4tchad");
    WriteDijetPRWO(options, "EF_j110_a4tchad");
    WriteDijetPRWO(options, "EF_j145_a4tchad");
    WriteDijetPRWO(options, "EF_j180_a4tchad");
    WriteDijetPRWO(options, "EF_j220_a4tchad");
    WriteDijetPRWO(options, "EF_j360_a4tchad");
    WriteDijetPRWO(options, "EF_j460_a4tchad");



    //if (chain->Exists("GRL")) WriteGRLObject(chain->String("GRL"));  

    ///----------------------------------------------------------------
    /// ProofAna global Config object
    ///----------------------------------------------------------------
    Config* confProofAna = new Config("ProofAna");
    
    confProofAna->Set("DEBUG"          , false        );  // "false", 0, "0" etc. also works
    confProofAna->Set("SAVETIMERS"     , false        );  // ProofAna timer histos in output file   
    confProofAna->Set("IDENTIFIER"     , identifier   );
    confProofAna->Set("DATASET"        , dataset      );
    confProofAna->Set("OUTPUTPATH"     , path         );
    confProofAna->Set("EVENTBUILDER"   , eventbuilder );
    confProofAna->Set("MERGE",true);     // enable dataset mode
   
    cout << "set eventbuilder to " << eventbuilder << endl;
 
    ///----------------------------------------------------------------
    /// Read information used in MC weighting, multi-dataset jobs
    ///----------------------------------------------------------------
    ReadDatasetInfo(dataset  , confProofAna  ); 
    //WriteGroomedPRWO(options , "EF_j240_a4tc_EFFS" );
    confProofAna->Write();
    options->Close();
    delete options;
 
 
    cout << "All setup, ready to go " << endl; 
    int runNevents=1; 
 

    // Decide to run local or on the cluster
    if(mode.CompareTo("local")==0) runLocal(dataset,treename,nentries);
    else{
        runProof(url,dataset,-1,treename);
    }
    gSystem->Unlink("options.root");

}


const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y%m%d.%H.%M", &tstruct);
    return buf;
}
