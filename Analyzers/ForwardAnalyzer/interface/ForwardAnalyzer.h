#ifndef FORWARDANALYZER_H
#define FORWARDANALYZER_H

// system include files
#include <string>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/Common/interface/Handle.h"

//newly added headers
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/ZdcSimpleRecAlgo.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalHFStatusBitFromRecHits.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalHFStatusBitFromDigis.h"
#include "DataFormats/METReco/interface/HcalCaloFlagLabels.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/HcalObjects/interface/HcalChannelStatus.h"
#include "CondFormats/HcalObjects/interface/HcalLongRecoParams.h"
#include "CondFormats/HcalObjects/interface/HcalLongRecoParam.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HBHEStatusBitSetter.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalTimingCorrector.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HBHETimeProfileStatusBitSetter.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HBHETimingShapedFlag.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalADCSaturationFlag.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HFTimingTrustFlag.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/EDCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

//TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// L1 bit analyzer
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

// new for NTracks
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"

//Root Classes
#include "TTree.h"
#include "TFile.h"

class ForwardAnalyzer : public edm::EDAnalyzer {
public:
	explicit ForwardAnalyzer(const edm::ParameterSet&);
	~ForwardAnalyzer();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);

	long runBegin,lumibegin,lumiend,evtNo;
	int run, event, lumi;
	std::string startTime;
	std::string *BranchNames;
        int DigiDataADC[180];
        float DigiDatafC[180];
	int BeamData[6];
    
	int Runno;
	edm::Service<TFileService> mFileServer;
	   
    TTree* ZDCDigiTree;
    TTree* CentTree;
    TTree* BeamTree;
   
    //New HF Information
    TTree* HFRecoTree;
    int HF_NumberofHits;
    float HF_iEta[2000];
    float HF_iPhi[2000];
    float HF_Depth[2000];
    float HF_Energy[2000];
    //End New HF information,
    
    // Get centrality
    edm::EDGetTokenT<reco::Centrality> CentralityTag_;
    edm::EDGetTokenT<int> CentralityBinTag_;
    
    int hiBin;
    
//    int hiNpix, hiNpixelTracks, hiNtracks, hiNtracksPtCut, hiNtracksEtaCut, hiNtracksEtaPtCut;
//    float hiHF, hiHFplus, hiHFminus, hiHFplusEta4, hiHFminusEta4, hiHFhit, hiHFhitPlus, hiHFhitMinus, hiEB, hiET, hiEE, hiEEplus, hiEEminus;
    int hiNpix, hiNpixelTracks, hiNtracks;
    float hiHF, hiHFhit;
    
    std::vector<int> selectedBins_;
    edm::Handle<int> cbin_;
    edm::EDGetTokenT<int> tag_;
    int CentHFtowers;
//    int Noff;

    // Get NTracks
    edm::InputTag trackSrc;
    edm::InputTag vtxCollection_;
//    double minvtx_;
//    double maxvtx_;
//    int Noffmin_;
//    int Noffmax_;
    int vs_sell;   // vertex collection size
    float vzr_sell;
    float vzErr_sell;
    double vtxval;
    int nEtaBins;

};
DEFINE_FWK_MODULE(ForwardAnalyzer);

#endif
