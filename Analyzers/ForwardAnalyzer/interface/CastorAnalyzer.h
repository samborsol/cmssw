#ifndef CASTORANALYZER_H
#define CASTORANALYZER_H

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
//#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "DataFormats/HcalDetId/interface/HcalCastorDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"



//Root Classes
#include "TTree.h"
#include "TFile.h"

#include <vector>
#include "Math/Vector3D.h"

using std::vector;

class CastorAnalyzer : public edm::EDAnalyzer {
public:
	explicit CastorAnalyzer(const edm::ParameterSet&);
	~CastorAnalyzer();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
private:

	edm::Service<TFileService> mFileServer;

   //Castor
   TTree* CastorRecoTree;
   int NumberOfHits;
   vector<float> Section;
   vector<float> Sector;
   vector<float> Module;
   vector<float> Energy;
};
DEFINE_FWK_MODULE(CastorAnalyzer);

#endif
