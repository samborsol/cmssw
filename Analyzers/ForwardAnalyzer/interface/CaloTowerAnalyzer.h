#ifndef CALOTOWERANALYZER_H
#define CALOTOWERANALYZER_H

// system include files
#include <string>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"



//TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Root Classes
#include "TTree.h"
#include "TFile.h"

using namespace reco;
using namespace std;

class CaloTowerAnalyzer : public edm::EDAnalyzer {
public:
	explicit CaloTowerAnalyzer(const edm::ParameterSet&);
	~CaloTowerAnalyzer();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void getCaloHits(edm::Handle<CaloTowerCollection>,vector<double>&,
			 vector<double>&,vector<double>&,vector<double>&,
			 vector<double>&,vector<double>&);

	//virtual void getCaloHits(edm::Handle<CaloTowerCollection>);

	edm::Service<TFileService> mFileServer;
	   
	string towerCollection;

	//Calo Hits
   TTree* CaloTowerTree;
   int Calo_NumberOfHits;
   vector<double> Energy,Et,Eta,Phi,EMEt,HADEt;
   /*   float Calo_Energy[10000];
   float Calo_Et[10000];
   float Calo_Eta[10000];
   float Calo_Phi[10000];
   float Calo_EMEt[10000];
   float Calo_HADEt[10000];*/


};
DEFINE_FWK_MODULE(CaloTowerAnalyzer);

#endif
