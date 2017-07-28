#include "Analyzers/ForwardAnalyzer/interface/CaloTowerAnalyzer.h"

using namespace edm;

CaloTowerAnalyzer::CaloTowerAnalyzer(const edm::ParameterSet& iConfig):towerCollection(iConfig.getParameter<string>("towerCollection")){}

CaloTowerAnalyzer::~CaloTowerAnalyzer(){}

void CaloTowerAnalyzer::beginJob(){
  mFileServer->file().SetCompressionLevel(9);
  mFileServer->file().cd();

  string tName(towerCollection+"Tree");
  CaloTowerTree = new TTree(tName.c_str(),tName.c_str());

  CaloTowerTree->Branch("Calo_NumberOfHits",&Calo_NumberOfHits,"Calo_NumberOfHits/I");
  CaloTowerTree->Branch("Energy",&Energy[0],"Energy[Calo_NumberOfHits]/D");
  CaloTowerTree->Branch("Et",&Et[0],"Et[Calo_NumberOfHits]/D");
  CaloTowerTree->Branch("Eta",&Eta[0],"Eta[Calo_NumberOfHits]/D");
  CaloTowerTree->Branch("Phi",&Phi[0],"Phi[Calo_NumberOfHits]/D");
  CaloTowerTree->Branch("EMEt",&EMEt[0],"EMEt[Calo_NumberOfHits]/D");
  CaloTowerTree->Branch("HADEt",&HADEt[0],"HADEt[Calo_NumberOfHits]/D");
}

void CaloTowerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Handle<CaloTowerCollection> calo_h;
  iEvent.getByLabel("towerMaker", calo_h);


  Energy.clear(); Et.clear();
    Eta.clear(); Phi.clear(); EMEt.clear();
    HADEt.clear();

  if(!calo_h.failedToGet()){getCaloHits(calo_h,Energy,Et,Eta,Phi,EMEt,HADEt);}
  //if(!calo_h.failedToGet()){getCaloHits(calo_h);}

  Calo_NumberOfHits=Energy.size();


    CaloTowerTree->SetBranchAddress("Calo_NumberOfHits",&Calo_NumberOfHits);
      CaloTowerTree->SetBranchAddress("Energy",&Energy[0]);
      CaloTowerTree->SetBranchAddress("Et",&Et[0]);
      CaloTowerTree->SetBranchAddress("Eta",&Eta[0]);
      CaloTowerTree->SetBranchAddress("Phi",&Phi[0]);
      CaloTowerTree->SetBranchAddress("EMEt",&EMEt[0]);
      CaloTowerTree->SetBranchAddress("HADEt",&HADEt[0]);


  CaloTowerTree->Fill();
}

void CaloTowerAnalyzer::getCaloHits(Handle<CaloTowerCollection> TowerCol, vector<double> &Energy,
                          vector<double> &Et, vector<double> &Eta,
                          vector<double> &Phi, vector<double> &EMEt,
                          vector<double> &HADEt){
//void CaloTowerAnalyzer::getCaloHits(Handle<CaloTowerCollection> TowerCol){
  for(CaloTowerCollection::const_iterator cal=(&*TowerCol)->begin();cal!=(&*TowerCol)->end();cal++){
    if((fabs(cal->eta())>3) && (fabs(cal->eta())<5.5))
    {
     Energy.push_back(cal->energy());
    Et.push_back(cal->et());
     Eta.push_back(cal->eta());
     Phi.push_back(cal->phi());
    EMEt.push_back(cal->emEt());
    HADEt.push_back(cal->hadEt());
    }
  }
}
