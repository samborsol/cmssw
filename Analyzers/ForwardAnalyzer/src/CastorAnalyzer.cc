#include "Analyzers/ForwardAnalyzer/interface/CastorAnalyzer.h"
#include <iostream>

using namespace edm;
using namespace std;

CastorAnalyzer::CastorAnalyzer(const ParameterSet& iConfig)
{}

CastorAnalyzer::~CastorAnalyzer(){}

void CastorAnalyzer::beginJob(){
  mFileServer->file().SetCompressionLevel(9);
  mFileServer->file().cd();

  CastorRecoTree = new TTree("CastorRecoTree","CASTOR RecHit Tree");


  CastorRecoTree->Branch("NumberOfHits",&NumberOfHits,"NumberOfHits/I");
  CastorRecoTree->Branch("Energy",&Energy[0],"Energy[NumberOfHits]/F");
  CastorRecoTree->Branch("Module",&Module[0],"Module[NumberOfHits]/F");
  CastorRecoTree->Branch("Section",&Section[0],"Section[NumberOfHits]/F");
  CastorRecoTree->Branch("Sector",&Sector[0],"Sector[NumberOfHits]/F");  
}

void CastorAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;


  Handle<CastorRecHitCollection> castor_recHits_h;
  ESHandle<HcalDbService> conditions;
  //iEvent.getByType(castor_recHits_h);  //getByType eliminated as of 2013, needs to be replaced with something new
  iEvent.getByLabel("castor_recHits_h",castor_recHits_h); //OK! This now compiles (Nov 19 2015) but it remains to be seen if it works.
  const CastorRecHitCollection *castor_recHits = castor_recHits_h.failedToGet()? 0 : &*castor_recHits_h;
  iSetup.get<HcalDbRecord>().get(conditions);

NumberOfHits=0;
Energy.clear();
Sector.clear();
Module.clear();
Section.clear();


  ////////////////////NEW CASTOR STUFF//////////////////////
  ///Castor RecHits

  if(castor_recHits){
    for (CastorRecHitCollection::const_iterator castor_iter=castor_recHits->begin();castor_iter!=castor_recHits->end();++castor_iter){
      Energy.push_back(castor_iter->energy());
      HcalCastorDetId id(castor_iter->detid().rawId());
      Section.push_back(id.section());//1=EM 2=HAD
      Module.push_back(id.module());//modules 1-2EM, 1-12HAD
      Sector.push_back(id.sector());//Sectors 1-16
     }//end of CASTOR RecHits Iterator
  }//end of if(Castor_RecHits)

NumberOfHits=Energy.size();

  CastorRecoTree->SetBranchAddress("NumberOfHits",&NumberOfHits);
  CastorRecoTree->SetBranchAddress("Energy",&Energy[0]);
  CastorRecoTree->SetBranchAddress("Sector",&Sector[0]);
  CastorRecoTree->SetBranchAddress("Module",&Module[0]);
  CastorRecoTree->SetBranchAddress("Section",&Section[0]);

CastorRecoTree->Fill();

}//end of Analyze
