// -*- C++ -*-
//
// Package:    ForwardAnalyzer
// Class:      ForwardAnalyzer
//
/**\class ForwardAnalyzer ForwardAnalyzer.cc Analyzers/ForwardAnalyzer/src/ForwardAnalyzer.cc

Description: A simple analyzer which can read in Reco or RAW data and makes trees with Reco information and Digi information as well for all of the forward detectors.

90% of this was shamelessly stolen from Jaime Gomez's (Maryland) analyzer, 90% of which was itself stolen from Pat Kenny's (Kansas) analyzer

Implementation:
This can be run out of the box with the test cfg or can be installed into a crab.cfg
*/

#include "Analyzers/ForwardAnalyzer/interface/ForwardAnalyzer.h"
#include <iostream>

//static const float HFQIEConst = 4.0;
//static const float HFQIEConst = 2.6;
//static const float EMGain  = 0.025514;
//static const float HADGain = 0.782828;

using namespace edm;
using namespace std;

ForwardAnalyzer::ForwardAnalyzer(const ParameterSet& iConfig) 
//CentralityTag_(consumes<reco::Centrality>(iConfig.getParameter<edm::InputTag>("CentralitySrc"))),
//CentralityBinTag_(consumes<int>(iConfig.getParameter<edm::InputTag>("CentralityBinSrc")))
//:
//selectedBins_(iConfig.getParameter<std::vector<int> >("selectedBins"))
{
    runBegin = -1;
    evtNo = 0;
    lumibegin = 0;
    lumiend = 0;
    startTime = "Not Avaliable";
    
    // Get NTracks
    using namespace edm;
    vtxCollection_  = iConfig.getParameter<edm::InputTag>("vtxCollection_");
    trackSrc = iConfig.getParameter<edm::InputTag>("trackSrc");
}

ForwardAnalyzer::~ForwardAnalyzer(){}

void ForwardAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
{
    using namespace edm;
    using namespace reco;
    using namespace std;

    ++evtNo;
    time_t a = (iEvent.time().value()) >> 32; // store event info
    event = iEvent.id().event();
    run = iEvent.id().run();
    lumi = iEvent.luminosityBlock();
    edm::ESHandle<HcalTopology> htopo;
    iSetup.get<IdealGeometryRecord>().get(htopo);
    
    // get NTracks and centrality
    //edm::Handle<int> cbin_;
    //iEvent.getByToken(CentralityBinTag_,cbin_);
    //hiBin = *cbin_; //HF tower centrality
    //CentHFtowers=hiBin;
//    std::cout<<"CentHFtowers = "<<CentHFtowers<<"\n"<<std::endl;
    
    //edm::Handle<reco::Centrality> centrality;
    //iEvent.getByToken(CentralityTag_, centrality);
    
    //hiNpix = centrality->multiplicityPixel();
    //hiNpixelTracks = centrality->NpixelTracks();
    //hiNtracks = centrality->Ntracks();
//    hiNtracksPtCut = centrality->NtracksPtCut();
//    hiNtracksEtaCut = centrality->NtracksEtaCut();
//    hiNtracksEtaPtCut = centrality->NtracksEtaPtCut();
    
    //hiHF = centrality->EtHFtowerSum();
//    hiHFplus = centrality->EtHFtowerSumPlus();
//    hiHFminus = centrality->EtHFtowerSumMinus();
//    hiHFplusEta4 = centrality->EtHFtruncatedPlus();
//    hiHFminusEta4 = centrality->EtHFtruncatedMinus();

    //hiHFhit = centrality->EtHFhitSum();
//    hiHFhitPlus = centrality->EtHFhitSumPlus();
//    hiHFhitMinus = centrality->EtHFhitSumMinus();
    
//    hiEEplus = centrality->EtEESumPlus();
//    hiEEminus = centrality->EtEESumMinus();
//    hiEE = centrality->EtEESum();
//    hiEB = centrality->EtEBSum();
//    hiET = centrality->EtMidRapiditySum();
    
    /*
    std::cout<<"hiBin = "<<hiBin<<"\n"<<std::endl;
    
    std::cout<<"hiNpix = "<<hiNpix<<std::endl;
    std::cout<<"hiNpixelTracks = "<<hiNpixelTracks<<std::endl;
    std::cout<<"hiNtracks = "<<hiNtracks<<std::endl;
    std::cout<<"hiNtracksPtCut = "<<hiNtracksPtCut<<std::endl;
    std::cout<<"hiNtracksEtaCut = "<<hiNtracksEtaCut<<std::endl;
    std::cout<<"hiNtracksEtaPtCut = "<<hiNtracksEtaPtCut<<"\n"<<std::endl;

    cout<<"hiHF = "<<hiHF<<std::endl;
    cout<<"hiHFplus = "<<hiHFplus<<std::endl;
    cout<<"hiHFminus = "<<hiHFminus<<std::endl;
    cout<<"hiHFplusEta4 = "<<hiHFplusEta4<<std::endl;
    cout<<"hiHFminusEta4 = "<<hiHFminusEta4<<std::endl;
    cout<<"hiHFhit = "<<hiHFhit<<std::endl;
    cout<<"hiHFhitPlus = "<<hiHFhitPlus<<std::endl;
    cout<<"hiHFhitMinus = "<<hiHFhitMinus<<"\n"<<std::endl;
    */

/*    int noff = getNoff( iEvent, iSetup);
    if ( (Noff < Noffmin_) or (Noff >= Noffmax_) ) {
        return;
    }
    Noff = noff;
//    std::cout<<"Noff = "<<Noff<<"\n"<<std::endl; */

    // Get vertex
    edm::Handle<reco::VertexCollection> vertexCollection3;
    iEvent.getByLabel(vtxCollection_,vertexCollection3);
//    if(vzr_sell<minvtx_ || vzr_sell>maxvtx_) return;
    const reco::VertexCollection * vertices3 = vertexCollection3.product();
    vs_sell = vertices3->size();
    if(vs_sell>0) {
        vzr_sell = vertices3->begin()->z();
        vzErr_sell = vertices3->begin()->zError();
    } else
        vzr_sell = -999.9;
    vtxval = vzr_sell;
    //CentTree->Fill();

    if (runBegin < 0) {         // parameters for the first event
        startTime = ctime(&a);
        lumibegin = lumiend = lumi;
        runBegin = iEvent.id().run();
        Runno=iEvent.id().run();
    }
    if (lumi < lumibegin)lumibegin = lumi;
    if (lumi > lumiend)lumiend = lumi;

    BeamData[0]=iEvent.bunchCrossing();
    BeamData[1]=lumi;
    BeamData[2]=run;
    BeamData[3]=event;
    
    BeamTree->Fill();

    Handle<ZDCDigiCollection> castorDigis;

    edm::Handle<HFRecHitCollection> hf_recHits_h;
    Handle<CastorRecHitCollection> castor_recHits_h;
    edm::ESHandle<HcalDbService> conditions;
    iSetup.get<HcalDbRecord>().get(conditions);

    iEvent.getByLabel("castor_recHits_h",castor_recHits_h);
  

    if (!(iEvent.getByLabel("hfreco",hf_recHits_h)))
        {
            std::cout<<"NO HF, I don't care!!!"<<std::endl;
   //   return;
        }

    if (!(iEvent.getByLabel("castorDigis",castorDigis)))
    {
        std::cout<<"No ZDC Digis gotten with label:"<<castorDigis<<std::endl;
    //    return;
    }
    const ZDCDigiCollection *zdc_digi = castorDigis.failedToGet()? 0 : &*castorDigis;
    const HFRecHitCollection *hf_recHits = hf_recHits_h.failedToGet()? 0 : &*hf_recHits_h;

    iSetup.get<HcalDbRecord>().get(conditions);

    if(zdc_digi){
        for(int i=0; i<180; i++){DigiDatafC[i]=0;DigiDataADC[i]=0;}
//std::cout<<"ForwardAnalyzer::analyze if(zdc_digi) is fine, now to loop over the digis"<<std::endl;
        for (ZDCDigiCollection::const_iterator j=zdc_digi->begin();j!=zdc_digi->end();j++){
            const ZDCDataFrame digi = (const ZDCDataFrame)(*j);
            int iSide      = digi.id().zside();
            int iSection   = digi.id().section();
            int iChannel   = digi.id().channel();
            int chid = (iSection-1)*5+(iSide+1)/2*9+(iChannel-1);
      
//      std::cout << "ForwardAnalyzer::analyze, looping over digis, at j: "<<j<<std::endl;
//      std::cout<<"ForwardAnalyzer::analyze looping over digis, iSide: "<<iSide<<", iSection: "<<iSection<<", iChannel: "<<iChannel<<", chid: "<<chid<<std::endl;
      //<<

            HcalZDCDetId cell = j->id();
//	const HcalCalibrations& calibrations=conditions->getHcalCalibrations(cell);
            const HcalQIECoder* qiecoder = conditions->getHcalCoder (cell);
            const HcalQIEShape* qieshape = conditions->getHcalShape (qiecoder);
//	HcalCoderDb coder (*qiecoder, *qieshape);

//      const HcalQIEShape* qieshape=conditions->getHcalShape();
//      const HcalQIECoder* qiecoder=conditions->getHcalCoder(digi.id());
            CaloSamples caldigi;
            HcalCoderDb coder(*qiecoder,*qieshape);

            coder.adc2fC(digi,caldigi);

            int fTS = digi.size();
            for (int i = 0; i < fTS; ++i) {
                DigiDatafC[i+chid*10] = caldigi[i];
                DigiDataADC[i+chid*10] = digi[i].adc();
            }
        }

        ZDCDigiTree->Fill();
    }


    if(hf_recHits){
        int hf_counter=0;
        HF_NumberofHits=hf_recHits->size();
        for (HFRecHitCollection::const_iterator HFiter=hf_recHits->begin();HFiter!=hf_recHits->end();++HFiter){
            HF_Energy[hf_counter]=HFiter->energy();
            HcalDetId id(HFiter->detid().rawId());
            HF_iEta[hf_counter]=id.ieta();
            HF_iPhi[hf_counter]=id.iphi();
            HF_Depth[hf_counter]=id.depth();
            ++hf_counter;
        } // end of HF RecHits Iterator
        HFRecoTree->Fill();
    } // end of if(hf_recHits)
    
    
} // end of Analyze


void ForwardAnalyzer::beginJob(){
    mFileServer->file().SetCompressionLevel(9);
    mFileServer->file().cd();
    
    hiBin = -1;

    string bnames[] = {"negEM1","negEM2","negEM3","negEM4","negEM5",
                        "negHD1","negHD2","negHD3","negHD4",
                        "posEM1","posEM2","posEM3","posEM4","posEM5",
                        "posHD1","posHD2","posHD3","posHD4"};
    BranchNames=bnames;
    
    ZDCDigiTree = new TTree("ZDCDigiTree","ZDC Digi Tree");
    CentTree = new TTree("CentTree","Centrality Tree");
    BeamTree = new TTree("BeamTree","Beam Tree");
    HFRecoTree = new TTree("HFRecTree","HF RecHit Tree");

    for(int i=0; i<18; i++){
        ZDCDigiTree->Branch((bnames[i]+"fC").c_str(),&DigiDatafC[i*10],(bnames[i]+"cFtsz[10]/F").c_str());
        ZDCDigiTree->Branch((bnames[i]+"ADC").c_str(),&DigiDataADC[i*10],(bnames[i]+"ADCtsz[10]/I").c_str());
    }
    
    //CentTree->Branch("Cent_HFtowers",&CentHFtowers,"Cent_HFtowers/I");
//    CentTree->Branch("NTrackOff",&Noff,"NTracksOff/I");
    //CentTree->Branch("Vtx",&vtxval,"Vtx/D");
    
    //CentTree->Branch("hiHF",&hiHF,"hiHF/F");
//    CentTree->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
//    CentTree->Branch("hiHFminus",&hiHFminus,"hiHFminus/F");
//    CentTree->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
//    CentTree->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F");
    
    //CentTree->Branch("hiHFhit",&hiHFhit,"hiHFhit/F");
//    CentTree->Branch("hiHFhitPlus",&hiHFhitPlus,"hiHFhitPlus/F");
//    CentTree->Branch("hiHFhitMinus",&hiHFhitMinus,"hiHFhitMinus/F");
    
//    CentTree->Branch("hiET",&hiET,"hiET/F");
//    CentTree->Branch("hiEE",&hiEE,"hiEE/F");
//    CentTree->Branch("hiEB",&hiEB,"hiEB/F");
//    CentTree->Branch("hiEEplus",&hiEEplus,"hiEEplus/F");
//    CentTree->Branch("hiEEminus",&hiEEminus,"hiEEminus/F");
    //CentTree->Branch("hiNpix",&hiNpix,"hiNpix/I");
    //CentTree->Branch("hiNpixelTracks",&hiNpixelTracks,"hiNpixelTracks/I");
    //CentTree->Branch("hiNtracks",&hiNtracks,"hiNtracks/I");
//    CentTree->Branch("hiNtracksPtCut",&hiNtracksPtCut,"hiNtracksPtCut/I");
//    CentTree->Branch("hiNtracksEtaCut",&hiNtracksEtaCut,"hiNtracksEtaCut/I");
//    CentTree->Branch("hiNtracksEtaPtCut",&hiNtracksEtaPtCut,"hiNtracksEtaPtCut/I");

    HFRecoTree->Branch("HF_NumberOfHits",&HF_NumberofHits,"HF_NumberOfHits/I");
    HFRecoTree->Branch("HF_iEta",&HF_iEta,"HF_iEta[HF_NumberOfHits]/F");
    HFRecoTree->Branch("HF_iPhi",&HF_iPhi,"HF_iPhi[HF_NumberOfHits]/F");
    HFRecoTree->Branch("HF_Depth",&HF_Depth,"HF_Depth[HF_NumberOfHits]/F");
    HFRecoTree->Branch("HF_Energy",&HF_Energy,"HF_Energy[HF_NumberOfHits]/F");

    BeamTree->Branch("BunchXing",&BeamData[0],"BunchXing/I");
    BeamTree->Branch("LumiBlock",&BeamData[1],"LumiBlock/I");
    BeamTree->Branch("Run",&BeamData[2],"Run/I");
    BeamTree->Branch("Event",&BeamData[3],"Event/I");
}
