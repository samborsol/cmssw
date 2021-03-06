#ifndef ECAL_PULSESHAPES_HANDLER_H
#define ECAL_PULSESHAPES_HANDLER_H

#include <vector>
#include <typeinfo>
#include <string>
#include <map>
#include <iostream>
#include <time.h>

#include "CondCore/PopCon/interface/PopConSourceHandler.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"

#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/DataRecord/interface/EcalPulseShapesRcd.h"

#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/Provenance/interface/Timestamp.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

namespace popcon
{

  class EcalPulseShapesHandler : public popcon::PopConSourceHandler<EcalPulseShapes>
    {
      
    public:
      EcalPulseShapesHandler(edm::ParameterSet const & );
      ~EcalPulseShapesHandler();
      bool checkPulseShape(EcalPulseShapes::Item* item);
      void fillSimPulseShape(EcalPulseShapes::Item* item, bool isbarrel);
      void getNewObjects();
      std::string id() const { return m_name;}

    private:
      const EcalPulseShapes * mypulseshapes;

      unsigned int m_firstRun ;
      unsigned int m_lastRun ;
      
      std::string m_gentag;
      std::string m_filename;
      std::string m_name;      
      std::vector<double> m_EBPulseShapeTemplate, m_EEPulseShapeTemplate;

    };
}
#endif
