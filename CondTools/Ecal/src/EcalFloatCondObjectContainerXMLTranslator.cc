#include <iostream>
#include <sstream>
#include <fstream>

#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include "FWCore/Concurrency/interface/Xerces.h"
#include "Utilities/Xerces/interface/XercesStrUtils.h"
#include <xercesc/util/XMLString.hpp>
#include <xercesc/sax/SAXException.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>

#include "CondTools/Ecal/interface/EcalFloatCondObjectContainerXMLTranslator.h"
#include "CondTools/Ecal/interface/DOMHelperFunctions.h"

using namespace XERCES_CPP_NAMESPACE;
using namespace xuti;
using namespace std;


int  
EcalFloatCondObjectContainerXMLTranslator::readXML(const string& filename,
					      EcalCondHeader&          header,
					      EcalFloatCondObjectContainer& record){

  cms::concurrency::xercesInitialize();

  XercesDOMParser* parser = new XercesDOMParser;
  parser->setValidationScheme( XercesDOMParser::Val_Never );
  parser->setDoNamespaces( false );
  parser->setDoSchema( false );
  
  parser->parse(filename.c_str());

  DOMDocument* xmlDoc = parser->getDocument();
  

  if (!xmlDoc) {
    std::cout << "EcalFloatCondObjectContainerXMLTranslator::Error parsing document" << std::endl;
    return -1;
  }

  DOMElement* elementRoot = xmlDoc->getDocumentElement();

  xuti::readHeader(elementRoot, header);

  // get the first cell node
  DOMNode * cellnode=getChildNode(elementRoot,Cell_tag);
  
  // loop on cell nodes
  while  (cellnode){

    float val=0;

    // read id
    DetId detid= readCellId(dynamic_cast<DOMElement*>(cellnode));
       
    // read value
    DOMNode * c_node = getChildNode(cellnode,Value_tag);
    GetNodeData(c_node,val);

    // fill record
    record[detid]=val;

    // get next cell
    cellnode= cellnode->getNextSibling();
    
    while (cellnode&& cellnode->getNodeType( ) != DOMNode::ELEMENT_NODE)      
      cellnode= cellnode->getNextSibling();
    
    
  }


  delete parser;
  cms::concurrency::xercesTerminate();
  return 0;
    
}



std::vector<float>
EcalFloatCondObjectContainerXMLTranslator::barrelfromXML(const string& filename){
  EcalCondHeader header;  
  EcalFloatCondObjectContainer record;
  readXML(filename,header,record);
  
  return record.barrelItems();
 
}


std::vector<float>
EcalFloatCondObjectContainerXMLTranslator::endcapfromXML(const string& filename){
  EcalCondHeader header;  
  EcalFloatCondObjectContainer record;
  readXML(filename,header,record);
  
  return record.endcapItems();
 
}



std::string 
EcalFloatCondObjectContainerXMLTranslator::dumpXML(       
				  const EcalCondHeader&   header,
				  const std::vector<float>& eb,
				  const std::vector<float>& ee){


  if (eb.size() != EBDetId::kSizeForDenseIndexing){
    std::cerr<<"Error in EcalFloatCondObjectContainerXMLTranslator::dumpXML, invalid Barrel array size: "
	     <<eb.size() << " should be "<<  EBDetId::kSizeForDenseIndexing<< std::endl;
    return std::string("");
  }

  if (ee.size() != EEDetId::kSizeForDenseIndexing){
    std::cerr<<"Error in EcalFloatCondObjectContainerXMLTranslator::dumpXML, invalid Endcap array size: "
	     <<ee.size() << " should be "<<  EEDetId::kSizeForDenseIndexing<< std::endl;
    return std::string("");
  }

  EcalFloatCondObjectContainer record;
  
  for (int cellid = 0; 
       cellid < EBDetId::kSizeForDenseIndexing; 
       ++cellid){// loop on EB cells
        
    uint32_t rawid = EBDetId::unhashIndex(cellid);
    record[rawid]   = eb[cellid];
  } 
  
  for (int cellid = 0; 
       cellid < EEDetId::kSizeForDenseIndexing; 
       ++cellid){// loop on EE cells
    
    
    
    if (EEDetId::validHashIndex(cellid)){  
      uint32_t rawid = EEDetId::unhashIndex(cellid);
   
      record[rawid]=ee[cellid];
    } // if
  }
  
  return dumpXML(header,record);


}

std::string 
EcalFloatCondObjectContainerXMLTranslator::dumpXML(       
				  const EcalCondHeader&   header,
				  const EcalFloatCondObjectContainer& record){

  cms::concurrency::xercesInitialize();

  unique_ptr<DOMImplementation> impl( DOMImplementationRegistry::getDOMImplementation(cms::xerces::uStr("LS").ptr()));
  
  DOMLSSerializer* writer = impl->createLSSerializer();
  if( writer->getDomConfig()->canSetParameter( XMLUni::fgDOMWRTFormatPrettyPrint, true ))
    writer->getDomConfig()->setParameter( XMLUni::fgDOMWRTFormatPrettyPrint, true );
  
  DOMDocumentType* doctype = impl->createDocumentType( cms::xerces::uStr("XML").ptr(), 0, 0 );
  DOMDocument* doc = impl->createDocument( 0, cms::xerces::uStr(EcalFloatCondObjectContainer_tag.c_str()).ptr(), doctype );
  DOMElement* root = doc->getDocumentElement();

  xuti::writeHeader(root, header);

  for (int cellid = EBDetId::MIN_HASH; 
       cellid < EBDetId::kSizeForDenseIndexing; 
       ++cellid){// loop on EB cells
    
    uint32_t rawid= EBDetId::unhashIndex(cellid);
    EcalFloatCondObjectContainer::const_iterator value_ptr=
      record.find(rawid);
    if (value_ptr==record.end()) continue; // cell absent from original record
    
    DOMElement* cellnode=writeCell(root,rawid);

    WriteNodeWithValue(cellnode,Value_tag,*value_ptr);
 

  } // loop on EB cells
  
  
  
  for (int cellid = 0; 
       cellid < EEDetId::kSizeForDenseIndexing; 
       ++cellid){// loop on EE cells
    
    if (!EEDetId::validHashIndex(cellid)) continue;

    uint32_t rawid= EEDetId::unhashIndex(cellid);
    EcalFloatCondObjectContainer::const_iterator value_ptr=
      record.find(rawid);
    if (value_ptr==record.end()) continue; // cell absent from original record
    

    DOMElement* cellnode= writeCell(root,rawid);
    WriteNodeWithValue(cellnode,Value_tag,*value_ptr);

    
  } // loop on EE cells
  
  
  std::string dump = cms::xerces::toString( writer->writeToString( root ));
  doc->release();
  doctype->release();
  writer->release();
  
  cms::concurrency::xercesTerminate();

  return dump;
}



int 
EcalFloatCondObjectContainerXMLTranslator::writeXML(const std::string& filename,         
					       const EcalCondHeader&   header,
					       const EcalFloatCondObjectContainer& record){

  cms::concurrency::xercesInitialize();

  std::fstream fs(filename.c_str(),ios::out);
  fs<< dumpXML(header,record);

  cms::concurrency::xercesTerminate();

  return 0;  
}
