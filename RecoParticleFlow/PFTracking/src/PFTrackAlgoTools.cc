#include "RecoParticleFlow/PFTracking/interface/PFTrackAlgoTools.h"

namespace PFTrackAlgoTools {

  double dPtCut(const reco::TrackBase::TrackAlgorithm& algo,const std::vector<double>& cuts,bool hltIterativeTracking = true){
    switch (algo) {
    case reco::TrackBase::ctf:
    case reco::TrackBase::duplicateMerge:
    case reco::TrackBase::initialStep:
    case reco::TrackBase::highPtTripletStep:
    case reco::TrackBase::lowPtQuadStep:
    case reco::TrackBase::lowPtTripletStep:
    case reco::TrackBase::pixelPairStep:
    case reco::TrackBase::jetCoreRegionalStep:
      return cuts[0];
    case reco::TrackBase::detachedQuadStep:
    case reco::TrackBase::detachedTripletStep:
      return cuts[1];
    case reco::TrackBase::mixedTripletStep:
      return cuts[2];
    case reco::TrackBase::pixelLessStep:
      return cuts[3];
    case reco::TrackBase::tobTecStep:
      return cuts[4];
    case reco::TrackBase::muonSeededStepInOut:
    case reco::TrackBase::muonSeededStepOutIn:
      return cuts[5];
    case reco::TrackBase::hltIter0:
    case reco::TrackBase::hltIter1:
    case reco::TrackBase::hltIter2:
      return cuts[0];
    case reco::TrackBase::hltIter3:
      return  hltIterativeTracking ? cuts[1] : cuts[0];
    case reco::TrackBase::hltIter4:
      return  hltIterativeTracking ? cuts[2] : cuts[0];
    case reco::TrackBase::hltIterX:
      return  cuts[0];
    default:
      return hltIterativeTracking ? cuts[6]:cuts[0];

    }
  }



  unsigned int  nHitCut(const reco::TrackBase::TrackAlgorithm& algo,const std::vector<unsigned int>& cuts,bool hltIterativeTracking = true){
    switch (algo) {
    case reco::TrackBase::ctf:
    case reco::TrackBase::duplicateMerge:
    case reco::TrackBase::initialStep:
    case reco::TrackBase::highPtTripletStep:
    case reco::TrackBase::lowPtQuadStep:
    case reco::TrackBase::lowPtTripletStep:
    case reco::TrackBase::pixelPairStep:
    case reco::TrackBase::jetCoreRegionalStep:
      return cuts[0];
    case reco::TrackBase::detachedQuadStep:
    case reco::TrackBase::detachedTripletStep:
      return cuts[1];
    case reco::TrackBase::mixedTripletStep:
      return cuts[2];
    case reco::TrackBase::pixelLessStep:
      return cuts[3];
    case reco::TrackBase::tobTecStep:
      return cuts[4];
    case reco::TrackBase::muonSeededStepInOut:
    case reco::TrackBase::muonSeededStepOutIn:
      return cuts[5];
    case reco::TrackBase::hltIter0:
    case reco::TrackBase::hltIter1:
    case reco::TrackBase::hltIter2:
      return cuts[0];
    case reco::TrackBase::hltIter3:
      return  hltIterativeTracking ? cuts[1] : cuts[0];
    case reco::TrackBase::hltIter4:
      return  hltIterativeTracking ? cuts[2] : cuts[0];
    case reco::TrackBase::hltIterX:
      return  cuts[0];
    default:
      return hltIterativeTracking ? cuts[6]:cuts[0];

    }
  }



  double  errorScale(const reco::TrackBase::TrackAlgorithm& algo,const std::vector<double>& errorScale){
    switch (algo) {
    case reco::TrackBase::ctf:
    case reco::TrackBase::duplicateMerge:
    case reco::TrackBase::initialStep:
    case reco::TrackBase::highPtTripletStep:
    case reco::TrackBase::lowPtQuadStep:
    case reco::TrackBase::lowPtTripletStep:
    case reco::TrackBase::pixelPairStep:
    case reco::TrackBase::jetCoreRegionalStep:
    case reco::TrackBase::muonSeededStepInOut:
    case reco::TrackBase::muonSeededStepOutIn:
    case reco::TrackBase::detachedQuadStep:
    case reco::TrackBase::detachedTripletStep:
    case reco::TrackBase::mixedTripletStep:
    case reco::TrackBase::hltIter0:
    case reco::TrackBase::hltIter1:
    case reco::TrackBase::hltIter2:
    case reco::TrackBase::hltIter3:
    case reco::TrackBase::hltIter4:
    case reco::TrackBase::hltIterX:
      return 1.0;
    case reco::TrackBase::pixelLessStep:
      return errorScale[0];
    case reco::TrackBase::tobTecStep:
      return errorScale[1];
    default:
      return 1E9;
    }
  }


bool isGoodForEGM(const reco::TrackBase::TrackAlgorithm& algo){


  switch (algo) {
  case reco::TrackBase::ctf:
  case reco::TrackBase::duplicateMerge:
  case reco::TrackBase::initialStep:
  case reco::TrackBase::highPtTripletStep:
  case reco::TrackBase::lowPtQuadStep:
  case reco::TrackBase::lowPtTripletStep:
  case reco::TrackBase::pixelPairStep:
  case reco::TrackBase::jetCoreRegionalStep:
  case reco::TrackBase::detachedQuadStep:
  case reco::TrackBase::detachedTripletStep:
  case reco::TrackBase::mixedTripletStep:
  case reco::TrackBase::muonSeededStepInOut:
  case reco::TrackBase::muonSeededStepOutIn:
  case reco::TrackBase::hltIter0:
  case reco::TrackBase::hltIter1:
  case reco::TrackBase::hltIter2:
  case reco::TrackBase::hltIter3:
  case reco::TrackBase::hltIter4:
  case reco::TrackBase::hltIterX:
    return true;
  default:
    return false;
  }

}

bool isGoodForEGMPrimary(const reco::TrackBase::TrackAlgorithm& algo){
  switch (algo) {
  case reco::TrackBase::ctf:
  case reco::TrackBase::duplicateMerge:
  case reco::TrackBase::cosmics:
  case reco::TrackBase::initialStep:
  case reco::TrackBase::highPtTripletStep:
  case reco::TrackBase::lowPtQuadStep:
  case reco::TrackBase::lowPtTripletStep:
  case reco::TrackBase::pixelPairStep:
  case reco::TrackBase::detachedQuadStep:
  case reco::TrackBase::detachedTripletStep:
  case reco::TrackBase::mixedTripletStep:
    return true;
  default:
    return false;
  }

}

bool isFifthStep(const reco::TrackBase::TrackAlgorithm& algo){
  switch (algo) {
  case reco::TrackBase::ctf:
  case reco::TrackBase::duplicateMerge:
  case reco::TrackBase::initialStep:
  case reco::TrackBase::highPtTripletStep:
  case reco::TrackBase::lowPtQuadStep:
  case reco::TrackBase::lowPtTripletStep:
  case reco::TrackBase::pixelPairStep:
  case reco::TrackBase::jetCoreRegionalStep:
  case reco::TrackBase::detachedQuadStep:
  case reco::TrackBase::detachedTripletStep:
  case reco::TrackBase::mixedTripletStep:
  case reco::TrackBase::pixelLessStep:
  case reco::TrackBase::muonSeededStepInOut:
  case reco::TrackBase::muonSeededStepOutIn:
  case reco::TrackBase::hltIter0:
  case reco::TrackBase::hltIter1:
  case reco::TrackBase::hltIter2:
  case reco::TrackBase::hltIter3:
  case reco::TrackBase::hltIter4:
  case reco::TrackBase::hltIterX:
    return false;
  case reco::TrackBase::tobTecStep:
    return true;
  default:
    return true;
  }

}


bool highQuality(const reco::TrackBase::TrackAlgorithm& algo){
  switch (algo) {
  case reco::TrackBase::initialStep:
  case reco::TrackBase::highPtTripletStep:
  case reco::TrackBase::lowPtQuadStep:
  case reco::TrackBase::lowPtTripletStep:
  case reco::TrackBase::pixelPairStep:
  case reco::TrackBase::detachedQuadStep:
  case reco::TrackBase::detachedTripletStep:
  case reco::TrackBase::duplicateMerge:
  case reco::TrackBase::jetCoreRegionalStep:
    return true;
  default:
    return false;
  }

}


bool nonIterative(const reco::TrackBase::TrackAlgorithm& algo){
  switch (algo) {
  case reco::TrackBase::undefAlgorithm:
  case reco::TrackBase::ctf:
  case reco::TrackBase::cosmics:
    return true;
  default:
    return false;

  }

}


bool step45(const reco::TrackBase::TrackAlgorithm& algo){
  switch (algo) {
  case reco::TrackBase::mixedTripletStep:
  case reco::TrackBase::pixelLessStep:
  case reco::TrackBase::tobTecStep:
    return true;
  default:
    return false;
  }

}


bool step5(const reco::TrackBase::TrackAlgorithm& algo){
  return (algo==reco::TrackBase::tobTecStep||algo==reco::TrackBase::pixelLessStep);
}

bool goodPtResolution(const reco::TrackRef& trackref, 
		      const std::vector<double>& DPtovPtCut,
		      const std::vector<unsigned>& NHitCut,
		      bool useIterTracking,
		      bool debug) {
  //recheck that the track is high purity!
  if (!trackref->quality(reco::TrackBase::highPurity))
    return false;
    
  const double P = trackref->p();
  const double Pt = trackref->pt();
  const double DPt = trackref->ptError();
  const unsigned int NHit = 
    trackref->hitPattern().trackerLayersWithMeasurement();
  const unsigned int NLostHit = 
    trackref->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  const unsigned int LostHits = trackref->numberOfLostHits();
  const double sigmaHad = std::sqrt(1.20*1.20/P+0.06*0.06) / (1.+LostHits);

  // Protection against 0 momentum tracks
  if ( P < 0.05 ) return false;
 
  if (debug) std::cout << " PFBlockAlgo: PFrecTrack->Track Pt= "
		       << Pt << " DPt = " << DPt << std::endl;


  double dptCut = dPtCut(trackref->algo(),DPtovPtCut,useIterTracking);
  unsigned int nhitCut = nHitCut(trackref->algo(),NHitCut,useIterTracking);

  if ( ( dptCut > 0. && 
	 DPt/Pt > dptCut*sigmaHad ) || 
       NHit < nhitCut ) { 
    if (debug) std::cout << " PFBlockAlgo: skip badly measured track"
			 << ", P = " << P 
			 << ", Pt = " << Pt 
			 << " DPt = " << DPt 
			 << ", N(hits) = " << NHit 
			 << " (Lost : " << LostHits << "/" << NLostHit << ")"
			 << ", Algo = " << trackref->algo()
			 << std::endl;
    if (debug) std::cout << " cut is DPt/Pt < " << dptCut * sigmaHad << std::endl;
    if (debug) std::cout << " cut is NHit >= " << nhitCut << std::endl;
    return false;
  }

  return true;
}

}
