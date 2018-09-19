/*!
 *  \file       PHG4EICTrackingBase.h
 *  \brief      Kalman Filter based on smeared truth PHG4Hit
 *  \details    Kalman Filter based on smeared truth PHG4Hit
 *  \author     Carlos Perez <carlos.perezlara@stonybrook.edu>
 */

#ifndef __PHG4EICTrackingBase_H__
#define __PHG4EICTrackingBase_H__

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>
#include <phgenfit/Measurement.h>

#include <TMatrixDSym.h>
#include <TVector3.h>

class TFile;
class TH1F;
class TH2F;
class TH3F;
class PHG4Particle;
namespace PHGenFit {
  class PlanarMeasurement;
  class Track;
  class Fitter;
}
namespace genfit {
  class GFRaveVertexFactory;
}
class SvtxTrack;

class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;
class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxClusterMap;
class SvtxEvalStack;
class TFile;
class TTree;

class PHG4EICTrackingBase: public SubsysReco {
 public:
  enum DETECTOR_TYPE {Vertical_Plane, Cylinder};
  PHG4EICTrackingBase(const std::string &name = "PHG4EICTrackingBaseEIC",
		   const std::string &file = "none");
  ~PHG4EICTrackingBase();
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  int End(PHCompositeNode *);
  void Verbosity(int verb) {verbosity = verb;}
  void Display() { fDisplay = true; }
  void Configure(const std::string* phg4hitsNames,
		 const DETECTOR_TYPE* phg4dettype, 
		 const float *radres,
		 const float *phires,
		 const float *lonres,
		 const float *eff,
		 const float *noise,
		 const int nlayer) {
    fHitsNames.clear();
    fDetectorType.clear();
    fDetectorRadRes.clear();
    fDetectorPhiRes.clear();
    fDetectorLonRes.clear();
    fDetectorHitFinEff.clear();
    fDetectorHitNoise.clear();
    for(int i=0;i!=nlayer;++i) {
      fHitsNames.push_back(phg4hitsNames[i]);
      fDetectorType.push_back(phg4dettype[i]);
      fDetectorRadRes.push_back(radres[i]);
      fDetectorPhiRes.push_back(phires[i]);
      fDetectorLonRes.push_back(lonres[i]);
      fDetectorHitFinEff.push_back(eff[i]);
      fDetectorHitNoise.push_back(noise[i]);
    }
  }
  void VertexSmear(float dxy, float dz) {
    fVertexSmearXY = dxy;
    fVertexSmearZ = dz;
  }
  void VertexIn() { fVertexIn = true; }
  void TrackSecondaries() { fPrimaryOnly = false; }

private:
  int CreateNodes(PHCompositeNode *);
  int GetNodes(PHCompositeNode *);
  void FindMeasurements(const PHG4Particle* particle,
			TVector3& vtx,
			std::map<double,PHGenFit::Measurement*>& meas_out);
  //std::vector<PHGenFit::Measurement*>& meas_out);
  PHGenFit::PlanarMeasurement* MeasurementVertical(const PHG4Hit *g4hit,
						   const double phires,
						   const double radres);
  PHGenFit::PlanarMeasurement* MeasurementOnion(const PHG4Hit *g4hit,
						const double phires,
						const double lonres);
  SvtxTrack* MakeSvtxTrack(const PHGenFit::Track *phgf_track_in, 
			   const PHG4Particle* particle,
			   const unsigned int nmeas,
			   const TVector3& vtx);

  PHG4TruthInfoContainer* fTruths;
  std::vector<PHG4HitContainer*> fHits;
  std::vector<std::string> fHitsNames;
  std::vector<DETECTOR_TYPE> fDetectorType;
  std::vector<float> fDetectorRadRes;
  std::vector<float> fDetectorPhiRes;
  std::vector<float> fDetectorLonRes;
  std::vector<float> fDetectorHitFinEff;
  std::vector<float> fDetectorHitNoise;
  SvtxTrackMap* fTrackMap;
  PHGenFit::Fitter* fMyFitter;
  std::string fFittingMethod;
  bool fDisplay;
  bool fVertexIn;
  double fVertexSmearXY;
  double fVertexSmearZ;
  bool fPrimaryOnly;

  std::string fFileOut;
  std::vector<TH1F*> fHistHitsEta;
  std::vector<TH1F*> fHistHitsPhi;
  std::vector<TH1F*> fHistHitsRad;
  std::vector<TH2F*> fHistHitsSme;
  std::vector<TH2F*> fHistHitsTru;

  TH2F *fHistCounter;
  TH1F *fHistVertexReco;
  TH1F *fHistParticlesEta;
  TH1F *fHistParticlesPhi;
  TH1F *fHistRecoPartEta;
  TH1F *fHistRecoPartPhi;
  TH2F *fHistTracksEta;
  TH2F *fHistTracksPhi;
  TH2F *fHistFailsEta;
  TH2F *fHistMomMagEta;
  TH2F *fHistMomMagPhi;
  TH2F *fHistMomEtaEta;
  TH2F *fHistMomPhiPhi;
  TH3F *fHistMomMagMap;
  TH3F *fHistMomEtaMap;
  TH3F *fHistMomPhiMap;
  TH2F *fHistParticleMap;
  TH2F *fHistTrackMap;
};

#endif /*__PHG4EICTrackingBase_H__*/
