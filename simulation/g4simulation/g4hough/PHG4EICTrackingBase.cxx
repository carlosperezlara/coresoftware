/*!
 *  \file       PHG4EICTrackingBase.cxx
 *  \brief      Quick MC-based Reconstruction Method
 *  \details    Quick MC-based Reconstruction Method
 *  \author     Carlos Perez <carlos.perezlara@stonybrook.edu>
 */
#include <cmath>
#include <map>
#include <utility>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <TMath.h>
#include <TMatrixF.h>
#include <TRandom.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TVector3.h>
#include <TRotation.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4VtxPointv1.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <phgeom/PHGeomUtility.h>
#include <phfield/PHFieldUtility.h>

#include <GenFit/AbsMeasurement.h>
#include <GenFit/EventDisplay.h>
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/StateOnPlane.h>
#include <GenFit/Track.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/Track.h>
#include <phgenfit/SpacepointMeasurement.h>

#include "SvtxTrackMap.h"
#include "SvtxTrackMap_v1.h"
#include "SvtxTrackState.h"
#include "SvtxTrackState_v1.h"
#include "SvtxTrack.h"
#include "SvtxTrack_FastSim.h"

#include "PHG4EICTrackingBase.h"

using namespace std;
//========================================
PHG4EICTrackingBase::PHG4EICTrackingBase(const std::string &name, const std::string &fileout) :
  SubsysReco(name),
  fTruths(NULL),
  fTrackMap(NULL),
  fMyFitter(NULL),
  fFittingMethod("KalmanFitterRefTrack"),
  //KalmanFitter, KalmanFitterRefTrack, DafSimple, DafRef 
  fDisplay(false),
  fVertexIn(true),
  fVertexSmearXY(50E-4),
  fVertexSmearZ(50E-4),
  fPrimaryOnly(true),
  fFileOut(fileout),
  fHistCounter(NULL),
  fHistVertexReco(NULL),
  fHistParticlesEta(NULL),
  fHistParticlesPhi(NULL),
  fHistRecoPartEta(NULL),
  fHistRecoPartPhi(NULL),
  fHistTracksEta(NULL),
  fHistTracksPhi(NULL),
  fHistFailsEta(NULL),
  fHistMomMagEta(NULL),
  fHistMomMagPhi(NULL),
  fHistMomEtaEta(NULL),
  fHistMomPhiPhi(NULL),
  fHistMomMagMap(NULL),
  fHistMomEtaMap(NULL),
  fHistMomPhiMap(NULL),
  fHistParticleMap(NULL),
  fHistTrackMap(NULL)
{
  fHistHitsEta.clear();
  fHistHitsPhi.clear();
  fHistHitsRad.clear();
}
//========================================
PHG4EICTrackingBase::~PHG4EICTrackingBase() {
  delete fMyFitter;
}
//========================================
int PHG4EICTrackingBase::Init(PHCompositeNode *topNode) {
  fHistCounter = new TH2F("HistCounter","Counter;Particles Per Event;Tracks Per Event",
			  10,-0.5,9.5, 10,-0.5,9.5);
  fHistVertexReco = new TH1F("HistVertexReco","Vertex; Distance3D( Reco - Truth )",100,0,0.1);
  fHistParticlesEta = new TH1F("HistParticlesEta","Particles Eta;ETA",100,-4.5,+4.5);
  fHistParticlesPhi = new TH1F("HistParticlesPhi","Particles Phi;PHI",100,-3.2,+3.2);
  fHistRecoPartEta = new TH1F("HistRecoPartEta","Reco Particles Eta;ETA",100,-4.5,+4.5);
  fHistRecoPartPhi = new TH1F("HistRecoPartPhi","Reco Particles Phi;PHI",100,-3.2,+3.2);
  fHistTracksEta = new TH2F("HistTracksEta","Tracks Eta;ETA",100,-4.5,+4.5,70,-0.5,69.5);
  fHistTracksPhi = new TH2F("HistTracksPhi","Tracks Phi;PHI",100,-3.2,+3.2,70,-0.5,69.5);
  fHistFailsEta=new TH2F("HistFailsEta","Particle Fails Eta;ETA",100,-4.5,+4.5,70,-0.5,69.5);
  fHistMomMagEta = new TH2F("HistMomMagEta","Momentum;ETA;dp / p",
					100,-4.5,+4.5,100,-0.5,+0.5);
  fHistMomMagPhi = new TH2F("HistMomMagPhi","Momentum;PHI;dp / p",
					100,-3.2,+3.2,100,-0.5,+0.5);
  fHistMomEtaEta = new TH2F("HistMomEtaEta","Momentum  ETA;TRUTH;RECO",
					100,-4.5,+4.5,100,-4.5,+4.5);
  fHistMomPhiPhi = new TH2F("HistMomPhiPhi","Momentum  PHI;TRUTH;RECO",
					100,-3.2,+3.2,100,-3.2,+3.2);
  fHistMomMagMap = new TH3F("HistMomMagMap","MomMag;P;Eta",
			    40,0.,40.,18,-4.5,+4.5, 60, -0.5,+0.5);
  fHistMomPhiMap = new TH3F("HistMomPhiMap","MomPhi;P;Eta",
			    40,0.,40.,18,-4.5,+4.5, 60, -0.5,+0.5);
  fHistMomEtaMap = new TH3F("HistMomEtaMap","MomEta;P;Eta",
			    40,0.,40.,18,-4.5,+4.5, 60, -0.5,+0.5);
  fHistParticleMap = new TH2F("HistParticleMap","Particle;P;Eta", 40,0.,40.,18,-4.5,+4.5);
  fHistTrackMap = new TH2F("HistTrackMap","Track;P;Eta", 40,0.,40.,18,-4.5,+4.5);
  float mrad[16] = {5,85, 20,20,80,25,90,30, 20,25,90,30,90,30,160,55};
  for(unsigned int i=0; i!=fHitsNames.size(); ++i) {
    fHistHitsEta.push_back( new TH1F( Form("HistHitsEta_%d",i),
				      Form("Det_%s;ETA",fHitsNames[i].c_str()),
				      100,-4.5,+4.5 ) );
    fHistHitsPhi.push_back( new TH1F( Form("HistHitsPhi_%d",i),
				      Form("Det_%s;PHI",fHitsNames[i].c_str()),
				      100,-3.2,+3.2 ) ); 
    fHistHitsRad.push_back( new TH1F( Form("HistHitsRad_%d",i),
				      Form("Det_%s;RAD",fHitsNames[i].c_str()),
				      100,0,mrad[i] ) );
    if(fDetectorType[i] == Cylinder) {
      fHistHitsSme.push_back( new TH2F( Form("HistHitsSme_%d",i),
					Form("Det_%s;RPhi;Z",fHitsNames[i].c_str()),
					50,-5*fDetectorPhiRes[i],+5*fDetectorPhiRes[i],
					50,-5*fDetectorLonRes[i],+5*fDetectorLonRes[i] ) );
      fHistHitsTru.push_back( new TH2F( Form("HistHitsTru_%d",i),
					Form("Det_%s;Phi;Rad",fHitsNames[i].c_str()),
					50,-3.2,+3.2,
					50,0,mrad[i] ) );
    }
    if(fDetectorType[i] == Vertical_Plane) {
      fHistHitsSme.push_back( new TH2F( Form("HistHitsSme_%d",i),
					Form("Det_%s;RPhi;Rad",fHitsNames[i].c_str()),
					50,-5*fDetectorPhiRes[i],+5*fDetectorPhiRes[i],
					50,-5*fDetectorRadRes[i],+5*fDetectorRadRes[i] ) );
      fHistHitsTru.push_back( new TH2F( Form("HistHitsTru_%d",i),
					Form("Det_%s;X;Y",fHitsNames[i].c_str()),
					50,-mrad[i],+mrad[i],
					50,-mrad[i],+mrad[i] ) );
    }
 }
  return Fun4AllReturnCodes::EVENT_OK;
}
//========================================
int PHG4EICTrackingBase::InitRun(PHCompositeNode *topNode) {
  CreateNodes(topNode);
  TGeoManager *tgeomanager = PHGeomUtility::GetTGeoManager(topNode);
  PHField *field = PHFieldUtility::GetFieldMapNode(nullptr, topNode);
  fMyFitter = PHGenFit::Fitter::getInstance(tgeomanager,
					    field,
					    fFittingMethod,
					    "RKTrackRep",
					    fDisplay);
  if(!fMyFitter) {
    std::cout << PHWHERE << " Cannot get instance of PHGenFit::Fitter" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  fMyFitter->set_verbosity(verbosity);
  GetNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}
//========================================
int PHG4EICTrackingBase::End(PHCompositeNode *topNode) {
  if(fDisplay) {
    fMyFitter->displayEvent();
  }
  if(!fFileOut.compare("none")) return Fun4AllReturnCodes::EVENT_OK;

  TFile *file = new TFile(fFileOut.c_str(),"RECREATE");
  fHistVertexReco->Write("Vertex");
  fHistCounter->Write("Counter");
  fHistParticlesPhi->Write("ParticlesPhi");
  fHistParticlesEta->Write("ParticlesEta");
  fHistRecoPartEta->Write("RecoEta");
  fHistRecoPartPhi->Write("RecoPhi");
  fHistTracksEta->Write("TracksEta");
  fHistTracksPhi->Write("TracksPhi");
  fHistFailsEta->Write("FailsEta");
  fHistMomMagEta->Write("MomMagEta");
  fHistMomMagPhi->Write("MomMagPhi");
  fHistMomEtaEta->Write("MomEtaEta");
  fHistMomPhiPhi->Write("MomPhiPhi");
  fHistMomMagMap->Write("MomMagMap");
  fHistMomPhiMap->Write("MomPhiMap");
  fHistMomEtaMap->Write("MomEtaMap");
  fHistParticleMap->Write("ParticleMap");
  fHistTrackMap->Write("TrackMap");
  for(unsigned int i=0; i!=fHitsNames.size(); ++i) {
    fHistHitsEta[i]->Write(Form("ETADet_%s",fHitsNames[i].c_str()));
    fHistHitsPhi[i]->Write(Form("PHIDet_%s",fHitsNames[i].c_str()));
    fHistHitsRad[i]->Write(Form("RADDet_%s",fHitsNames[i].c_str()));
    fHistHitsSme[i]->Write(Form("SMEDet_%s",fHitsNames[i].c_str()));
    fHistHitsTru[i]->Write(Form("TRUDet_%s",fHitsNames[i].c_str()));
  }
  file->Close();
  std::cout << "PHG4EICTrackingBase => " << fFileOut.c_str() << std::endl;
  delete file;
  return Fun4AllReturnCodes::EVENT_OK;
}
//========================================
int PHG4EICTrackingBase::process_event(PHCompositeNode *topNode) {
  fTrackMap->empty();
  vector<PHGenFit::Track*> rf_tracks;

  //=== RECONSTRUCT VERTEX
  PHG4VtxPoint *MyVertex = fTruths->GetPrimaryVtx(fTruths->GetPrimaryVertexIndex());
  double dvx = gRandom->Gaus(0, fVertexSmearXY);
  double dvy = gRandom->Gaus(0, fVertexSmearXY);
  double dvz = gRandom->Gaus(0, fVertexSmearZ);
  double dist = TMath::Sqrt( dvx*dvx + dvy*dvy + dvz*dvz );
  fHistVertexReco->Fill( dist );
  MyVertex->set_x(MyVertex->get_x()+ dvx);
  MyVertex->set_y(MyVertex->get_y()+ dvy);
  MyVertex->set_z(MyVertex->get_z()+ dvz);
  TVector3 MyRootVertex;
  MyRootVertex.SetX(MyVertex->get_x());
  MyRootVertex.SetY(MyVertex->get_y());
  MyRootVertex.SetZ(MyVertex->get_z());
  
  //=== LOOPING OVER PARTICLES
  PHG4TruthInfoContainer::ConstRange irange;
  if(fPrimaryOnly){
    irange = fTruths->GetPrimaryParticleRange();
  } else{
    irange = fTruths->GetParticleRange();
  }
  int numprim = 0;
  int numtrks = 0;
  for(PHG4TruthInfoContainer::ConstIterator ithis = irange.first;
      ithis != irange.second;
      ++ithis) {
    ++numprim;
    PHG4Particle* particle = ithis->second;
    TVector3 gp( particle->get_px(), particle->get_py(), particle->get_pz() );
    if(gp.Mag()<2.0) {
      std::cout << "Particle " << particle->get_pid() <<  " found with momentum: " << gp.Mag() << " GeV/c. Skipping" << std::endl;
      continue;
    }
    if(TMath::Abs(gp.Eta())>4.5) {
      std::cout << "Particle " << particle->get_pid() <<  " found with eta: " << gp.Eta() << ". Skipping" << std::endl;
      continue;
    }
    std::cout << "Particle " << particle->get_pid();
    std::cout <<  " | Momentum: " << gp.Mag() << " GeV/c";
    std::cout <<  " | Eta: " << gp.Eta() << " Phi: " << gp.Phi() << std::endl;    
    fHistParticlesEta->Fill( gp.Eta() );
    fHistParticlesPhi->Fill( gp.Phi() );
    fHistParticleMap->Fill( gp.Mag(), gp.Eta() );
    std::map<double,PHGenFit::Measurement*> sorted;
    
    if(fVertexIn) { //=== ADD A VERTEX MEASUREMENT
      PHGenFit::Measurement* meas = NULL;
      TMatrixDSym cov(3);
      cov.Zero();
      cov(0, 0) = fVertexSmearXY * fVertexSmearXY;
      cov(1, 1) = fVertexSmearXY * fVertexSmearXY;
      cov(2, 2) = fVertexSmearZ * fVertexSmearZ;
      meas = new PHGenFit::SpacepointMeasurement(MyRootVertex, cov);
      sorted[0.0] = meas;
      //measurements.push_back(meas);
    }
    //=== ADDING COMPATIBLE MEASUREMENTS
    FindMeasurements(particle, MyRootVertex, sorted/*measurements*/);
    std::vector<PHGenFit::Measurement*> measurements;
    //std::cout << " SORTING ";
    for(const auto& kv: sorted) {
      //std::cout << kv.first << " ";
      measurements.push_back( kv.second );
    }
    //std::cout << std::endl;
    int nmeasurements = measurements.size();
    if(nmeasurements < 3) {
      if(verbosity >= 2) {
	std::cout << "  Found measurements.size() < 3; skipping\n";
      }
      fHistFailsEta->Fill( gp.Eta(), nmeasurements );
      // === clean and skip
      for(unsigned int im=0; im<measurements.size(); im++) {
	delete measurements[im]->getMeasurement(); 
	delete measurements[im]; 
      }
      continue;
    }
    if(verbosity >= 20)
      std::cout << "  measurements.size() " << measurements.size() << std::endl;
    //=== SEEDING
    TVector3 seed_pos(0,0,0);//MyVertex->get_x(), MyVertex->get_y(), MyVertex->get_z());
    TVector3 seed_mom(0, 0, 1); //unit vector
    TVector3 True_mom(particle->get_px(), particle->get_py(), particle->get_pz());
    const double onedeg = 3*TMath::Pi()/180.0; // 1deg
    const double mag5percent = 0.10*True_mom.Mag();
    seed_mom.SetPhi(   gRandom->Gaus(True_mom.Phi(),   onedeg));
    seed_mom.SetTheta( gRandom->Gaus(True_mom.Theta(), onedeg));
    seed_mom.SetMag(   gRandom->Gaus(True_mom.Mag(),   mag5percent));
    TMatrixDSym seed_cov(6);
    seed_cov.ResizeTo(6, 6);
    if(0) {
      double covXY = fVertexSmearXY*fVertexSmearXY;
      double covZ = fVertexSmearZ*fVertexSmearZ;
      double covPX = 2*(True_mom.X() - seed_mom.X())*(True_mom.X() - seed_mom.X());
      double covPY = 2*(True_mom.Y() - seed_mom.Y())*(True_mom.Y() - seed_mom.Y());
      double covPZ = 2*(True_mom.Z() - seed_mom.Z())*(True_mom.Z() - seed_mom.Z());
      seed_cov[0][0] = covXY;
      seed_cov[1][1] = covXY;
      seed_cov[2][2] = covZ;
      seed_cov[3][3] = covPX;
      seed_cov[4][4] = covPY;
      seed_cov[5][5] = covPZ;
    } else {
      for(int i=0; i!=3; ++i)
	seed_cov[i][i] = 1;//2.5E-6;
      for(int i=3; i!=6; ++i)
	seed_cov[i][i] = 10;
    }
    //=== FITTING
    genfit::AbsTrackRep *rep = new genfit::RKTrackRep(particle->get_pid());
    PHGenFit::Track* track = new PHGenFit::Track(rep, seed_pos, seed_mom, seed_cov);
    rf_tracks.push_back(track);
    track->addMeasurements(measurements);
    int fitting_err=-1;
    int ntrials = 5;
    //std::cout << " **** FITING PART ETA " << True_mom.Eta() << " PHI " << True_mom.Phi() << std::endl;
    while(ntrials--) {
      fitting_err = fMyFitter->processTrack(track, false);
      if(fitting_err==0) break;
      if(verbosity>=0) {
	std::cout << "Failed with error " << fitting_err;
	std::cout << "... tying again" << ntrials << endl;
      }
    }
    if(fitting_err != 0) {
      std::cout << "Yiaks! Give up after several tries. Fitting_err = " << fitting_err << ". Skipping..." << "\n";
      fHistFailsEta->Fill( gp.Eta(), nmeasurements );
      continue; 
    }
    if(verbosity >= 20)
      std::cout << "  DONE FITTING fitting_err " << fitting_err << std::endl;

    //=== BUILDING TRACK
    TVector3 vtx( MyVertex->get_x(), MyVertex->get_y(), MyVertex->get_z()  );
    SvtxTrack *MyTrack = MakeSvtxTrack(track, particle,
				       measurements.size(), vtx);
    if(MyTrack) {
      TVector3 gt( MyTrack->get_px(), MyTrack->get_py(), MyTrack->get_pz() );
      fHistRecoPartEta->Fill( gp.Eta() );
      fHistRecoPartPhi->Fill( gp.Phi() );
      fHistTracksEta->Fill( gt.Eta(), nmeasurements );
      fHistTracksPhi->Fill( gt.Phi(), nmeasurements );
      fHistMomMagEta->Fill( gp.Eta(), gt.Mag()/gp.Mag()-1 );
      fHistMomMagPhi->Fill( gp.Phi(), gt.Mag()/gp.Mag()-1 );
      fHistMomEtaEta->Fill( gp.Eta(), gt.Eta() );
      fHistMomPhiPhi->Fill( gp.Phi(), gt.Phi() );
      fHistMomMagMap->Fill( gp.Mag(), gp.Eta(), gt.Mag()/gp.Mag()-1 );
      fHistMomPhiMap->Fill( gp.Mag(), gp.Eta(), gt.Phi()/gp.Phi()-1 );
      fHistMomEtaMap->Fill( gp.Mag(), gp.Eta(), gt.Eta()/gp.Eta()-1 );
      fHistTrackMap->Fill( gp.Mag(), gp.Eta() );
      fTrackMap->insert(MyTrack);
      delete MyTrack;
      ++numtrks;
    }
    if(verbosity >= 20)
      std::cout << "  ALL DONE" << std::endl;
  } /*ithis*/

  if(fDisplay){
    vector<genfit::Track*> rf_gf_tracks;
    for(std::vector<PHGenFit::Track*>::iterator it = rf_tracks.begin();
  	it != rf_tracks.end(); ++it) 
      rf_gf_tracks.push_back((*it)->getGenFitTrack()); 	  
    fMyFitter->getEventDisplay()->addEvent(rf_gf_tracks);
  } else{
    for (std::vector<PHGenFit::Track*>::iterator it = rf_tracks.begin();
  	 it != rf_tracks.end(); ++it)
      delete (*it);
    rf_tracks.clear();
  }

  fHistCounter->Fill(float(numprim),float(numtrks));
  return Fun4AllReturnCodes::EVENT_OK;
}
//========================================
int PHG4EICTrackingBase::CreateNodes(PHCompositeNode *topNode) {
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cerr << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHNodeIterator iter_dst(dstNode);
  PHCompositeNode *tb_node = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode",
									       "Svtx"));
  if(!tb_node) {
    tb_node = new PHCompositeNode("Svtx");
    dstNode->addNode(tb_node);
    if(verbosity > 0)
      cout << "Svtx node added" << endl;
  }
  fTrackMap = new SvtxTrackMap_v1;
  PHIODataNode<PHObject>* tracks_node = new PHIODataNode<PHObject>(fTrackMap,
								   "SvtxTrackMap",
								   "PHObject");
  tb_node->addNode(tracks_node);
  if (verbosity > 0)
    cout << "SvtxTrackMap node added" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
//========================================
int PHG4EICTrackingBase::GetNodes(PHCompositeNode *topNode) {
  fTruths = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if(!fTruths) {
    cout << PHWHERE << "G4TruthInfo node not found on node tree" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  for(unsigned int i=0; i!=fHitsNames.size(); ++i) {
    PHG4HitContainer *phg4hit = findNode::getClass<PHG4HitContainer>(topNode, fHitsNames[i].c_str());
    if(!phg4hit) {
      cout << PHWHERE << fHitsNames[i].c_str()
	   << " node not found on node tree" << endl;
      //return Fun4AllReturnCodes::ABORTRUN;
      fHits.push_back(NULL);
    } else {
      fHits.push_back(phg4hit);
    }
  }
  fTrackMap = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
  if(!fTrackMap) {
    cout << PHWHERE << "SvtxTrackMap node not found on node tree" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
//========================================
void PHG4EICTrackingBase::FindMeasurements(const PHG4Particle* particle,
					TVector3& vtx,
					std::map<double,PHGenFit::Measurement*>& meas_out) {
  //std::vector<PHGenFit::Measurement*>& meas_out/*,
  //std::vector<TVector3>& foundorder*/) {
  if(!particle) return;
  for(unsigned int idet = 0; idet != fHitsNames.size(); idet++) {
    if(!fHits[idet]) {
      std::cout << "No fHits["<< idet <<"] found!" << std::endl;
      continue;
    }    
    int dettype = fDetectorType[idet];
    float detradres = fDetectorRadRes[idet];
    float detphires = fDetectorPhiRes[idet];
    float detlonres = fDetectorLonRes[idet];
    float dethiteff = fDetectorHitFinEff[idet];
    float detnoise = fDetectorHitNoise[idet];
    if(verbosity>=20) {
      std::cout << "   fHits["<< idet <<"] ==> " << fHitsNames[idet].c_str() << std::endl;
      std::cout << "    dettype " << dettype << std::endl;
      std::cout << "    detphires " << detphires << std::endl;
      std::cout << "    detlonres " << detlonres << std::endl;
      std::cout << "    detradres " << detradres << std::endl;
      std::cout << "    detnoise " << detnoise << std::endl;
      std::cout << "    dethiteff " << dethiteff << std::endl;
    }
    int nhits = 0;
    int nhitsasoc = 0;
    for(PHG4HitContainer::LayerIter ilayer = fHits[idet]->getLayers().first;
	ilayer != fHits[idet]->getLayers().second;
	++ilayer) {
      for(PHG4HitContainer::ConstIterator ihit = fHits[idet]->getHits(*ilayer).first;
	  ihit != fHits[idet]->getHits(*ilayer).second; ++ihit) {
	PHG4Hit *hit = ihit->second;
	if(!hit) {
	  std::cout << "   No PHG4Hit Found!" << std::endl;
	  continue;
	}
	TVector3 hvec(hit->get_avg_x(),hit->get_avg_y(),hit->get_avg_z());
	double distance = TMath::Sqrt( (hvec.X() - vtx.X())*(hvec.X() - vtx.X()) +
				       (hvec.Y() - vtx.Y())*(hvec.Y() - vtx.Y()) +
				       (hvec.Z() - vtx.Z())*(hvec.Z() - vtx.Z()) );
	double hvecrad = TMath::Sqrt( hvec.X()*hvec.X() +
				      hvec.Y()*hvec.Y() );
	fHistHitsRad[idet]->Fill(hvecrad);
	fHistHitsEta[idet]->Fill(hvec.Eta());
	fHistHitsPhi[idet]->Fill(hvec.Phi());
	++nhits;
	if(hit->get_trkid() != particle->get_track_id() &&
	   gRandom->Uniform(0, 1) < detnoise) continue;
	if(gRandom->Uniform(0, 1) > dethiteff) continue;
	//adding measurement
	++nhitsasoc;
	PHGenFit::Measurement *meas = NULL;
	if(dettype == Vertical_Plane){
	  TVector3 pos(hit->get_avg_x(), hit->get_avg_y(), hit->get_avg_z());
	  TVector3 v(pos.X(), pos.Y(), 0);
	  v = 1 / v.Mag() * v; //========== r unit
	  TVector3 u = -v.Cross(TVector3(0, 0, 1));
	  u = 1 / u.Mag() * u; //== phi unit
	  if(true) {
	    double u_smear = gRandom->Gaus(0, detphires);
	    double v_smear = gRandom->Gaus(0, detradres);
	    pos.SetX(hit->get_avg_x() + u_smear * u.X() + v_smear * v.X());
	    pos.SetY(hit->get_avg_y() + u_smear * u.Y() + v_smear * v.Y());
	    if(verbosity >= 40) {
	      std::cout << "{ " << pos.X() << ", " << pos.Y();
	      std::cout << ", " << pos.Z() << " }, ";
	    }
	    fHistHitsSme[idet]->Fill(u_smear,v_smear);
	  }
	  fHistHitsTru[idet]->Fill(hvec.X(),hvec.Y());
	  meas = new PHGenFit::PlanarMeasurement(pos, u, v, detphires, detradres);
	  //foundorder.push_back(pos);
	} else if (dettype == Cylinder){
	  TVector3 pos(hit->get_avg_x(), hit->get_avg_y(), hit->get_avg_z());
	  TVector3 v(0, 0, 1); //=========================== z unit
	  TVector3 u = v.Cross(TVector3(pos.X(), pos.Y(), 0));
	  u = 1 / u.Mag() * u; //= phi unit
	  if(true) {
	    double u_smear = gRandom->Gaus(0, detphires);
	    double v_smear = gRandom->Gaus(0, detlonres);
	    pos.SetX(hit->get_avg_x() + u_smear * u.X() + v_smear * v.X());
	    pos.SetY(hit->get_avg_y() + u_smear * u.Y() + v_smear * v.Y());
	    if(verbosity >= 40) {
	      std::cout << "{ " << pos.X() << ", " << pos.Y();
	      std::cout << ", " << pos.Z() << " }, ";
	    }
	    fHistHitsSme[idet]->Fill(u_smear,v_smear);
	  }
	  fHistHitsTru[idet]->Fill(hvec.Phi(),hvecrad);
	  meas = new PHGenFit::PlanarMeasurement(pos, u, v, detphires, detlonres);
	  //foundorder.push_back(pos);
	} else {
	  std::cout << "   Type not implemented!" << std::endl;
	  return;
	}
	meas_out[distance] = meas;
	//meas->getMeasurement()->Print(); //DEBUG
      } /*ihit*/
    } /*ilayer*/
    if(verbosity>20)
      std::cout << "   Total Number Of Hits " << nhits << std::endl;
    //if(nhitsasoc>0)
    //  std::cout << fHitsNames[idet].c_str() <<"(" << nhitsasoc<< ")" << " ";
  } /*idet*/
  std::cout << std::endl;
return;
}
//========================================
SvtxTrack* PHG4EICTrackingBase::MakeSvtxTrack(const PHGenFit::Track* phgf_track,
					   const PHG4Particle* particle,
					   const unsigned int nmeas,
					   const TVector3& vtx) {
  TVector3 mct(particle->get_px(),particle->get_py(),particle->get_pz());

  //=== FIRST, EXTRAPOLATE PARAMETERS TO PRIMARY VERTEX
  //=== (along Z now, should change to DCA direction)
  unique_ptr<genfit::MeasuredStateOnPlane> gf_state(new genfit::MeasuredStateOnPlane());
  //if(!fVertexIn) {
  TVector3 dca = TVector3(0., 0., 1.);
  phgf_track->extrapolateToLine( *gf_state, vtx, dca);
  //there is a problem for the extraporation for
  //certain mid-rapidity tracks to the prim vertex 
  //for these cases a direction flip is being used
  //for the time being
  TVector3 trk = gf_state->getMom();

  double costheta = mct.Dot(trk)/mct.Mag()/trk.Mag();
  if( costheta<-0.1 ) {
    if(verbosity>200)
      std::cout << " ****FALSE Z!!! trying again " << std::endl;
    phgf_track->extrapolateToPlane( *gf_state, vtx, dca, 0);
  }
  trk = gf_state->getMom();
  costheta = mct.Dot(trk)/mct.Mag()/trk.Mag();
  if( costheta<-0.1 ) {
    //if( (gf_state->getMom()).Z()*particle_pz < 0 ) {
    if(verbosity>200)
      std::cout << " ****FALSE Z!!! flipping now " << std::endl;
    trk.SetX( -1*trk.X() );
    trk.SetY( -1*trk.Y() );
    trk.SetZ( -1*trk.Z() );
  }

  //=== SECOND, COPY TRACK PARAMETERS
  double chi2 = phgf_track->get_chi2();
  double ndf = phgf_track->get_ndf();
  TVector3 pos = gf_state->getPos();
  TMatrixDSym cov = gf_state->get6DCov();
  double dca2d = gf_state->getState()[3];
  double dca3d = sqrt(dca2d * dca2d + gf_state->getState()[4] * gf_state->getState()[4]);

  if(verbosity > 200) {
    std::cout << "  VERTEX " << vtx.X() << " " << vtx.Y() << " " << vtx.Z() << std::endl;
    std::cout << "  TRACKMOM " << trk.X() << " " << trk.Y() << " " << trk.Z() << std::endl;
    std::cout << "  TRACKPOS " << pos.X() << " " << pos.Y() << " " << pos.Z() << std::endl;
    std::cout << "  CHI2 " << chi2 << " ";
    std::cout << "  NDF " << ndf << std::endl;
    std::cout << "  DCA2 " << dca2d << " ";
    std::cout << "  DCA3 " << dca3d << std::endl;
    std::cout << "  PARTICLE " << particle->get_px() << " " << particle->get_py() << " " << particle->get_pz() << std::endl;
    std::cout << std::endl;
  }

  SvtxTrack_FastSim *out_track = new SvtxTrack_FastSim();
  out_track->set_truth_track_id( particle->get_track_id() );
  out_track->set_dca2d(dca2d);
  out_track->set_dca2d_error(gf_state->getCov()[3][3]);
  out_track->set_dca(dca3d);
  out_track->set_chisq(chi2);
  out_track->set_ndf(ndf);
  out_track->set_charge(phgf_track->get_charge());
  out_track->set_num_measurements(nmeas);
  out_track->set_px(trk.Px());
  out_track->set_py(trk.Py());
  out_track->set_pz(trk.Pz());
  out_track->set_x(pos.X());
  out_track->set_y(pos.Y());
  out_track->set_z(pos.Z());
  for (int i = 0; i != 6; ++i)
    for (int j = i; j != 6; ++j)
      out_track->set_error(i, j, cov[i][j]);
  return ((SvtxTrack*)out_track);
}
//========================================
PHGenFit::PlanarMeasurement* PHG4EICTrackingBase::MeasurementVertical(const PHG4Hit* g4hit,
								   const double phi_resolution,
								   const double r_resolution) {
  PHGenFit::PlanarMeasurement* meas = NULL;
  TVector3 pos(g4hit->get_avg_x(), g4hit->get_avg_y(), g4hit->get_avg_z());
  TVector3 v(pos.X(), pos.Y(), 0); v = 1 / v.Mag() * v; //========== r unit
  TVector3 u = v.Cross(TVector3(0, 0, 1)); u = 1 / u.Mag() * u; //== phi unit  
  if(true) {
    double u_smear = gRandom->Gaus(0, phi_resolution);
    double v_smear = gRandom->Gaus(0, r_resolution);
    pos.SetX(g4hit->get_avg_x() + u_smear * u.X() + v_smear * v.X());
    pos.SetY(g4hit->get_avg_y() + u_smear * u.Y() + v_smear * v.Y());
    if(verbosity >= 40) {
      std::cout << "{ " << pos.X() << ", " << pos.Y();
      std::cout << ", " << pos.Z() << " }, ";
    }
  }
  meas = new PHGenFit::PlanarMeasurement(pos, u, v, phi_resolution, r_resolution);
  return meas;
}
//========================================
PHGenFit::PlanarMeasurement* PHG4EICTrackingBase::MeasurementOnion(const PHG4Hit* g4hit,
								const double phi_resolution,
								const double z_resolution) {
  PHGenFit::PlanarMeasurement* meas = NULL;
  TVector3 pos(g4hit->get_avg_x(), g4hit->get_avg_y(), g4hit->get_avg_z());
  TVector3 v(0, 0, pos.Z()); v = 1 / v.Mag() * v; //=========================== z unit
  TVector3 u = v.Cross(TVector3(pos.X(), pos.Y(), 0)); u = 1 / u.Mag() * u; //= phi unit
  if(true) {
    double u_smear = gRandom->Gaus(0, phi_resolution);
    double v_smear = gRandom->Gaus(0, z_resolution);
    pos.SetX(g4hit->get_avg_x() + u_smear * u.X() + v_smear * v.X());
    pos.SetY(g4hit->get_avg_y() + u_smear * u.Y() + v_smear * v.Y());
    if(verbosity >= 40) {
      std::cout << "{ " << pos.X() << ", " << pos.Y();
      std::cout << ", " << pos.Z() << " }, ";
    }
  }
  meas = new PHGenFit::PlanarMeasurement(pos, u, v, phi_resolution, z_resolution);
  return meas;
}
