//
// GlueXSensitiveDetectorCGEM - class implementation
//
// author: richard.t.jones at uconn.edu
// version: october 11, 2016

#include "GlueXSensitiveDetectorCGEM.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"

#include <CLHEP/Random/RandPoisson.h>
#include <Randomize.hh>

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include <JANA/JApplication.h>

#include <stdlib.h>
#include <math.h>

G4Mutex GlueXSensitiveDetectorCGEM::fMutex = G4MUTEX_INITIALIZER;
int GlueXSensitiveDetectorCGEM::instanceCount = 0;

// Cutoff on the total number of allowed hits
int GlueXSensitiveDetectorCGEM::MAX_HITS = 1000;
// Minimum hit time difference for two hits on the same cell
double GlueXSensitiveDetectorCGEM::TWO_HIT_TIME_RESOL = 10.0*ns;
// Minimum energy deposition for a hit
double GlueXSensitiveDetectorCGEM::THRESH_MEV = 0.0;

GlueXSensitiveDetectorCGEM::GlueXSensitiveDetectorCGEM(const G4String& name)
 : G4VSensitiveDetector(name), 
   fLayersMap(0), fPointsMap(0)
{
   collectionName.insert("CGEMLayersCollection");
   collectionName.insert("CGEMPointsCollection");

   // The rest of this only needs to happen once, the first time an object
   // of this type is instantiated for this configuration of geometry and
   // fields. If the geometry or fields change in such a way as to modify
   // the drift-time properties of hits in the CGEM, you must delete all old
   // objects of this class and create new ones.

   G4AutoLock barrier(&fMutex);
   if (instanceCount++ == 0) {
   }
}

GlueXSensitiveDetectorCGEM::GlueXSensitiveDetectorCGEM(
                     const GlueXSensitiveDetectorCGEM &src)
  : G4VSensitiveDetector(src), 
    fLayersMap(src.fLayersMap), fPointsMap(src.fPointsMap)
{
   G4AutoLock barrier(&fMutex);
   ++instanceCount;
}

GlueXSensitiveDetectorCGEM &GlueXSensitiveDetectorCGEM::operator=(const
                                         GlueXSensitiveDetectorCGEM &src)
{
   G4AutoLock barrier(&fMutex);
   *(G4VSensitiveDetector*)this = src;
   fLayersMap = src.fLayersMap;
   fPointsMap = src.fPointsMap;
   return *this;
}

GlueXSensitiveDetectorCGEM::~GlueXSensitiveDetectorCGEM() 
{
   G4AutoLock barrier(&fMutex);
   --instanceCount;
}

void GlueXSensitiveDetectorCGEM::Initialize(G4HCofThisEvent* hce)
{
   fLayersMap = new
               GlueXHitsMapCGEMlayer(SensitiveDetectorName, collectionName[0]);
   fPointsMap = new
               GlueXHitsMapCGEMpoint(SensitiveDetectorName, collectionName[0]);
   G4SDManager *sdm = G4SDManager::GetSDMpointer();
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fLayersMap);
   hce->AddHitsCollection(sdm->GetCollectionID(collectionName[0]), fPointsMap);
}

G4bool GlueXSensitiveDetectorCGEM::ProcessHits(G4Step* step, 
                                              G4TouchableHistory* ROhist)
{

   double dEsum = step->GetTotalEnergyDeposit();
   if (dEsum == 0)
      return false;

   const G4ThreeVector &pin = step->GetPreStepPoint()->GetMomentum();
   const G4ThreeVector &xin = step->GetPreStepPoint()->GetPosition();
   const G4ThreeVector &xout = step->GetPostStepPoint()->GetPosition();
   double Ein = step->GetPreStepPoint()->GetTotalEnergy();
   double tin = step->GetPreStepPoint()->GetGlobalTime();
   double tout = step->GetPostStepPoint()->GetGlobalTime();
   G4ThreeVector x = (xin + xout) / 2;
   G4ThreeVector dx = xout - xin;
   double t = (tin + tout) / 2;
   double dr = dx.mag();
   double dEdx = (dr > 1e-3*cm)? dEsum/dr : 0;

   const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
   const G4AffineTransform &local_from_global = touch->GetHistory()
                                                     ->GetTopTransform();
   G4ThreeVector xlocal = local_from_global.TransformPoint(x);
   G4ThreeVector xinlocal = local_from_global.TransformPoint(xin);
   G4ThreeVector xoutlocal = local_from_global.TransformPoint(xout);

   // For particles that range out inside the active volume, the
   // "out" time may sometimes be set to something enormously high.
   // This screws up the hit. Check for this case here by looking
   // at tout and making sure it is less than 1 second. If it's
   // not, then just use tin for "t".

   if (tout > 1.0*s)
      t = tin;
   
   int layer = GetIdent("layer", touch);
   if (layer == 0) {
     printf("hitCGEM: CGEM package number evaluates to zero! "
	    "THIS SHOULD NEVER HAPPEN! drop this particle.\n");
     return false;
   }
   
   // Normally numeric identifiers start at 1, eg. layer, package, module
   // but if it is an index counting from zero, add the "No" suffix.
   int layerNo = layer - 1;
   
   G4Track *track = step->GetTrack();
   int trackID = track->GetTrackID();
   int pdgtype = track->GetDynamicParticle()->GetPDGcode();
   int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
   int dradius = 0;
   // Post the hit to the points list in the
   // order of appearance in the event simulation.
 
   GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*) track->GetUserInformation();
   int itrack = trackinfo->GetGlueXTrackID();
   if (trackinfo->GetGlueXHistory() == 0 && itrack > 0) {
     G4int key = fPointsMap->entries();
     GlueXHitCGEMpoint* lastPoint = (*fPointsMap)[key - 1];
     // Limit CGEM truthPoints to one per layer
     //if (lastPoint == 0 || lastPoint->track_ != trackID ||
     //lastPoint->layer_ != layerNo)
     //{
     GlueXHitCGEMpoint newPoint(layer);
     newPoint.primary_ = (track->GetParentID() == 0);
     newPoint.track_ = trackID;
     newPoint.x_cm = xout[0]/cm;
     newPoint.y_cm = xout[1]/cm;
     newPoint.z_cm = xout[2]/cm;
     newPoint.t_ns = tout/ns;
     newPoint.px_GeV = pin[0]/GeV;
     newPoint.py_GeV = pin[1]/GeV;
     newPoint.pz_GeV = pin[2]/GeV;
     newPoint.E_GeV = Ein/GeV;
     dradius = sqrt(pow(xout[0], 2) + pow(xout[1], 2));
     newPoint.dradius_cm = dradius/cm;
     newPoint.dEdx_GeV_cm = dEdx/(GeV/cm);
     newPoint.ptype_G3 = g3type;
     newPoint.trackID_ = itrack;
     fPointsMap->add(key, newPoint);
     trackinfo->SetGlueXHistory(2);
     //  }
   }
   
   // Post the hit to the hits map, ordered by layer index
   if (dEsum > 0) {
     int key = GlueXHitCGEMlayer::GetKey(layer);
     GlueXHitCGEMlayer *gap = (*fLayersMap)[key];
     if (gap == 0) {
       GlueXHitCGEMlayer newgap(layer);
       fLayersMap->add(key, newgap);
       gap = (*fLayersMap)[key];
     }

      // Add the hit to the hits vector, maintaining strict time ordering

      int merge_hit = 0;
      std::vector<GlueXHitCGEMlayer::hitinfo_t>::iterator hiter;
      
      for (hiter = gap->hits.begin(); hiter != gap->hits.end(); ++hiter) {
	if (fabs(hiter->t_ns*ns - t) < TWO_HIT_TIME_RESOL) {
	    merge_hit = 1;
            break;
	} else if (hiter->t_ns*ns > t) {
	  break;
	}
      }
      
      if (merge_hit) {
	// Use the time from the earlier hit but add the charge
	hiter->dE_MeV += dEsum/MeV;
	if (hiter->t_ns*ns > t) {
	  //hiter->layer = layer;
	  hiter->t_ns = t/ns;
	  //hiter->itrack_ = itrack;
	  //hiter->ptype_G3 = g3type;
	  hiter->t0_ns = tin/ns;
	  hiter->x_cm = x[0]/cm;
	  hiter->y_cm = x[1]/cm;
	  hiter->z_cm = x[2]/cm;
	}
      } else if ((int)gap->hits.size() < MAX_HITS)	{
         // create new hit 
         hiter = gap->hits.insert(hiter, GlueXHitCGEMlayer::hitinfo_t());
         hiter->dE_MeV = dEsum/MeV;
	 hiter->t_ns = t/ns;
         //hiter->itrack_ = itrack;
         //hiter->ptype_G3 = g3type;
         hiter->t0_ns = tin/ns;
         hiter->x_cm = x[0]/cm;
	 hiter->y_cm = x[1]/cm;
	 hiter->z_cm = x[2]/cm;
      } else {
	G4cerr << "GlueXSensitiveDetectorSTC::ProcessHits error: "
	       << "max hit count " << MAX_HITS << " exceeded, truncating!"
	       << G4endl;
      }
   }
   return true;
}

void GlueXSensitiveDetectorCGEM::EndOfEvent(G4HCofThisEvent*)
{
   std::map<int,GlueXHitCGEMlayer*> *layers = fLayersMap->GetMap();
   std::map<int,GlueXHitCGEMpoint*> *points = fPointsMap->GetMap();
   if (layers->size() == 0 && points->size() == 0)
      return;
   std::map<int,GlueXHitCGEMpoint*>::iterator piter;
   std::map<int,GlueXHitCGEMlayer*>::iterator siter;

   if (verboseLevel > 1) { 
      G4cout << G4endl
             << "--------> Hits Collection: in this event there are "
             << points->size() << " truth points in the CGEM: "
             << G4endl;
      for (piter = points->begin(); piter != points->end(); ++piter)
         piter->second->Print();

      G4cout << G4endl
             << "--------> Layer Hits Collection: in this event there are "
             << layers->size() << " truth points in the CGEM: "
             << G4endl;
      for (siter = layers->begin(); siter != layers->end(); ++siter)
         siter->second->Print();
   }

   // pack hits into ouptut hddm record
 
   G4EventManager* mgr = G4EventManager::GetEventManager();
   G4VUserEventInformation* info = mgr->GetUserInformation();
   hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
   if (record == 0) {
      G4cerr << "GlueXSensitiveDetectorCGEM::EndOfEvent error - "
             << "hits seen but no output hddm record to save them into, "
             << "cannot continue!" << G4endl;
      exit(1);
   }

   if (record->getPhysicsEvents().size() == 0) 
      record->addPhysicsEvents();
   if (record->getHitViews().size() == 0) 
      record->getPhysicsEvent().addHitViews();
   hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
   if (hitview.getCylindricalGEMs().size() == 0)
      hitview.addCylindricalGEMs();
   hddm_s::CylindricalGEM &cGEM = hitview.getCylindricalGEM();

   // Collect and output the stcTruthHits
   for (siter = layers->begin(); siter != layers->end(); ++siter) {
   
     std::vector<GlueXHitCGEMlayer::hitinfo_t> &hits = siter->second->hits;
     // apply a pulse height threshold cut
     for (unsigned int ih=0; ih < hits.size(); ++ih) {
       if (hits[ih].dE_MeV <= THRESH_MEV) {
	 hits.erase(hits.begin() + ih);
	 --ih;
       }
     }
     if (hits.size() > 0) {
       hddm_s::CgemLayerList layer = cGEM.addCgemLayers(1);
       layer(0).setLayer(siter->second->layer_);
       for (int ih=0; ih < (int)hits.size(); ++ih) {
	 hddm_s::CgemTruthHitList thit = layer(0).addCgemTruthHits(1);
	 thit(0).setDE(hits[ih].dE_MeV*MeV/GeV);
	 thit(0).setT(hits[ih].t_ns);
	 thit(0).setX(hits[ih].x_cm);
	 thit(0).setY(hits[ih].y_cm);
	 thit(0).setZ(hits[ih].z_cm);
       }
     }
   }
   
   // Collect and output the CGemTruthPoints

   int last_layer = -1;
   hddm_s::CgemTruthPoint *last_point = 0;
   for (piter = points->begin(); piter != points->end(); ++piter) {
     hddm_s::CgemTruthPointList point = cGEM.addCgemTruthPoints(1);
     //cout << " in cgem " << piter->second->z_cm << endl;
     point(0).setE(piter->second->E_GeV);
     point(0).setDEdx(piter->second->dEdx_GeV_cm);
     point(0).setDradius(piter->second->dradius_cm);
     point(0).setPrimary(piter->second->primary_);
     point(0).setPtype(piter->second->ptype_G3);
     point(0).setPx(piter->second->px_GeV);
     point(0).setPy(piter->second->py_GeV);
     point(0).setPz(piter->second->pz_GeV);
     point(0).setT(piter->second->t_ns);
     point(0).setX(piter->second->x_cm);
     point(0).setY(piter->second->y_cm);
     point(0).setZ(piter->second->z_cm);
     point(0).setTrack(piter->second->track_);
     hddm_s::TrackIDList tid = point(0).addTrackIDs();
     tid(0).setItrack(piter->second->trackID_);
     last_layer = piter->second->layer_;
     last_point = &point(0);
   }

}

int GlueXSensitiveDetectorCGEM::GetIdent(std::string div, 
                                        const G4VTouchable *touch)
{

   const HddsG4Builder* bldr = GlueXDetectorConstruction::GetBuilder();
   std::map<std::string, std::vector<int> >::const_iterator iter;
   std::map<std::string, std::vector<int> > *identifiers;
   int max_depth = touch->GetHistoryDepth();
   for (int depth = 0; depth < max_depth; ++depth) {
      G4VPhysicalVolume *pvol = touch->GetVolume(depth);
      G4LogicalVolume *lvol = pvol->GetLogicalVolume();
      int volId = fVolumeTable[lvol];
      if (volId == 0) {
         volId = bldr->getVolumeId(lvol);
         fVolumeTable[lvol] = volId;
      }
      identifiers = &Refsys::fIdentifierTable[volId];
      if ((iter = identifiers->find(div)) != identifiers->end()) {
         int copyNum = touch->GetCopyNumber(depth);
         copyNum += (dynamic_cast<G4PVPlacement*>(pvol))? -1 : 0;
         return iter->second[copyNum];
      }
   }

   return -1;
}
