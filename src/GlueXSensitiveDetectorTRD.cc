//
// GlueXSensitiveDetectorTRD - class implementation
//

#include "GlueXSensitiveDetectorTRD.hh"
#include "GlueXDetectorConstruction.hh"
#include "GlueXPrimaryGeneratorAction.hh"
#include "GlueXUserEventInformation.hh"
#include "GlueXUserTrackInformation.hh"
#include "GlueXUserOptions.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4TransportationManager.hh"
#include "G4ParallelWorldProcess.hh"

#include <JANA/JApplication.h>

int GlueXSensitiveDetectorTRD::instanceCount = 0;
G4Mutex GlueXSensitiveDetectorTRD::fMutex = G4MUTEX_INITIALIZER;

GlueXSensitiveDetectorTRD::GlueXSensitiveDetectorTRD(const G4String& name)
  : G4VSensitiveDetector(name)
{
  // The rest of this only needs to happen once, the first time an object
  // of this type is instantiated for this configuration of geometry and
  // fields. If the geometry or fields change in such a way as to modify
  // the drift-time properties of hits in the TRD, you must delete all old
  // objects of this class and create new ones.

  G4AutoLock barrier(&fMutex);
  if (instanceCount++ == 0) {
    extern int run_number;
    extern jana::JApplication *japp;
    if (japp == 0) {
      G4cerr << "Error in GlueXSensitiveDetector constructor - "
         << "jana global DApplication object not set, "
         << "cannot continue." << G4endl;
      exit(-1);
    }
    jana::JCalibration *jcalib = japp->GetJCalibration(run_number);
    if (japp == 0) {   // dummy
      jcalib = 0;
      G4cout << "TRD: ALL parameters loaded from ccdb" << G4endl;
    }
  }

  //GlueXUserOptions *user_opts = GlueXUserOptions::GetInstance();
}

GlueXSensitiveDetectorTRD::GlueXSensitiveDetectorTRD(const GlueXSensitiveDetectorTRD &src)
      : G4VSensitiveDetector(src)
{
  G4AutoLock barrier(&fMutex);
  ++instanceCount;
}

GlueXSensitiveDetectorTRD &GlueXSensitiveDetectorTRD::operator=(const
                                  GlueXSensitiveDetectorTRD &src)
{
  G4AutoLock barrier(&fMutex);
  *(G4VSensitiveDetector*)this = src;
  return *this;
}

GlueXSensitiveDetectorTRD::~GlueXSensitiveDetectorTRD() 
{
  G4AutoLock barrier(&fMutex);
  --instanceCount;
}

void GlueXSensitiveDetectorTRD::Initialize(G4HCofThisEvent* hce)
{
}

G4bool GlueXSensitiveDetectorTRD::ProcessHits(G4Step* step, 
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

  const G4VTouchable* touch = step->GetPreStepPoint()->GetTouchable();
  const G4AffineTransform &local_from_global = touch->GetHistory()->GetTopTransform();
  G4ThreeVector xlocal = local_from_global.TransformPoint(x);
  
  // For particles that range out inside the active volume, the
  // "out" time may sometimes be set to something enormously high.
  // This screws up the hit. Check for this case here by looking
  // at tout and making sure it is less than 1 second. If it's
  // not, then just use tin for "t".

  if (tout > 1.0*s) t = tin;

  // Post the hit to the points list in the
  // order of appearance in the event simulation.

  G4Track *track = step->GetTrack();
  GlueXUserTrackInformation *trackinfo = (GlueXUserTrackInformation*)
    track->GetUserInformation();
  int itrack = trackinfo->GetGlueXTrackID();
  G4String volname = touch->GetVolume()->GetName();  
 
  // sensitive gas volume: TRDG
  if (itrack > 0) { //trackinfo->GetGlueXHistory() == 0 && xin.dot(pin) > 0) {
      G4cout << "Found TRD sensitive gas with charged particle" << G4endl;
      int pdgtype = track->GetDynamicParticle()->GetPDGcode();
      int g3type = GlueXPrimaryGeneratorAction::ConvertPdgToGeant3(pdgtype);
      
      int plane = GetIdent("plane", touch);

      GlueXHitTRDpoint trdhit;
      trdhit.E_GeV = Ein/GeV;
      trdhit.t_ns = t/ns;
      trdhit.x_cm = x[0]/cm;
      trdhit.y_cm = x[1]/cm;
      trdhit.z_cm = x[2]/cm;
      trdhit.px_GeV = pin[0]/GeV;
      trdhit.py_GeV = pin[1]/GeV;
      trdhit.pz_GeV = pin[2]/GeV;
      trdhit.pdg = g3type;
      trdhit.plane = plane; // from HDDS geometry
      trdhit.track = itrack; // track id of the charged particle
      fHitsPlane.push_back(trdhit);
      return true;
  }

  return true;
}

void GlueXSensitiveDetectorTRD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) { 
    G4cout << G4endl
       << "--------> Hits Collection: in this event there are "
       << fHitsPlane.size() << " bar hits:"
       << G4endl;
    for(unsigned int h=0; h<fHitsPlane.size(); h++)
      fHitsPlane[h].Print();
  }

  // pack hits into ouptut hddm record
 
  G4EventManager* mgr = G4EventManager::GetEventManager();
  G4VUserEventInformation* info = mgr->GetUserInformation();
  hddm_s::HDDM *record = ((GlueXUserEventInformation*)info)->getOutputRecord();
  if (record == 0) {
    G4cerr << "GlueXSensitiveDetectorTRD::EndOfEvent error - "
           << "hits seen but no output hddm record to save them into, "
       << "cannot continue!" << G4endl;
    exit(1);
  }

  if (record->getPhysicsEvents().size() == 0) record->addPhysicsEvents();
  if (record->getHitViews().size() == 0) record->getPhysicsEvent().addHitViews();
  hddm_s::HitView &hitview = record->getPhysicsEvent().getHitView();
  if (hitview.getTRDs().size() == 0) hitview.addTRDs();
  hddm_s::TRD &trd = hitview.getTRD();

  // Collect and output the PlaneHits
  for(unsigned int h=0; h<fHitsPlane.size(); h++){
    hddm_s::TrdTruthPlaneHitList bhit = trd.addTrdTruthPlaneHits(1);
    bhit(0).setE(fHitsPlane[h].E_GeV);
    bhit(0).setT(fHitsPlane[h].t_ns);
    bhit(0).setX(fHitsPlane[h].x_cm);
    bhit(0).setY(fHitsPlane[h].y_cm);
    bhit(0).setZ(fHitsPlane[h].z_cm);
    bhit(0).setPx(fHitsPlane[h].px_GeV);
    bhit(0).setPy(fHitsPlane[h].py_GeV);
    bhit(0).setPz(fHitsPlane[h].pz_GeV);
    bhit(0).setPdg(fHitsPlane[h].pdg);
    bhit(0).setPlane(fHitsPlane[h].plane);
    bhit(0).setTrack(fHitsPlane[h].track);
  }

  fHitsPlane.clear();
}

int GlueXSensitiveDetectorTRD::GetIdent(std::string div, 
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
