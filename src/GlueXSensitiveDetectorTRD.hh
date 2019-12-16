//
// GlueXSensitiveDetectorTRD - class header
//
// author: richard.t.jones at uconn.edu
// version: november 28, 2016
//          may 7, 2017 R. Dzhygadlo - added adjustments for
//                      GlueXHitTRDpoint GlueXHitTRDPmt    
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.

#ifndef GlueXSensitiveDetectorTRD_h
#define GlueXSensitiveDetectorTRD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitTRDpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorTRD : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorTRD(const G4String& name);
   GlueXSensitiveDetectorTRD(const GlueXSensitiveDetectorTRD &right);
   GlueXSensitiveDetectorTRD &operator=(const GlueXSensitiveDetectorTRD &right);
   virtual ~GlueXSensitiveDetectorTRD();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* ROhist);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

   static double GetDetectionEfficiency(double energy_GeV);
  
 private:
  std::vector<GlueXHitTRDpoint> fHitsPlane;
  
  std::map<G4LogicalVolume*, int> fVolumeTable;
  
  static int instanceCount;
  static G4Mutex fMutex;
};

#endif
