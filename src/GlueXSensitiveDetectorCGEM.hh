/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       GlueX and KLF Collaborations                    * 
*                                                                         *                                                                                                                               
* Author: The GlueX and KLF Collaborations                                *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/


#ifndef GlueXSensitiveDetectorCGEM_h
#define GlueXSensitiveDetectorCGEM_h 1

#include "G4VSensitiveDetector.hh"
#include "G4AutoLock.hh"

#include "GlueXHitCGEMhole.hh"
#include "GlueXHitCGEMpoint.hh"

class G4Step;
class G4HCofThisEvent;

class GlueXSensitiveDetectorCGEM : public G4VSensitiveDetector
{
 public:
   GlueXSensitiveDetectorCGEM(const G4String& name);
   GlueXSensitiveDetectorCGEM(const GlueXSensitiveDetectorCGEM &right);
   GlueXSensitiveDetectorCGEM &operator=(const GlueXSensitiveDetectorCGEM &right);
   virtual ~GlueXSensitiveDetectorCGEM();
  
   virtual void Initialize(G4HCofThisEvent* hitCollection);
   virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* unused);
   virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   int GetIdent(std::string div, const G4VTouchable *touch);

 private:
   GlueXHitsMapCGEMhole* fHoleHitsMap;
   GlueXHitsMapCGEMpoint* fPointsMap;

   std::map<G4LogicalVolume*, int> fVolumeTable;

   static int MAX_HITS;
   static double TWO_HIT_TIME_RESOL;
   static double THRESH_KEV;

   static int instanceCount;
   static G4Mutex fMutex;
};

#endif
