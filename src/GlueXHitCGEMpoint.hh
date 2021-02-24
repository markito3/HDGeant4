/**************************************************************************                                                                                                                           
* HallD software                                                          * 
* Copyright(C) 2020       GlueX and KLF Collaborations                    * 
*                                                                         *                                                                                                                               
* Author: The GlueX and KLF Collaborations                                *                                                                                                                                
* Contributors: Igal Jaegle                                               *                                                                                                                               
*                                                                         *                                                                                                                               
* This software is provided "as is" without any warranty.                 *
**************************************************************************/
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitCGEMpoint_h
#define GlueXHitCGEMpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitCGEMpoint : public G4VHit
{
 public:
   GlueXHitCGEMpoint() {}
   GlueXHitCGEMpoint(int layer);
   GlueXHitCGEMpoint(const GlueXHitCGEMpoint &src);
   int operator==(const GlueXHitCGEMpoint &right) const;
   GlueXHitCGEMpoint &operator+=(const GlueXHitCGEMpoint &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   G4int layer_;        // CGEM layer, count from upstream starting at 1

   G4double E_GeV;       // energy of primary track at point
   G4double dEdx_GeV_cm; // dE/dx (GeV/cm) of this track inside straw
   G4double dradius_cm;  // track distance (cm) of closest approach to wire
   G4bool primary_;      // true if track belongs to from a primary particle
   G4int ptype_G3;       // G3 type of particle making this track
   G4double px_GeV;      // momentum (GeV/c) of track at point, x component
   G4double py_GeV;      // momentum (GeV/c) of track at point, y component
   G4double pz_GeV;      // momentum (GeV/c) of track at point, z component
   G4double x_cm;        // global x coordinate of track at point (cm)
   G4double y_cm;        // global y coordinate of track at point (cm)
   G4double z_cm;        // global z coordinate of track at point (cm)
   G4double t_ns;        // time of track crossing at point (ns)
   G4int track_;         // Geant4 track ID of particle making this track
   G4int trackID_;       // GlueX-assigned track ID of particle making this track

   G4int GetKey() const { return (track_ << 20) + int(t_ns * 100); }
};

typedef G4THitsMap<GlueXHitCGEMpoint> GlueXHitsMapCGEMpoint;

extern G4ThreadLocal G4Allocator<GlueXHitCGEMpoint>* GlueXHitCGEMpointAllocator;

inline void* GlueXHitCGEMpoint::operator new(size_t)
{
   if (!GlueXHitCGEMpointAllocator)
      GlueXHitCGEMpointAllocator = new G4Allocator<GlueXHitCGEMpoint>;
   return (void *) GlueXHitCGEMpointAllocator->MallocSingle();
}

inline void GlueXHitCGEMpoint::operator delete(void *aHit)
{
   GlueXHitCGEMpointAllocator->FreeSingle((GlueXHitCGEMpoint*) aHit);
}

#endif
