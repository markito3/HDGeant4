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

#ifndef GlueXHitCGEMhole_h
#define GlueXHitCGEMhole_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitCGEMhole : public G4VHit
{
 public:
   GlueXHitCGEMhole() {}
   GlueXHitCGEMhole(G4int layer, G4int hole);
   GlueXHitCGEMhole(const GlueXHitCGEMhole &src);
   int operator==(const GlueXHitCGEMhole &right) const;
   GlueXHitCGEMhole &operator+=(const GlueXHitCGEMhole &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int layer_;           // layer number, from 1 upstream to down
   G4int hole_;            // hole number, from 1 low to high u

   struct hitinfo_t {
      G4double dE_keV;     // energy loss (keV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double dx_cm;      // track length in active region
      G4int itrack_;       // number of track creating the hit
      G4double z0_cm;      // track global coordinate (cm)
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(layer_, hole_); }
   static G4int GetKey(G4int layer, G4int hole) {
      return (layer << 10) + hole + 1;
   }
};

typedef G4THitsMap<GlueXHitCGEMhole> GlueXHitsMapCGEMhole;

extern G4ThreadLocal G4Allocator<GlueXHitCGEMhole>* GlueXHitCGEMholeAllocator;

inline void* GlueXHitCGEMhole::operator new(size_t)
{
   if (!GlueXHitCGEMholeAllocator)
      GlueXHitCGEMholeAllocator = new G4Allocator<GlueXHitCGEMhole>;
   return (void *) GlueXHitCGEMholeAllocator->MallocSingle();
}

inline void GlueXHitCGEMhole::operator delete(void *aHit)
{
   GlueXHitCGEMholeAllocator->FreeSingle((GlueXHitCGEMhole*) aHit);
}

#endif
