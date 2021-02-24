//
// GlueXHitCGEMlayer@ - class header
//
// author: richard.t.jones at uconn.edu
// version: october 21, 2016
//
// In the context of the Geant4 event-level multithreading model,
// this class is "thread-local", ie. has thread-local state. Its
// allocator is designed to run within a worker thread context.
// This class is final, do NOT try to derive another class from it.

#ifndef GlueXHitCGEMlayer_h
#define GlueXHitCGEMlayer_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitCGEMlayer : public G4VHit
{
 public:
   GlueXHitCGEMlayer() {}
   GlueXHitCGEMlayer(G4int layer);
   GlueXHitCGEMlayer(const GlueXHitCGEMlayer &src);
   int operator==(const GlueXHitCGEMlayer &right) const;
   GlueXHitCGEMlayer &operator+=(const GlueXHitCGEMlayer &right);

   void *operator new(size_t);
   void operator delete(void *aHit);

   void Draw() const;
   void Print() const;

   // no reason to hide hit data

   G4int layer_;          // counter number, from 1 at/past phi=0

   struct hitinfo_t {
      G4double dE_MeV;     // energy deposition (MeV)
      G4double t_ns;       // pulse leading-edge time (ns)
      G4double itrack_;    // track index of first particle making this hit
      G4double ptype_G3;   // G3 type of first particle making this hit
      G4double t0_ns;      // time of passage of the track making this hit
      G4double x_cm;       // z coordinate of the hit in global refsys
      G4double y_cm;       // z coordinate of the hit in global refsys
      G4double z_cm;       // z coordinate of the hit in global refsys
   };
   std::vector<hitinfo_t> hits;

   G4int GetKey() const { return GetKey(layer_); }
   static G4int GetKey(G4int layer) {
      return layer;
   }
};

typedef G4THitsMap<GlueXHitCGEMlayer> GlueXHitsMapCGEMlayer;

extern G4ThreadLocal G4Allocator<GlueXHitCGEMlayer>* GlueXHitCGEMlayerAllocator;

inline void* GlueXHitCGEMlayer::operator new(size_t)
{
   if (!GlueXHitCGEMlayerAllocator)
      GlueXHitCGEMlayerAllocator = new G4Allocator<GlueXHitCGEMlayer>;
   return (void *) GlueXHitCGEMlayerAllocator->MallocSingle();
}

inline void GlueXHitCGEMlayer::operator delete(void *aHit)
{
   GlueXHitCGEMlayerAllocator->FreeSingle((GlueXHitCGEMlayer*) aHit);
}

#endif
