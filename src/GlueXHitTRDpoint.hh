//
// GlueXHitTRDpoint - class implementation
//

#ifndef GlueXHitTRDpoint_h
#define GlueXHitTRDpoint_h 1

#include "G4VHit.hh"
#include "G4THitsMap.hh"
#include "G4Allocator.hh"

class GlueXHitTRDpoint : public G4VHit
{
public:
  GlueXHitTRDpoint() {}
  GlueXHitTRDpoint(const GlueXHitTRDpoint &src);

  void *operator new(size_t);
  void operator delete(void *aHit);

  void Draw() const;
  void Print() const;

  G4double E_GeV;      // track particle total energy (GeV)
  G4double t_ns;       // pulse leading-edge time (ns)
  G4double x_cm;       // x coordinate where ch. track hits the radiator
  G4double y_cm;        
  G4double z_cm;        
  G4double px_GeV;     // px component of the track momentum
  G4double py_GeV;     
  G4double pz_GeV;
  G4int pdg;           // PDG of the particle
  G4int plane;         // index of the TRD plane
  G4int track;         // index of the MC track 
};

typedef G4THitsMap<GlueXHitTRDpoint> GlueXHitsMapTRDpoint;

extern G4ThreadLocal G4Allocator<GlueXHitTRDpoint>* GlueXHitTRDpointAllocator;

inline void* GlueXHitTRDpoint::operator new(size_t)
{
  if (!GlueXHitTRDpointAllocator)
    GlueXHitTRDpointAllocator = new G4Allocator<GlueXHitTRDpoint>;
  return (void *) GlueXHitTRDpointAllocator->MallocSingle();
}

inline void GlueXHitTRDpoint::operator delete(void *aHit)
{
  GlueXHitTRDpointAllocator->FreeSingle((GlueXHitTRDpoint*) aHit);
}

#endif
