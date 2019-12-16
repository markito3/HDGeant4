//
// GlueXHitTRDpoint - class implementation
//

#include "GlueXHitTRDpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitTRDpoint>* GlueXHitTRDpointAllocator = 0;

GlueXHitTRDpoint::GlueXHitTRDpoint(const GlueXHitTRDpoint &src)
{
  E_GeV = src.E_GeV;
  t_ns = src.t_ns;
  x_cm = src.x_cm;
  y_cm = src.y_cm;
  z_cm = src.z_cm;
  px_GeV = src.px_GeV;
  py_GeV = src.py_GeV;
  pz_GeV = src.pz_GeV;
  pdg = src.pdg;
  plane = src.plane;
  track = src.track;
}

void GlueXHitTRDpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitTRDpoint::Print() const
{
   G4cout << "GlueXHitTRDpoint:" << G4endl
          << "   track = " << track << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
	  << "   t = " << t_ns << " ns" << G4endl
          << "   pdg = " << pdg << G4endl
	  << "   plane = " << plane << G4endl
          << "   px = " << px_GeV << " GeV/c" << G4endl
          << "   py = " << py_GeV << " GeV/c" << G4endl
          << "   pz = " << pz_GeV << " GeV/c" << G4endl
          << "   x = " << x_cm << " cm" << G4endl
          << "   y = " << y_cm << " cm" << G4endl
          << "   z = " << z_cm << " cm" << G4endl
          << G4endl;
}
