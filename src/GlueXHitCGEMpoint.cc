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

#include "GlueXHitCGEMpoint.hh"

G4ThreadLocal G4Allocator<GlueXHitCGEMpoint>* GlueXHitCGEMpointAllocator = 0;

GlueXHitCGEMpoint::GlueXHitCGEMpoint(G4int layer)
 : G4VHit(),
   layer_(layer)
{}

GlueXHitCGEMpoint::GlueXHitCGEMpoint(const GlueXHitCGEMpoint &src)
{
   layer_ = src.layer_;
   E_GeV = src.E_GeV;
   dEdx_GeV_cm = src.dEdx_GeV_cm;
   dradius_cm = src.dradius_cm;
   primary_ = src.primary_;
   ptype_G3 = src.ptype_G3;
   px_GeV = src.px_GeV;
   py_GeV = src.py_GeV;
   pz_GeV = src.pz_GeV;
   x_cm = src.x_cm;
   y_cm = src.y_cm;
   z_cm = src.z_cm;
   t_ns = src.t_ns;
   track_ = src.track_;
   trackID_ = src.trackID_;
}

int GlueXHitCGEMpoint::operator==(const GlueXHitCGEMpoint &right) const
{
   if (E_GeV          != right.E_GeV       ||
       dEdx_GeV_cm    != right.dEdx_GeV_cm ||
       dradius_cm     != right.dradius_cm  ||
       primary_       != right.primary_    ||
       ptype_G3       != right.ptype_G3    ||
       px_GeV         != right.px_GeV      ||
       py_GeV         != right.py_GeV      ||
       pz_GeV         != right.pz_GeV      ||
       x_cm           != right.x_cm        ||
       y_cm           != right.y_cm        ||
       z_cm           != right.z_cm        ||
       t_ns           != right.t_ns        ||
       track_         != right.track_      ||
       trackID_       != right.trackID_    )
   {
      return 0;
   }
   return 1;
}

GlueXHitCGEMpoint &GlueXHitCGEMpoint::operator+=(const GlueXHitCGEMpoint &right)
{
   G4cerr << "Error in GlueXHitCGEMpoint::operator+= - "
          << "illegal attempt to merge two TruthPoint objects in the fdc!"
          << G4endl;
   return *this;
}

void GlueXHitCGEMpoint::Draw() const
{
   // not yet implemented
}

void GlueXHitCGEMpoint::Print() const
{
   G4cout << "GlueXHitCGEMpoint:" << G4endl
          << "   track = " << track_ << G4endl
          << "   trackID = " << trackID_ << G4endl
          << "   E = " << E_GeV << " GeV" << G4endl
          << "   dEdx = " << dEdx_GeV_cm << " GeV/cm" << G4endl
          << "   dradius = " << dradius_cm << " cm" << G4endl
          << "   primary = " << primary_ << G4endl
          << "   ptype = " << ptype_G3 << G4endl
          << "   px = " << px_GeV << " GeV/c" << G4endl
          << "   py = " << py_GeV << " GeV/c" << G4endl
          << "   pz = " << pz_GeV << " GeV/c" << G4endl
          << "   x = " << x_cm << " cm" << G4endl
          << "   y = " << y_cm << " cm" << G4endl
          << "   z = " << z_cm << " cm" << G4endl
          << "   t = " << t_ns << " ns" << G4endl
          << G4endl;
}

void printallhits(GlueXHitsMapCGEMpoint *hitsmap)
{
   std::map<int, GlueXHitCGEMpoint*> *map = hitsmap->GetMap();
   std::map<int, GlueXHitCGEMpoint*>::const_iterator iter;
   G4cout << "G4THitsMap " << hitsmap->GetName() << " with " << hitsmap->entries()
          << " entries:" << G4endl;
   for (iter = map->begin(); iter != map->end(); ++iter) {
     G4cout << "  key=" << iter->first << " ";
     iter->second->Print();
   }
}
