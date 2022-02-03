/// \file NovoScintHit.cc
/// \brief Implementation of the NovoScintHit class
//
//
#include "NovoScintHit.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"


G4ThreadLocal G4Allocator<NovoScintHit>* NovoScintHitAllocator=0;

NovoScintHit::NovoScintHit()
:	fEdep(0.),
	fTime(0),
	fPos(0.),
	fPhysVol(0),
	fIsPrimary(false)
{
}


NovoScintHit::NovoScintHit(G4VPhysicalVolume* pVol)
: fPhysVol(pVol)
{
}


NovoScintHit::~NovoScintHit()
{
}

const NovoScintHit& NovoScintHit::operator=(const NovoScintHit &right)
{
  fEdep = right.fEdep;
  fPos = right.fPos;
  fTime = right.fTime;
  fPhysVol = right.fPhysVol;
  fIsPrimary = right.fIsPrimary;
  return *this;
}

G4int NovoScintHit::operator==(const NovoScintHit&) const
{
  return false;
  //returns false because there currently isnt need to check for equality yet
}

void NovoScintHit::Draw()
{
}

void NovoScintHit::Print()
{
}
