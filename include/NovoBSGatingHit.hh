/// \file NovoBSGatingHit.hh
/// \brief Definition of the NovoBSGatingHit class
/// In progress....


#ifndef NovoBSGatingHit_h
#define NovoBSGatingHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"


#include "tls.hh"

class NovoBSGatingHit : public G4VHit
{
	public:
		NovoBSGatingHit(); //constructor
		NovoBSGatingHit(G4VPhysicalVolume* pVol); //constructor ??
		virtual ~NovoBSGatingHit(); //Deconstructor
		NovoBSGatingHit(const NovoBSGatingHit &right); //??
		const NovoBSGatingHit& operator=(const NovoBSGatingHit &right);
		G4int operator==(const NovoBSGatingHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void *aHit);

		virtual void Draw();
		virtual void Print();

		inline void SetEdep(G4double de) { fEdep = de; }
		inline void AddEdep(G4double de) { fEdep += de; }
		inline G4double GetEdep() const { return fEdep; }

		inline void SetPos(G4ThreeVector xyz) { fPos = xyz; }
		inline G4ThreeVector GetPos() const { return fPos; }

		inline const G4VPhysicalVolume * GetPhysV() { return fPhysVol; }

		inline G4bool 	GetIsPrimary () const {return fIsPrimary;}
		inline void 	SetIsPrimary (G4bool isPrimary) {fIsPrimary=isPrimary;}

	private:
		G4double fEdep;
		G4ThreeVector fPos;
		const G4VPhysicalVolume* fPhysVol;
		G4bool 	fIsPrimary;
};

typedef G4THitsCollection<NovoBSGatingHit> NovoBSGatingHitsCollection;

extern G4ThreadLocal G4Allocator<NovoBSGatingHit>* NovoBSGatingHitAllocator; //??

inline void* NovoBSGatingHit::operator new(size_t)
{
  if(!NovoBSGatingHitAllocator)
      NovoBSGatingHitAllocator = new G4Allocator<NovoBSGatingHit>;
  return (void *) NovoBSGatingHitAllocator->MallocSingle();
}

inline void NovoBSGatingHit::operator delete(void *aHit)
{
  NovoBSGatingHitAllocator->FreeSingle((NovoBSGatingHit*) aHit);
}


#endif
