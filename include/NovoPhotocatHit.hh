/// \file NovoPhotocatHit.hh
/// \brief Definition of the NovoPhotoCatHit class
//
//
#ifndef NovoPhotocatHit_h
#define NovoPhotocatHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"

#include "tls.hh"
#include <vector>

class G4VTouchable;

class NovoPhotocatHit : public G4VHit
{
	public:
		NovoPhotocatHit();
		virtual ~NovoPhotocatHit();
		NovoPhotocatHit(const NovoPhotocatHit &right);

		const NovoPhotocatHit& operator=(const NovoPhotocatHit &right);
		G4int operator==(const NovoPhotocatHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void *aHit);
 
		virtual void Draw();
		virtual void Print();

		inline void SetDrawit(G4bool b){fDrawit=b;}
		inline G4bool GetDrawit(){return fDrawit;}


		inline void IncPhotonCount(){fPhotons++;}
		inline G4int GetPhotonCount(){return fPhotons;}

		inline void SetPhotocatNumber(G4int n) { fPhotocatNumber = n; }
		inline G4int GetPhotocatNumber() { return fPhotocatNumber; }

		inline void SetPhotocatPhysVol(G4VPhysicalVolume* physVol){this->fPhysVol=physVol;}
		inline G4VPhysicalVolume* GetPhotocatPhysVol(){return fPhysVol;}

		inline void SetPhotocatPos(G4double x,G4double y,G4double z)
		{
			//~ G4cout << "SetPhotocatPos: " << x << " " << y << " " << z << G4endl; 
			fPos=G4ThreeVector(x,y,z);
		}
 
		inline G4ThreeVector GetPhotocatPos(){return fPos;}
		
		inline void SetTime(G4double t){fTimestamps.push_back(t);}
		inline G4DataVector GetTimestamps(){return fTimestamps;}
		
		inline void IncEdep(G4double de){fEdep += de;}
		inline G4double GetEdep(){return fEdep;}
		
		inline void SetPhotonE(G4double e){fPhotonE.push_back(e);}
		inline G4DataVector GetPhotonE(){return fPhotonE;}

	private:
		G4int fPhotocatNumber;
		G4int fPhotons;
		G4ThreeVector fPos;
		G4VPhysicalVolume* fPhysVol;
		G4bool fDrawit;
		G4DataVector fTimestamps;
		G4double fEdep;
		G4DataVector fPhotonE;
};

typedef G4THitsCollection<NovoPhotocatHit> NovoPhotocatHitsCollection;

extern G4ThreadLocal G4Allocator<NovoPhotocatHit>* NovoPhotocatHitAllocator;

inline void* NovoPhotocatHit::operator new(size_t)
{
	if(!NovoPhotocatHitAllocator)
		NovoPhotocatHitAllocator = new G4Allocator<NovoPhotocatHit>;
	return (void *) NovoPhotocatHitAllocator->MallocSingle();
}

inline void NovoPhotocatHit::operator delete(void *aHit)
{
	NovoPhotocatHitAllocator->FreeSingle((NovoPhotocatHit*) aHit);
}

#endif
