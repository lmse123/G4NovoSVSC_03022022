/// \file NovoOpticalPhysics.hh
/// \brief Definition of the NovoOpticalPhysics class
//


#ifndef NovoOpticalPhysics_h
#define NovoOpticalPhysics_h 1


#include "globals.hh"

#include "G4OpWLS.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"

#include "G4OpMieHG.hh"
#include "G4OpRayleigh.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4VPhysicsConstructor.hh"


class NovoOpticalPhysics : public G4VPhysicsConstructor
{
	public:

		NovoOpticalPhysics(const G4String& name ="Optical");
		virtual ~NovoOpticalPhysics();

		virtual void ConstructParticle();
		virtual void ConstructProcess();

		G4OpWLS* GetWLSProcess() {return fWLSProcess;}
		G4Cerenkov* GetCerenkovProcess() {return fCerenkovProcess;}
		G4Scintillation* GetScintillationProcess() {return fScintProcess;}
		G4OpAbsorption* GetAbsorptionProcess() {return fAbsorptionProcess;}
		G4OpRayleigh* GetRayleighScatteringProcess() {return fRayleighScattering;}
		G4OpMieHG* GetMieHGScatteringProcess() {return fMieHGScatteringProcess;}
		G4OpBoundaryProcess* GetBoundaryProcess() { return fBoundaryProcess;}

		void SetNbOfPhotonsCerenkov(G4int);

	private:

		G4OpWLS*             fWLSProcess;
		G4Cerenkov*          fCerenkovProcess;
		G4Scintillation*     fScintProcess;
		G4OpAbsorption*      fAbsorptionProcess;
		G4OpRayleigh*        fRayleighScattering;
		G4OpMieHG*           fMieHGScatteringProcess;
		G4OpBoundaryProcess* fBoundaryProcess;
 
    

};


#endif
