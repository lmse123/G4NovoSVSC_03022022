/*
 * G4NovoSVSC:
 * Event action: Actions happening at the beginning and end of each event.
 * Header file
 * Author: Kyrre Skjerdal (kyrre.skjerdal@hvl.no)
 */


/// \file NovoEventAction.hh
/// \brief Definition of the NovoEventAction class


#ifndef NovoEventAction_h
#define NovoEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4Accumulable.hh"
#include "G4ThreeVector.hh"
#include "NovoDetectorConstruction.hh"
#include <vector>

//class NovoDetectorConstruction;
class NovoRunAction;

/// Event action class
///

class NovoEventAction : public G4UserEventAction
{
	public:
		NovoEventAction(NovoDetectorConstruction* det);
		virtual ~NovoEventAction();
		virtual void BeginOfEventAction(const G4Event* event);
		virtual void EndOfEventAction(const G4Event* event);

		void SetEventVerbose(G4int v){fVerbose=v;}

		//void SetPCThreshold(G4int t){fPCThreshold=t;}

		void SetForceDrawPhotons(G4bool b){fForcedrawphotons=b;}
		void SetForceDrawNoPhotons(G4bool b){fForcenophotons=b;}


		inline G4int GetEventID() const {return fEventID;}
		std::vector<G4double>& GetTimestampsN() {return fTimestampsN;}
		std::vector<G4double>& GetTimestampsP() {return fTimestampsP;}

	private:
		//~ NovoDetectorConstruction* fDetector;
		//~ NovoRunAction* fRunAction;
		G4int fEventID;

		//G4int fSaveThreshold;

		G4int fScintCollID;
		G4int fPCCollID;
		G4int fBSGatingCollID;

		G4int fVerbose;
		std::vector<G4double> fTimestampsN;
		std::vector<G4double> fTimestampsP;

		//G4int fPCThreshold;

		G4bool fForcedrawphotons;
		G4bool fForcenophotons;
		G4ThreeVector CalculatePosition(G4ThreeVector vec);
		G4ThreeVector CalculatePositionSigmoid(G4ThreeVector vec);
		G4ThreeVector CalculatePositionTime(G4ThreeVector vec);
		G4ThreeVector CalculatePositionTimeSigmoid(G4ThreeVector vec);
		NovoDetectorConstruction *fDetector;
		//NovoRunAction* fRunAction;
};


#endif
