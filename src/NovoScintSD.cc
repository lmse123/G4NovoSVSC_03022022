/// \file NovoScintSD.cc
/// \brief Implementation of the NovoScintSD class
//
//
#include "NovoScintSD.hh"
#include "NovoScintHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NovoScintSD::NovoScintSD(G4String name)
  : G4VSensitiveDetector(name)
{
	fScintCollection = NULL;
	collectionName.insert("scintCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

NovoScintSD::~NovoScintSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NovoScintSD::Initialize(G4HCofThisEvent* hitsCE)
{
	fScintCollection = new NovoScintHitsCollection(SensitiveDetectorName,collectionName[0]);
	//A way to keep all the hits of this event in one place if needed
	static G4int hitsCID = -1;
	if(hitsCID<0)
  {
		hitsCID = GetCollectionID(0);
	}
	hitsCE->AddHitsCollection( hitsCID, fScintCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool NovoScintSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ){
	// G4double edep = aStep->GetTotalEnergyDeposit();
	// if(edep==0.) return false; //No edep so dont count as hit
  //
	// G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
	// G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
	// G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();
  //
	// G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
  //
	// //Get the average position of the hit
	// G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
	// pos/=2.;
	// G4double globaltime = thePrePoint->GetGlobalTime();
	// NovoScintHit* scintHit = new NovoScintHit(thePrePV);
  //
	// scintHit->SetEdep(edep);
	// scintHit->SetPos(pos);
	// scintHit->SetTime(globaltime);
	// fScintCollection->insert(scintHit);
  //
	// return true;
	return false;
}

G4bool NovoScintSD::ProcessHits_constStep(const G4Step* aStep,G4TouchableHistory* ){

  G4int trackID = aStep->GetTrack()->GetTrackID();
  G4String particleSubType = aStep->GetTrack()->GetParticleDefinition()->GetParticleSubType(); // electron: e
  G4cout << particleSubType<< G4endl;
  // Hit information
  G4double edep = aStep->GetTotalEnergyDeposit(); // Hit energy
  G4int parentID = aStep->GetTrack()->GetParentID(); //ParentID

  // Filter hits
	//if(particleSubType!="e") {return false; }// ignore step if not electron
	//if(parentID!=1) {return false; }// Only save particles created by a primary (?)
	//if(parentID!=0) {return false; }// Only save primary hits
	if(edep==0.) return false; //No edep so dont count as hit

  G4bool isPrimary = (parentID == 0 ) ? true : false;

  //Get the average position of the hit
	G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
	G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
	G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
	pos/=2.;

	G4double globaltime = thePrePoint->GetGlobalTime();
  // Find pre step volume
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(thePrePoint->GetTouchable());
  G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();

  // Create hit and set attributes
  NovoScintHit* scintHit = new NovoScintHit(thePrePV);
	scintHit->SetEdep(edep);
	scintHit->SetPos(pos);
	scintHit->SetTime(globaltime);
	scintHit->SetIsPrimary(isPrimary);

  // Add hit to hit collection
	fScintCollection->insert(scintHit);

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NovoScintSD::EndOfEvent(G4HCofThisEvent* ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NovoScintSD::clear() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NovoScintSD::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void NovoScintSD::PrintAll() {}
