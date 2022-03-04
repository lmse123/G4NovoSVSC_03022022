/*
 * G4NovoSVSC:
 * Event action
 * Source file
 * Author: Kyrre Skjerdal (kyrre.skjerdal@hvl.no)
 */


/// \file NovoEventAction.cc
/// \brief Implementation of the NovoEventAction class


#include "NovoEventAction.hh"
#include "NovoRunAction.hh"
#include "NovoAnalysis.hh"
#include "NovoBSGatingHit.hh"
#include "NovoScintHit.hh"
#include "NovoPhotocatHit.hh"
#include "NovoUserEventInformation.hh"

#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4VAnalysisManager.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"

#include <algorithm>
#include <vector>

NovoEventAction::NovoEventAction(NovoDetectorConstruction* det)
:	G4UserEventAction(),
	fEventID(0),
	fScintCollID(-1),
	fPCCollID(-1),
	fBSGatingCollID(-1),
	fVerbose(1),
	fTimestampsN(0),
	fTimestampsP(0),
	fForcedrawphotons(false),
	fForcenophotons(false),
	fDetector(det),
	fScintHitPosXs(0), //set to 0?
	fScintHitPosYs(0), //set to 0?
	fScintHitPosZs(0) //set to 0?
	//fRunAction(runAction)
{
}


NovoEventAction::~NovoEventAction()
{
}

void NovoEventAction::BeginOfEventAction(const G4Event* anEvent)
{
	G4EventManager::
    GetEventManager()->SetUserInformation(new NovoUserEventInformation);

	fScintHitPosXs.clear();
	fScintHitPosYs.clear();
	fScintHitPosZs.clear();
	fTimestampsN.clear();
	fTimestampsP.clear();
  fVerbose = 0;
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	if(fBSGatingCollID<0){
		fBSGatingCollID=SDman->GetCollectionID("BSGatingCollection");
	}
	if(fScintCollID<0){
		fScintCollID=SDman->GetCollectionID("scintCollection");
	}
	if(fPCCollID<0){
		fPCCollID=SDman->GetCollectionID("photocatHitCollection");
	}
	// her kan du lage en ntuple for lagrind av info. create ntuple, double
/*	G4int evNum = anEvent->GetEventID();
	//G4cout << "Event number: " << evNum << G4endl;
	if(evNum < 100){
		G4cout << "Event number: " << evNum << G4endl;
	} else if (evNum < 1000 && evNum%100 == 0){
		G4cout << "Event number: " << evNum << G4endl;
	} else if (evNum < 10000 && evNum%1000 == 0){
		G4cout << "Event number: " << evNum << G4endl;
	} else if(evNum%10000 == 0){
		G4cout << "Event number: " << evNum << G4endl;
	}
*/
}

void NovoEventAction::EndOfEventAction(const G4Event* anEvent)
{
	NovoUserEventInformation* eventInformation = (NovoUserEventInformation*)anEvent->GetUserInformation();

	NovoBSGatingHitsCollection* bsHC = 0;
	NovoScintHitsCollection* scintHC = 0;
	NovoPhotocatHitsCollection* pcHC = 0;
	G4HCofThisEvent* hitsCE = anEvent->GetHCofThisEvent();
	auto analysisManager = G4AnalysisManager::Instance();

	// Get hit collections of sensitive detectors (SDs)
	if(hitsCE){
		if(fBSGatingCollID >= 0){
			bsHC = (NovoBSGatingHitsCollection*)(hitsCE->GetHC(fBSGatingCollID));
		}
		if(fScintCollID >= 0){
			scintHC = (NovoScintHitsCollection*)(hitsCE->GetHC(fScintCollID));
		}
		if(fPCCollID >=0 ){
			pcHC = (NovoPhotocatHitsCollection*)(hitsCE->GetHC(fPCCollID));
		}
	}

// Process the event's scintillator hit collection
	if(scintHC){
		int n_hit = scintHC->entries();
		G4ThreeVector  eWeightPos(0.);
		G4double edep;
		G4double edepMax=0;

		for(int i=0;i<n_hit;i++){ //gather info on hits in scintillator
			edep=(*scintHC)[i]->GetEdep();
			eventInformation->IncEDep(edep); //sum up the edep
			eWeightPos += (*scintHC)[i]->GetPos()*edep;//calculate energy weighted pos
			if(edep>edepMax){
				edepMax=edep;//store max energy deposit
				G4ThreeVector posMax=(*scintHC)[i]->GetPos();
				eventInformation->SetPosMax(posMax,edep);
			}
		}
		if(eventInformation->GetEDep()==0.){
			if(fVerbose>0)G4cout<<"No hits in the scintillator this event."<<G4endl;
		}
		else{
			//Finish calculation of energy weighted position
			eWeightPos/=eventInformation->GetEDep();
			eventInformation->SetEWeightPos(eWeightPos);
			if(fVerbose>0){
				G4cout << "\tEnergy weighted position of hits in Scintillator : " << eWeightPos/mm << G4endl;
				G4cout << "\tTotal energy deposition in scintillator : " << eventInformation->GetEDep() / MeV << " (MeV)" << G4endl;
			}
			analysisManager->FillH1(11, eventInformation->GetEDep());
			//~ analysisManager->FillNtupleDColumn(1, 16, eventInformation->GetEDep());
		}

	}

	// Process the event's BS gating hit collection
	if(bsHC){
		G4int n_scintHit = scintHC->entries();
		G4int n_BSGatingHit = bsHC->entries();
		G4double ScintEDep = eventInformation->GetEDep(); // tot energy dep in scint. (0 by default)
		G4double BSGatingEDep = 0; // tot energy dep in BS detector
		G4double edep = 0;
		G4bool p = 0;

		int BSGating_isPrimary = 0;
		for(int i=0;i<n_BSGatingHit;i++){ //gather info on hits in CeBr3
			// Calculate total edep in BS gating detector
			edep=(*bsHC)[i]->GetEdep();
			p=(*bsHC)[i]->GetIsPrimary();
			BSGatingEDep = BSGatingEDep+edep;
			//G4cout << "i: "<<i<<" ; isP: "<< p<<" ;edep: "<< BSGatingEDep << G4endl;

			NovoBSGatingHit * BSGatingHit = static_cast<NovoBSGatingHit*>(bsHC->GetHit(i));
			//G4cout<< BSGatingHit->GetIsPrimary()<<G4endl;
			if (BSGatingHit->GetIsPrimary()== true) BSGating_isPrimary++; //count primaries in BS Gating collection
		}

		int Scint_isPrimary = 0;
		//std::vector<G4ThreeVector> test;
		//G4ThreeVector hitPos; // hit coordinate vector
		for (int i = 0; i < n_scintHit; i++){

			NovoScintHit * scintHit = static_cast<NovoScintHit*>(scintHC->GetHit(i));
			G4double hitPosX=(*scintHC)[i]->GetPos().x();
			G4double hitPosY=(*scintHC)[i]->GetPos().y();
			G4double hitPosZ=(*scintHC)[i]->GetPos().z();
			fScintHitPosXs.push_back(hitPosX);
			fScintHitPosYs.push_back(hitPosY);
			fScintHitPosZs.push_back(hitPosZ);

			G4cout<< fScintHitPosXs[i]<<":"<< fScintHitPosYs[i]<<":"<< fScintHitPosZs[i]<<":" <<G4endl;
			//G4cout<< scintHit->GetIsPrimary()<<G4endl;
			if (scintHit->GetIsPrimary()== true) Scint_isPrimary++; //count primaries in scint collection
		}

		//if(n_BSGatingHit>0) G4cout << "gating hits: " <<n_BSGatingHit<< "," <<n_scintHit<<G4endl;
		//G4cout << "primary Scint hits: "<<Scint_isPrimary << G4endl;
		//G4cout << "primary BS hits: "<<BSGating_isPrimary << G4endl;

		G4int n_BSGating_tupleNo = 1; // <-- TODO: update if scint ntuples >1.
		analysisManager->FillNtupleDColumn(n_BSGating_tupleNo,0,n_scintHit);
		analysisManager->FillNtupleDColumn(n_BSGating_tupleNo,1,n_BSGatingHit);
		analysisManager->FillNtupleDColumn(n_BSGating_tupleNo,2,ScintEDep);
		analysisManager->FillNtupleDColumn(n_BSGating_tupleNo,3,BSGatingEDep);
		analysisManager->FillNtupleDColumn(n_BSGating_tupleNo,4,Scint_isPrimary);
		analysisManager->FillNtupleDColumn(n_BSGating_tupleNo,5,BSGating_isPrimary);
		// TODO: scint position. No this is already stored
		analysisManager->AddNtupleRow(n_BSGating_tupleNo);
	}

	// Process the event's photocathode hit collection
	if(pcHC){
		G4int nHitPcs = pcHC->entries();
		if(nHitPcs > 0){
			//NovoRunAction::AddNevScint(1);
			for(G4int i = 0; i < nHitPcs; i++){
				G4double photoncount = (*pcHC)[i]->GetPhotonCount();
				G4int pcnumber = (*pcHC)[i]->GetPhotocatNumber();
				G4int ntupleNo = pcnumber/2;

				G4ThreeVector pos = (*pcHC)[i]->GetPhotocatPos();
				std::vector<G4double> timestamps = (*pcHC)[i]->GetTimestamps();
				sort(timestamps.begin(), timestamps.end());
				if(pcnumber%2 == 0){
					// N
					fTimestampsN = timestamps;// (N)egative z-axis size
					analysisManager->FillNtupleDColumn(ntupleNo, 0, photoncount);
					analysisManager->FillH1(6, photoncount);
					analysisManager->FillNtupleDColumn(ntupleNo, 2, fTimestampsN[0]);
					analysisManager->FillNtupleDColumn(ntupleNo, 4, fTimestampsN[4]);
					analysisManager->FillNtupleDColumn(ntupleNo, 6, pos.x());
					analysisManager->FillNtupleDColumn(ntupleNo, 8, pos.y());
					analysisManager->FillNtupleDColumn(ntupleNo, 10, pos.z());
					//~ G4cout << "pcnumber: " << pcnumber << " "<<"ntupleNo: " << ntupleNo << G4endl;
					//~ G4cout << "pos: " << pos << G4endl;
					//~ G4cout << "photoncount: " << photoncount << G4endl;

				}
				else if(pcnumber%2 != 0){
					// P
					fTimestampsP = timestamps; // (P)ositive z-axis size
					analysisManager->FillNtupleDColumn(ntupleNo, 1, photoncount);
					analysisManager->FillH1(7, photoncount);
					analysisManager->FillNtupleDColumn(ntupleNo, 3, fTimestampsP[0]);
					analysisManager->FillNtupleDColumn(ntupleNo, 5, fTimestampsP[4]);
					analysisManager->FillNtupleDColumn(ntupleNo, 7, pos.x());
					analysisManager->FillNtupleDColumn(ntupleNo, 9, pos.y());
					analysisManager->FillNtupleDColumn(ntupleNo, 11, pos.z());
					//~ G4cout << "pcnumber: " << pcnumber << " "<<"ntupleNo: " << ntupleNo << G4endl;
					//~ G4cout << "pos: " << pos << G4endl;
					//~ G4cout << "photoncount: " << photoncount << G4endl;
				}
			}
		}
	}
	// ??
	for(G4int i = 0; i < fDetector->GetNScint(); i++){
		analysisManager->AddNtupleRow(i);
	}
}


G4ThreeVector NovoEventAction::CalculatePosition(G4ThreeVector vec)
{
	G4ThreeVector pos(0., 0., 0.);
	if (fDetector->GetScintZ() == 30.0*cm){
		G4double p0 = 1.02500e-02;
		G4double p1 = 3.38410e+00;
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(p0 + vec.z()*p1);
	}
	if (fDetector->GetScintZ() == 20.0*cm){
		G4double p0 = -6.18240e-02;
		G4double p1 = 7.26294e+00;
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(p0 + vec.z()*p1);
	}
	return pos;
}

G4ThreeVector NovoEventAction::CalculatePositionSigmoid(G4ThreeVector vec)
{
	G4ThreeVector pos(0., 0., 0.);
	//~ G4cout << "fDetector->GetScintZ(): " << fDetector->GetScintZ()/cm << G4endl;
	if (fDetector->GetScintZ() == 30.0*cm){
		G4double p0 = 4.53401e+00;
		G4double p1 = 1.01916e-02;
		G4double z = (p0*vec.z())/(1+std::abs(p1*vec.z()));
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(z);
	}
	if (fDetector->GetScintZ() == 20.0*cm){
		G4double p0 = 1.36367e+01;
		G4double p1 = -8.50906e-02;
		G4double z = (p0*vec.z())/(1+std::abs(p1*vec.z()));
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(z);
	}
	return pos;
}

G4ThreeVector NovoEventAction::CalculatePositionTime(G4ThreeVector vec)
{
	G4ThreeVector pos(0., 0., 0.);
	if (fDetector->GetScintZ() == 30.0*cm){
		G4double p0 = -1.21310e-02;
		G4double p1 = 2.31182e+00;
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(p0 + vec.z()*p1);
	}
	if (fDetector->GetScintZ() == 20.0*cm){
		G4double p0 = 6.45101e-02;
		G4double p1 = 4.68424e+00;
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(p0 + vec.z()*p1);
	}
	return pos;
}

G4ThreeVector NovoEventAction::CalculatePositionTimeSigmoid(G4ThreeVector vec)
{
	G4ThreeVector pos(0., 0., 0.);
	//~ G4cout << "fDetector->GetScintZ(): " << fDetector->GetScintZ()/cm << G4endl;
	if (fDetector->GetScintZ() == 30.0*cm){
		G4double p0 = 2.89218e+00;
		G4double p1 = 5.09484e-03;
		G4double z = (p0*vec.z())/(1+std::abs(p1*vec.z()));
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(z);
	}
	if (fDetector->GetScintZ() == 20.0*cm){
		G4double p0 = 7.97598e+00;
		G4double p1 = 4.41081e-02;
		G4double z = (p0*vec.z())/(1+std::abs(p1*vec.z()));
		pos.setX(vec.x());
		pos.setY(vec.y());
		pos.setZ(z);
	}
	return pos;
}
