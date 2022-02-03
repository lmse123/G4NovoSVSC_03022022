/*
 * G4NovoSVSC:
 * Run action
 * Source file
 * Author: Kyrre Skjerdal (kyrre.skjerdal@hvl.no)
 */


/// \file NovoRunAction.cc
/// \brief Implementation of the NovoRunAction class

#include "NovoRunAction.hh"
#include "NovoPrimaryGeneratorAction.hh"
#include "NovoAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include <string>

NovoRunAction::NovoRunAction(NovoDetectorConstruction* det, NovoEventAction* ev)
:	G4UserRunAction(),
	fNevScint(0),
	fDetector(det),
	fEventAction(ev)
{
	//auto analysisManager = G4AnalysisManager::Instance();
	// Create analysis manager
	// The choice of analysis technology is done via selectin of a namespace
	// in NovoAnalysis.hh
	auto analysisManager = G4AnalysisManager::Instance();
	G4cout << "Using " << analysisManager->GetType() << G4endl;

	// Create directories
	analysisManager->SetHistoDirectoryName("histograms");
	analysisManager->SetNtupleDirectoryName("ntuples");
	analysisManager->SetVerboseLevel(1);

	// Create histograms
	analysisManager->CreateH1("nReflectionsHitPC", "Number of reflections", 400, 0, 400, "none", "n refl"); // H1: 0
	analysisManager->CreateH1("Timestamps1", "Timestamps of hits 1", 5000, 0.0*ns, 500.0*ns, "ns", "time"); // H1: 1
	analysisManager->CreateH1("Timestamps2", "Timestamps of hits 2", 5000, 0.0*ns, 500.0*ns, "ns", "time"); // H1: 2
	analysisManager->CreateH1("ReconPos", "Reconstructed position", 300, -150*mm, 150*mm, "mm", "pos"); // H1: 3
	analysisManager->CreateH1("AvgTime1", "Avg times for PC 1", 50, 0.0*ns, 100.0*ns, "ns", "time"); // H1: 4
	analysisManager->CreateH1("AvgTime2", "Avg times for PC 2", 50, 0.0*ns, 100.0*ns, "ns", "time"); // H1: 5
	analysisManager->CreateH1("PhotonCount1", "Photon count for PC 1", 2500, 0.0, 5000.0, "none", "n photon"); // H1: 6
	analysisManager->CreateH1("PhotonCount2", "Photon count for PC 2", 2500, 0.0, 5000.0, "none", "n photon"); // H1: 7
	analysisManager->CreateH1("CorrReconPos", "Corrected reconstructed position", 500, -250*mm, 250*mm, "mm", "pos"); // H1: 8
	analysisManager->CreateH1("ReconPosTime", "Reconstructed position", 300, -150*mm, 150*mm, "mm", "pos"); // H1: 9
	analysisManager->CreateH1("TimestampSum", "Sum time first photon", 100, 0.0*ns, 10.0*ns, "ns", "time"); // H1: 10
	analysisManager->CreateH1("EDep", "Energy deposited in scintilator", 100, 0.0*MeV, 5.0*MeV, "MeV", "E"); // H1: 11
	analysisManager->CreateH1("ESource", "Energy of source particle", 100, 0.0*MeV, 2.5*MeV, "MeV", "E"); // H1: 12
	analysisManager->CreateH1("SourceDirX", "Source X-direction", 100, -1, 1); // H1: 13
	analysisManager->CreateH1("SourceDirY", "Source Y-direction", 100, -1, 1); // H1: 14
	analysisManager->CreateH1("SourceDirZ", "Source Z-direction", 100, -1, 1); // H1: 15

	analysisManager->CreateH2("XYsource", "XY-dist of source", 100, -5.0*cm, 5.0*cm, 100, 5.0*cm, 25.0*cm, "cm", "cm"); //H2: 0
	analysisManager->CreateH2("XZsource", "XZ-dist of source", 100, -5.0*cm, 5.0*cm, 100, -5.0*cm, 5.0*cm, "cm", "cm"); //H2: 1
	analysisManager->CreateH2("YZsource", "YZ-dist of source", 100, 5.0*cm, 25.0*cm, 100, -5.0*cm, 5.0*cm, "cm", "cm"); //H2: 2

	G4int nScintX = fDetector->GetNScintX();
	G4int nScintY = fDetector->GetNScintY();
	G4int ntupleNo = 0;
	analysisManager->SetFirstNtupleId(0);
	analysisManager->SetFirstNtupleColumnId (0);
	for(G4int x = 0; x < nScintX; x++){
		for(G4int y = 0; y < nScintY; y++){
			G4String ntuplename = "Ntuple-" + std::to_string(x) + '-' + std::to_string(y);
			G4String ntupletitle = "Event info " + std::to_string(x) + ' ' + std::to_string(y);
			G4cout << "ntupleNo: " << ntupleNo << G4endl;
			analysisManager->CreateNtuple(ntuplename, ntupletitle);
			analysisManager->CreateNtupleDColumn(ntupleNo, "PhotonCountN"); //0
			analysisManager->CreateNtupleDColumn(ntupleNo, "PhotonCountP"); //1
			analysisManager->CreateNtupleDColumn(ntupleNo, "FirstTimestampN"); //2
			analysisManager->CreateNtupleDColumn(ntupleNo, "FirstTimestampP"); //3
			analysisManager->CreateNtupleDColumn(ntupleNo, "FifthTimestampN"); //4
			analysisManager->CreateNtupleDColumn(ntupleNo, "FifthTimestampP"); //5
			analysisManager->CreateNtupleDColumn(ntupleNo, "PosX_N"); //6
			analysisManager->CreateNtupleDColumn(ntupleNo, "PosX_P"); //7
			analysisManager->CreateNtupleDColumn(ntupleNo, "PosY_N"); //8
			analysisManager->CreateNtupleDColumn(ntupleNo, "PosY_P"); //9
			analysisManager->CreateNtupleDColumn(ntupleNo, "PosZ_N"); //10
			analysisManager->CreateNtupleDColumn(ntupleNo, "PosZ_P"); //11
			analysisManager->CreateNtupleDColumn(ntupleNo, "TimestampsVecN", fEventAction->GetTimestampsN()); // 12
			analysisManager->CreateNtupleDColumn(ntupleNo, "TimestampsVecP", fEventAction->GetTimestampsP()); //13
			analysisManager->FinishNtuple(ntupleNo); //?
			ntupleNo++;
		}
	}

	// My ntuple 17/12/2021
	G4int n_BSGating_tupleNo = 1; // <-- TODO: update if scint ntuples >1.
	analysisManager->CreateNtuple("BSGating", "Even info");
	analysisManager->CreateNtupleDColumn(n_BSGating_tupleNo, "nScintHits"); //0
	analysisManager->CreateNtupleDColumn(n_BSGating_tupleNo, "nBSGatingHits"); // 1
	analysisManager->CreateNtupleDColumn(n_BSGating_tupleNo, "ScintEDep"); // 2
	analysisManager->CreateNtupleDColumn(n_BSGating_tupleNo, "BSGatingEDep"); // 3
	analysisManager->CreateNtupleDColumn(n_BSGating_tupleNo, "Scint_isPrimary"); // 4
	analysisManager->CreateNtupleDColumn(n_BSGating_tupleNo, "BSGating_isPrimary"); // 5
	analysisManager->FinishNtuple(n_BSGating_tupleNo); //?
}



NovoRunAction::~NovoRunAction()
{
	delete G4AnalysisManager::Instance();
}

void NovoRunAction::BeginOfRunAction(const G4Run*)
{
	// inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);

	// Get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();

	// reset accumulables to their initial values
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Reset();

	// Open ouput file
	G4String filename = "NovoOutput";
	analysisManager->OpenFile(filename);

}

void NovoRunAction::EndOfRunAction(const G4Run* run)
{
	G4int nofEvents = run->GetNumberOfEvent();
	if ( nofEvents == 0) return;

	// Print
	//
	if (IsMaster()) {
		G4cout
			<< G4endl
			<< "--------------------End of Global Run-----------------------";
	}
	else {
		G4cout
			<< G4endl
			<< "--------------------End of Local Run------------------------";
	}
	G4cout
		<< G4endl
		<< " The run consists of " << nofEvents << " events"
		<< " where " << fNevScint.GetValue() << " detected scintillation photons in both photodetectors. "
		<< G4endl;

	// Get analysis manager
	auto analysisManager = G4AnalysisManager::Instance();
	analysisManager->Write();
	analysisManager->CloseFile();

}

void NovoRunAction::AddNevScint(G4int n)
{
	fNevScint += n;
}
