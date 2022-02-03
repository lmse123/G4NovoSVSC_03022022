/*
 * G4Rangeverifcation:
 * Primary generator action
 * Source file
 * Author: Kyrre Skjerdal (kyrre.skjerdal@hvl.no)
 */

/// \file NovoPrimaryGeneratorAction.hh
/// \brief Implementation of the NovoPrimaryGeneratorAction class


#include "NovoPrimaryGeneratorAction.hh"
#include "NovoDetectorConstruction.hh"
#include "NovoAnalysis.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4VAnalysisManager.hh"

NovoPrimaryGeneratorAction::NovoPrimaryGeneratorAction(NovoDetectorConstruction* det)
:G4VUserPrimaryGeneratorAction(),
 fParticleGun(0),
 fDetector(det)
{
  fParticleGun = new G4GeneralParticleSource();
}

NovoPrimaryGeneratorAction::~NovoPrimaryGeneratorAction()
{
  delete fParticleGun;
}

void NovoPrimaryGeneratorAction::SetDefaultPrimaryParticle()
{

}

void NovoPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fParticleGun->GeneratePrimaryVertex(anEvent);
  fEbeamCumul += fParticleGun->GetParticleEnergy();
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillH1(12, fParticleGun->GetParticleEnergy());
  G4ThreeVector pos = fParticleGun->GetParticlePosition();
  analysisManager->FillH2(0, pos.x(), pos.y());
  analysisManager->FillH2(1, pos.x(), pos.z());
  analysisManager->FillH2(2, pos.y(), pos.z());
  G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  analysisManager->FillH1(13, dir.x());
  analysisManager->FillH1(14, dir.y());
  analysisManager->FillH1(15, dir.z());
}
