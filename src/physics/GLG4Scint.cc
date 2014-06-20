/**
 * @file GLG4Scint.cc
 *
 * For GLG4Scint class, providing advanced scintillation process.
 * Distantly based on an extensively modified version of G4Scintillation.cc.
 * 
 * This file is part of the GenericLAND software library.
 *
 * @author Glenn Horton-Smith (Tohoku) 28-Jan-1999
 */

#include <sstream>
#include <Randomize.hh>
#include <G4UnitsTable.hh>
#include <G4ios.hh>
#include <G4Timer.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIdirectory.hh>
#include <G4TrackFastVector.hh>
#include <G4IonTable.hh>
#include <G4hZiegler1985Nuclear.hh>
#include <G4hZiegler1985p.hh>
#include <G4hIonEffChargeSquare.hh>
#include <G4hParametrisedLossModel.hh>
#include <G4PSTARStopping.hh>
#include <G4AtomicShells.hh>
#include <G4ParticleTable.hh>
#include <G4Event.hh>
#include <G4EventManager.hh>
#include <TF1.h>
#include <RAT/PhotonThinning.hh>
#include <RAT/TrackInfo.hh>
#include <RAT/Sampling.hh>
#include <RAT/EventInfo.hh>
#include <RAT/Log.hh>
#include <GLG4Scint.hh>

G4std::vector<GLG4Scint*> GLG4Scint::sScintProcessVector;
G4UIdirectory* GLG4Scint::sUIDirectory = NULL;
GLG4DummyProcess GLG4Scint::sScintProcess("Scintillation", fUserDefined);
GLG4DummyProcess GLG4Scint::sReemissionProcess("Reemission", fUserDefined);

G4bool GLG4Scint::sDoScintillation = true;
G4bool GLG4Scint::sDoReemission = true;
G4int GLG4Scint::sMaxTracksPerStep = 180000;
G4double GLG4Scint::sMeanPhotonsPerSecondary = 1.0; 
G4bool GLG4Scint::sUseUserQF = false;
G4double GLG4Scint::sQuenchingFactor = 1.0;

G4double GLG4Scint::sTotEdep = 0.0;
G4double GLG4Scint::sTotEdepQuenched = 0.0;
G4double GLG4Scint::sTotEdepTime = 0.0;
G4ThreeVector GLG4Scint::sScintCentroidSum(0.0, 0.0, 0.0);


GLG4Scint::GLG4Scint(const G4String& tablename, G4double lowerMassLimit)
    : fLowerMassLimit(lowerMassLimit) {
  fPhysicsTable = MyPhysicsTable::FindOrBuild(tablename);
  fPhysicsTable->IncUsedBy();

#ifdef G4DEBUG
  fPhysicsTable->Dump();
#endif

  // Register this instance in the global ordered list
  if (sScintProcessVector.empty() ||
      fLowerMassLimit >= sScintProcessVector.back()->fLowerMassLimit) {
    sScintProcessVector.push_back(this);
  }
  else {
    G4std::vector<GLG4Scint*>::iterator it;
    for (it=sScintProcessVector.begin();
         it!=sScintProcessVector.end();
         it++) {
      if (fLowerMassLimit < (*it)->fLowerMassLimit) {
         sScintProcessVector.insert(it, this);
         break;
      }
    }
  }

  // Create UI commands if necessary
  if (sUIDirectory == NULL) {
    sUIDirectory = new G4UIdirectory("/glg4scint/");
    sUIDirectory->SetGuidance("Scintillation process control");

    G4UIcommand* cmd;
    cmd = new G4UIcommand("/glg4scint/on", this);
    cmd->SetGuidance("Turn on scintillation");

    cmd = new G4UIcommand("/glg4scint/off", this);
    cmd->SetGuidance("Turn off scintillation");

    cmd = new G4UIcommand("/glg4scint/reemission", this);
    cmd->SetGuidance("Turn on/off reemission of absorbed opticalphotons");
    cmd->SetParameter(new G4UIparameter("status", 's', false));

    cmd = new G4UIcommand("/glg4scint/maxTracksPerStep", this);
    cmd->SetGuidance("Set maximum number of opticalphoton tracks per step\n"
                     "(If more real photons are needed, "
                     "weight of tracked particles is increased.)\n" );
    cmd->SetParameter(new G4UIparameter("maxTracksPerStep", 'i', false));

    cmd = new G4UIcommand("/glg4scint/meanPhotonsPerSecondary", this);
    cmd->SetGuidance("Set mean number of \"real\" photons per secondary\n");
    cmd->SetParameter(new G4UIparameter("meanPhotonsPerSecondary",
                                        'd', false));

    cmd = new G4UIcommand("/glg4scint/dump", this);
    cmd->SetGuidance("Dump tables");

    cmd = new G4UIcommand("/glg4scint/setQF", this);
    cmd->SetGuidance("Set a constant quenching factor, default is 1");
    cmd->SetParameter(new G4UIparameter("QuenchingFactor", 'd', false));
  }
  
  RAT::debug << "GLG4Scint[" << tablename << "]" << " created" << newline;
}


GLG4Scint::~GLG4Scint()  {
  fPhysicsTable->DecUsedBy();
  G4std::vector<GLG4Scint*>::iterator it;
  for (it=sScintProcessVector.begin();
       it!=sScintProcessVector.end();
       it++) {
    if (*it == this) {
      sScintProcessVector.erase(it);
      break;
    }
  }
}


G4VParticleChange*
GLG4Scint::PostPostStepDoIt(const G4Track& track, const G4Step& step) {
  // Note: Below this anonymous namespace is label PostStepDoIt_DONE
  {
    // Prepare to generate an event, organizing to
    // check for things that cause an early exit.
    aParticleChange.Initialize(track);

    const G4Material* material = track.GetMaterial();

    const MyPhysicsTable::Entry* physicsEntry =
      fPhysicsTable->GetEntry(material->GetIndex());

    bool flagReemission = false;

    if (track.GetDefinition() == G4OpticalPhoton::OpticalPhoton()) {
      // Add the original parent track and creation step to PhotonIDParentStep
      G4Event* event =
        G4EventManager::GetEventManager()->GetNonconstCurrentEvent();

      RAT::EventInfo* eventInfo =
        dynamic_cast<RAT::EventInfo*>(event->GetUserInformation());

      RAT::TrackInfo* currentTrackInfo =
        dynamic_cast<RAT::TrackInfo*>(track.GetUserInformation());

      if (eventInfo && eventInfo->StorePhotonIDs) {
        // Only occurs on first step
        if (track.GetCurrentStepNumber() == 1) {
          eventInfo->PhotonIDParentStep[track.GetTrackID()].push_back(track.GetParentID());
          eventInfo->PhotonIDParentStep[track.GetTrackID()].push_back(currentTrackInfo->GetCreatorStep()-1);
        }
      }

     flagReemission =
       (GetDoReemission() &&
        track.GetTrackStatus() == fStopAndKill &&
        step.GetPostStepPoint()->GetStepStatus() != fGeomBoundary);

      // If not reemission, there is nothing left to do for optical photons
      if (!flagReemission) {
        goto PostStepDoIt_DONE;
      }
    }

    G4double totalEnergyDeposit = step.GetTotalEnergyDeposit();
  
    if (totalEnergyDeposit <= 0.0 && !flagReemission) {
      goto PostStepDoIt_DONE;
    }

    if (!physicsEntry) {
      goto PostStepDoIt_DONE;
    }

    // Finds E-dependent QF, unless the user provided an E-independent one
    if (!GetUseUserQF() && physicsEntry->fQuenchingArray != NULL) {
      // This interpolates or uses first/last value if out of range
      SetQuenchingFactor(physicsEntry->fQuenchingArray->Value(track.GetVertexKineticEnergy()));
    }
    else {
      SetQuenchingFactor(1.0);
    }

    // Retrieve the light yield or scintillation integral for this material  
    G4double scintillationYield = physicsEntry->fLightYield;

    G4PhysicsOrderedFreeVector* scintillationIntegral =
      physicsEntry->fScintillationSpectrumIntegral;

    G4PhysicsOrderedFreeVector* reemissionIntegral =
      physicsEntry->fReemissionSpectrumIntegral;
    
    if (!scintillationIntegral && !reemissionIntegral) {
      goto PostStepDoIt_DONE;
    }

    // If no LY defined Max Scintillation Integral == ScintillationYield
    if (!scintillationYield) {
      if (!scintillationIntegral) {
        scintillationYield = 0;
      }
      else {
        scintillationYield = scintillationIntegral->GetMaxValue();
      }
    }

    // Set positions, directions, etc.
    G4StepPoint* preStepPoint = step.GetPreStepPoint();
    G4StepPoint* postStepPoint = step.GetPostStepPoint();

    G4ThreeVector x0 = preStepPoint->GetPosition();
    G4ThreeVector p0 = preStepPoint->GetMomentumDirection();
    G4double t0 = preStepPoint->GetGlobalTime();

    // Figure out how many photons we want to make
    G4int numSecondaries;
    G4double weight;
    //G4double reemissionProb = 0;
    //G4int numComponents = -1;
    //G4int absorberIndex = -1;

    if (flagReemission) {
      // Generate reemission photons
      G4MaterialPropertiesTable* mptScint =
        material->GetMaterialPropertiesTable();

      // Check if there are multiple components
      //if (mpt_scint->ConstPropertyExists("COMPONENTS")) {
      //  RAT::Log::Die("GLG4Scint: COMPONENTS not yet implemented");
      //  numComponents = mpt_scint->GetConstProperty("NUM_COMP");
      //}

      G4MaterialPropertyVector* reemissionProbVector =
        mptScint->GetProperty("REEMISSION_PROB");

      if (!reemissionProbVector) {
        goto PostStepDoIt_DONE;
      }

      G4double reemissionProb =
        reemissionProbVector->Value(track.GetKineticEnergy());

      numSecondaries = (G4int)(CLHEP::RandPoisson::shoot(reemissionProb));

      if (numSecondaries == 0) {
        goto PostStepDoIt_DONE;
      }

      weight = track.GetWeight();
    }
    else {
      // Generate primary scintillation photons

      // Apply Birks' law
      G4double birksConstant = physicsEntry->fBirksConstant;
      G4double quenchedTotalEnergyDeposit = totalEnergyDeposit;
      if (birksConstant != 0) {
        G4double dEdx = totalEnergyDeposit / step.GetStepLength();
        quenchedTotalEnergyDeposit /= (1.0 + birksConstant * dEdx);
      }

      // Track total edep, quenched edep
      sTotEdep += totalEnergyDeposit;
      sTotEdepQuenched += quenchedTotalEnergyDeposit;
      sTotEdepTime = t0;
      sScintCentroidSum +=
        quenchedTotalEnergyDeposit * (x0 + p0 * (0.5 * step.GetStepLength()));

      // Now we are done if we are not actually making photons here
      if (!GetDoScintillation()) {
        goto PostStepDoIt_DONE;
      }

      // Calculate MeanNumPhotons
      G4double meanNumPhotons =
        (scintillationYield *
         GetQuenchingFactor() *
         quenchedTotalEnergyDeposit *
         (1.0 + birksConstant * (physicsEntry->fRef_dEdx)));

      if (meanNumPhotons <= 0) {
        goto PostStepDoIt_DONE;
      }
      
      // Randomize number of TRACKS (not photons).
      //
      // This gets statistics right for number of PE after applying
      // boolean random choice to final absorbed track (change from
      // old method of applying binomial random choice to final absorbed
      // track, which did want poissonian number of photons divided
      // as evenly as possible into tracks)
      // Note for weight=1, there's no difference between tracks and photons.
      G4double meanNumTracks =
        (meanNumPhotons /
         GetMeanPhotonsPerSecondary() /
         RAT::PhotonThinning::GetFactor());

      G4double resolutionScale = physicsEntry->fResolutionScale;
      if (meanNumTracks > 12.0) {
        numSecondaries=
          (G4int)(CLHEP::RandGauss::shoot(meanNumTracks,
                                          resolutionScale
                                          * sqrt(meanNumTracks)));
      }
      else {
        if (resolutionScale > 1.0) {
          G4double scale = sqrt(resolutionScale * resolutionScale - 1.0);
          meanNumTracks =
            CLHEP::RandGauss::shoot(meanNumTracks, scale * meanNumTracks);
        }
        numSecondaries =
          (G4int)(CLHEP::RandPoisson::shoot(meanNumTracks));
      }

      weight = GetMeanPhotonsPerSecondary();
      if (numSecondaries > GetMaxTracksPerStep()) {
        // It's probably better to just set meanPhotonsPerSecondary to
        // a big number if you want a small number of secondaries, but
        // this feature is retained for backwards compatibility.
        weight = weight * numSecondaries / GetMaxTracksPerStep();
        numSecondaries = GetMaxTracksPerStep();
      }
    }
    
    // If there are no photons produced, return unchanged particle and no
    // secondaries.
    if (numSecondaries <= 0) {
      aParticleChange.SetNumberOfSecondaries(0);
      goto PostStepDoIt_DONE;
    }

    // Making at least one secondary... notify the proper authorities.
    aParticleChange.SetNumberOfSecondaries(numSecondaries);

    if (!flagReemission) {
      if (track.GetTrackStatus() == fAlive) {
        aParticleChange.ProposeTrackStatus(fSuspend);
      }
    }

    // Now look up waveform information we need to add the secondaries
    G4PhysicsOrderedFreeVector* waveformIntegral;
    if (flagReemission) {
      waveformIntegral = physicsEntry->fReemissionTimeIntegral;
    }
    else {
      waveformIntegral = physicsEntry->fScintillationTimeIntegral;
    }

    for (G4int iSecondary=0; iSecondary<numSecondaries; iSecondary++) {
      // Determine photon momentum
      G4double sampledMomentum;

      if (!flagReemission && scintillationIntegral != NULL) {
        // Normal scintillation
        G4double CIIvalue =
          G4UniformRand() * scintillationIntegral->GetMaxValue();
        sampledMomentum = scintillationIntegral->GetEnergy(CIIvalue);
      }
      else {
        // Reemission
        G4double CIIvalue =
          (G4UniformRand() *
           reemissionIntegral->Value(track.GetKineticEnergy()));

        if (CIIvalue == 0.0) {
          // Return unchanged particle and no secondaries  
          aParticleChange.SetNumberOfSecondaries(0);
          goto PostStepDoIt_DONE;
        }

        sampledMomentum = reemissionIntegral->GetEnergy(CIIvalue);
        aParticleChange.ProposeLocalEnergyDeposit(track.GetKineticEnergy() -
                                                  sampledMomentum);

        if (sampledMomentum > track.GetKineticEnergy()) {
           goto PostStepDoIt_DONE;
        }
      }

      // Generate random photon direction
      G4double cost = 1.0 - 2.0 * G4UniformRand();
      G4double sint = sqrt(1.0 - cost * cost);  // Fix G4Scint bug

      G4double phi = 2 * M_PI * G4UniformRand();
      G4double sinp = sin(phi);
      G4double cosp = cos(phi);

      G4double px = sint * cosp;
      G4double py = sint * sinp;
      G4double pz = cost;

      // Create photon momentum direction vector
      G4ParticleMomentum photonMomentum(px, py, pz);

      // Determine polarization of new photon
      G4double sx = cost * cosp;
      G4double sy = cost * sinp; 
      G4double sz = -sint;

      G4ThreeVector photonPolarization(sx, sy, sz);

      G4ThreeVector perp = photonMomentum.cross(photonPolarization);

      phi = 2 * M_PI * G4UniformRand();
      sinp = sin(phi);
      cosp = cos(phi);

      photonPolarization = cosp * photonPolarization + sinp * perp;
      photonPolarization = photonPolarization.unit();

      // Generate a new photon
      G4DynamicParticle* scintillationPhoton =
        new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), 
                              photonMomentum);

      scintillationPhoton->SetPolarization(photonPolarization.x(),
                                           photonPolarization.y(),
                                           photonPolarization.z());

      scintillationPhoton->SetKineticEnergy(sampledMomentum);

      // Generate new G4Track object
      G4ThreeVector secondaryPosition;
      G4double deltaTime;
      if (flagReemission) {
        deltaTime = postStepPoint->GetGlobalTime() - t0;
        secondaryPosition = postStepPoint->GetPosition();
      }
      else {
        G4double delta = G4UniformRand() * step.GetStepLength();
        secondaryPosition = x0 + delta * p0;

        // Start deltaTime based on where on the track it happened
        deltaTime = (delta /
                     ((preStepPoint->GetVelocity() +
                       postStepPoint->GetVelocity()) / 2.0));
      }

      // Delay for scintillation time
      if (waveformIntegral) {
        G4double wfSample = G4UniformRand() * waveformIntegral->GetMaxValue();
        G4double sampledDelayTime = waveformIntegral->GetEnergy(wfSample);
        deltaTime += sampledDelayTime;
      }

      // Set secondary time
      G4double secondaryTime = t0 + deltaTime;

      // Create secondary track                
      G4Track* secondaryTrack =
        new G4Track(scintillationPhoton,
                    secondaryTime,
                    secondaryPosition);

      secondaryTrack->SetWeight(weight);
      secondaryTrack->SetParentID(track.GetTrackID());
      RAT::TrackInfo* trackInfo = new RAT::TrackInfo();
      
      trackInfo->SetCreatorStep(track.GetCurrentStepNumber());

      if (flagReemission) {
        secondaryTrack->SetCreatorProcess(&sReemissionProcess);
        trackInfo->SetCreatorProcess(sReemissionProcess.GetProcessName());
      }
      else {
        secondaryTrack->SetCreatorProcess(&sScintProcess);
        trackInfo->SetCreatorProcess(sScintProcess.GetProcessName());
      }

      secondaryTrack->SetUserInformation(trackInfo);

      // Add the secondary to the ParticleChange object
      aParticleChange.SetSecondaryWeightByProcess(true); // recommended
      aParticleChange.AddSecondary(secondaryTrack);
      
      // AddSecondary() overrides our setting of the secondary track weight
      // in Geant4.3.1 & earlier (and also later, at least until Geant4.7?).
      // Maybe not required if SetWeightByProcess(true) called,
      // but we do both, just to be sure
      secondaryTrack->SetWeight(weight);
    }
  }

PostStepDoIt_DONE:
#ifdef G4VERBOSE
  RAT::debug << "\n Exiting from GLG4Scint::DoIt -- NumberOfSecondaries = " 
             << aParticleChange.GetNumberOfSecondaries() << RAT::newline;
#endif

  return &aParticleChange;
}


G4VParticleChange* GLG4Scint::GenericPostPostStepDoIt(const G4Step *step) {
  G4Track* track = step->GetTrack();
  G4double mass = track->GetDynamicParticle()->GetMass();
  G4std::vector<GLG4Scint*>::iterator it = sScintProcessVector.begin();
  for (int i=sScintProcessVector.size(); (i--)>1;) {
    it++;
    if (mass < (*it)->fLowerMassLimit) {
      return (*(--it))->PostPostStepDoIt(*track, *step);
    }
  }
  return (*it)->PostPostStepDoIt(*track, *step);
}


/////////////////////////////////////
// GLG4Scint::MyPhysicsTable
/////////////////////////////////////

GLG4Scint::MyPhysicsTable* GLG4Scint::MyPhysicsTable::sHead = NULL;


GLG4Scint::MyPhysicsTable::MyPhysicsTable() : fNext(NULL), fUserCount(0) {}


GLG4Scint::MyPhysicsTable::~MyPhysicsTable() {
  if (fUserCount != 0) {
    G4cerr << "Error, GLG4Scint::MyPhysicsTable is being deleted with "
           << "user count = " << fUserCount << G4endl;
    return;
  }
}


void GLG4Scint::MyPhysicsTable::Dump() const {
  G4cout << " GLG4Scint::MyPhysicsTable {"  << G4endl
         << "  fName=" << fName << G4endl
         << "  fData.size()=" << fData.size() << G4endl
         << "  fUserCount=" << fUserCount << G4endl;

  for (size_t i=0; i<fData.size(); i++) {
    G4cout << "  fData[" << i << "]= { // "
           << (*G4Material::GetMaterialTable())[i]->GetName() << G4endl;
    G4cout << "   fScintillationSpectrumIntegral=";
    if (fData[i].fScintillationSpectrumIntegral)
      (fData[i].fScintillationSpectrumIntegral)->DumpValues();
    else
      G4cout << "NULL" << G4endl;
    
    G4cout << "   fReemissionSpectrumIntegral=";
    if (fData[i].fReemissionSpectrumIntegral)
      (fData[i].fReemissionSpectrumIntegral)->DumpValues();
    else
      G4cout << "NULL" << G4endl;
      
    G4cout << "   fScintillationTimeIntegral=";
    if (fData[i].fScintillationTimeIntegral)
      (fData[i].fScintillationTimeIntegral)->DumpValues();
    else
      G4cout << "NULL" << G4endl;

    G4cout << "   fReemissionTimeIntegral=";
    if (fData[i].fReemissionTimeIntegral)
      (fData[i].fReemissionTimeIntegral)->DumpValues();
    else
      G4cout << "NULL" << G4endl;

    G4cout << "   resolutionScale=" << fData[i].fResolutionScale
           << "   birksConstant=" << fData[i].fBirksConstant
           << "   ref_dE_dx=" << fData[i].fRef_dEdx << G4endl
           << "   light yield=" << fData[i].fLightYield << G4endl;

    G4cout << "Quenching = " << G4endl;
    if (fData[i].fQuenchingArray != NULL)
      fData[i].fQuenchingArray->DumpValues();
    else
      G4cout << "NULL" << G4endl << "  }" << G4endl;
  }

  G4cout << " }" << G4endl;
}


GLG4Scint::MyPhysicsTable*
GLG4Scint::MyPhysicsTable::FindOrBuild(const G4String& name) {
  // Head should always exist and should always be the default (name=="")
  if (sHead == NULL) {
    sHead = new MyPhysicsTable;
    sHead->Build("");
  }

  MyPhysicsTable* pTable = sHead;
  while (pTable) {
    if (name == *(pTable->fName))
      return pTable;
    pTable = pTable->fNext;
  }

  pTable = new MyPhysicsTable;
  pTable->Build(name);
  pTable->fNext = sHead->fNext;  // Always keep head pointing to default
  sHead->fNext = pTable;

  return pTable;
}


void GLG4Scint::MyPhysicsTable::Build(const G4String& newname) {
  fData.clear();

  // Name in the physics list, i.e. "" or "heavy" or "alpha" etc.
  // This is a suffix on material property vectors in RATDB
  fName = "";

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  // Create new physics tables
  fData.resize(G4Material::GetNumberOfMaterials());
  for (size_t i=0; i<fData.size(); i++) {
    const G4Material* material = (*materialTable)[i];
    fData[i].Build(fName, i, material->GetMaterialPropertiesTable());
  }
}

/////////////////////////////////////
// GLG4Scint::MyPhysicsTable::Entry
/////////////////////////////////////

void GLG4Scint::MyPhysicsTable::Entry::Initialize() {
  fQuenchingArray = NULL;
  fScintillationSpectrumIntegral = NULL;
  fReemissionSpectrumIntegral = NULL;
  fScintillationTimeIntegral = NULL;
  fReemissionTimeIntegral = NULL;
  fOwnSpectrumIntegral = false;
  fOwnTimeIntegral = false;
  fResolutionScale = 1.0;
  fBirksConstant = 0.0;
  fRef_dEdx = 0.0;
  fLightYield = 0.0;
}

void GLG4Scint::MyPhysicsTable::Entry::Destroy() {
  if (fOwnSpectrumIntegral) {
   delete fScintillationSpectrumIntegral;
   delete fReemissionSpectrumIntegral;
  }

  if (fOwnTimeIntegral) {
    delete fScintillationTimeIntegral;
    delete fReemissionTimeIntegral;
  }

  delete fQuenchingArray;
}


G4PhysicsOrderedFreeVector*
GLG4Scint::MyPhysicsTable::Entry::BuildTimeIntegral(
    G4MaterialPropertyVector* waveformData) {
  G4PhysicsOrderedFreeVector* integral = NULL;

  // Do we have time-series or decay-time data?
  if (waveformData->GetMinLowEdgeEnergy() >= 0.0) {
    // We have digitized waveform (time-series) data
    integral = RAT::Integrate_MPV_to_POFV(waveformData);
  }
  else {
    // We have decay-time data. Sanity-check user's values:
    if (waveformData->Energy(waveformData->GetVectorLength() - 1) > 0.0) {
      G4cerr << "GLG4Scint::MyPhysicsTable::Entry::Build():  "
             << "SCINTWAVEFORM "
             << " has both positive and negative X values.  "
             << " Undefined results will ensue!" << G4endl;
    }

    G4double maxtime = -3.0 * (waveformData->GetMinLowEdgeEnergy());
    G4double mintime = -1.0 * (waveformData->GetMaxLowEdgeEnergy());
    G4double binwidth = mintime / 100;
    size_t nbins = (size_t) (maxtime / binwidth) + 1;
    G4double* tval = new G4double[nbins];
    G4double* ival = new G4double[nbins];
    for (size_t i=0; i<nbins; i++) {
      tval[i] = i * maxtime / nbins;
      ival[i] = 0.0;
    }

    for (unsigned int j=0; j<waveformData->GetVectorLength(); j++) {
      G4double ampl = (*waveformData)[j];
      G4double decy = waveformData->Energy(j);
      for (size_t i=0; i<nbins; i++) {
        ival[i] += ampl * (1.0 - exp(tval[i] / decy));
      }
    }

    for (size_t i=0; i<nbins; i++) {
      ival[i] /= ival[nbins-1];
    }

    integral = new G4PhysicsOrderedFreeVector(tval, ival, nbins);

    // G4PhysicsOrderedFreeVector makes its own copy of any array passed
    // to its constructor
    delete[] tval;
    delete[] ival;
  }

  return integral;
}


void
GLG4Scint::MyPhysicsTable::Entry::Build(const G4String& name,
                                        int materialIndex,
                                        G4MaterialPropertiesTable* matProps) {
  // Reset to defaults
  Destroy();
  Initialize();

  // Exit, leaving default values, if no material properties
  if (!matProps) {
    return;
  }
  
  // Retrieve vector of scintillation andreemission wavelength intensity
  // for the material from the optical properties table
  std::stringstream propname;

  propname.str("");
  propname << "SCINTILLATION" << name;
  G4MaterialPropertyVector* scintillationSpectrumVector = 
    matProps->GetProperty(propname.str().c_str());

  propname.str("");
  propname << "REEMISSION" << name;
  G4MaterialPropertyVector* reemissionSpectrumVector =
    matProps->GetProperty(propname.str().c_str());

  if (scintillationSpectrumVector && !reemissionSpectrumVector) {
    G4cout << "Warning! Found a scintillator without Re-emission spectrum";
    G4cout << " (probably a scintillator without WLS)" << G4endl;
    G4cout << "I will assume that for this material this spectrum is equal ";
    G4cout << "to the primary scintillation spectrum." << G4endl;
    reemissionSpectrumVector = scintillationSpectrumVector;
  }
     
  if (scintillationSpectrumVector) {
    if (matProps->ConstPropertyExists("LIGHT_YIELD"))
      fLightYield=matProps->GetConstProperty("LIGHT_YIELD");
    else { 
      G4cout << "Warning! Found a scintillator without LIGHT_YIELD parameter.";
      G4cout << "\nI will assume that for this material this parameter is ";
      G4cout << "implicit in the scintillation integral." << G4endl;

      // If no light yield, it's no scintillator
      scintillationSpectrumVector = NULL;
    }

    // Find the integral
    if (!scintillationSpectrumVector) {
      fScintillationSpectrumIntegral = NULL;
    }
    else {
      fScintillationSpectrumIntegral =
        RAT::Integrate_MPV_to_POFV(scintillationSpectrumVector);
    }

    fReemissionSpectrumIntegral =
      RAT::Integrate_MPV_to_POFV(reemissionSpectrumVector);   
    fOwnSpectrumIntegral = true;
  }
  else {
    // Use default integral (possibly null)
    fScintillationSpectrumIntegral =
      MyPhysicsTable::GetDefault()->GetEntry(materialIndex)->fScintillationSpectrumIntegral;
    fReemissionSpectrumIntegral = fScintillationSpectrumIntegral;
    fOwnSpectrumIntegral = false;
  }

  // Retrieve vector of scintillation and reemission time profiles
  // for the material from the optical properties table
  propname.str("");
  propname << "SCINTWAVEFORM" << name;
  G4MaterialPropertyVector* scintillationWaveform = 
    matProps->GetProperty(propname.str().c_str());

  propname.str("");
  propname << "REEMITWAVEFORM" << name;
  G4MaterialPropertyVector* reemissionWaveform = 
    matProps->GetProperty(propname.str().c_str());

  if (scintillationWaveform && !reemissionWaveform) {
    G4cout << "GLG4Scint: Warning: Using the primary scintillation timing "
           << "for reemission" << G4endl;
    reemissionWaveform = scintillationWaveform;
  }

  if (scintillationWaveform) {
    // User-specified scintillation time profile
    fScintillationTimeIntegral = BuildTimeIntegral(scintillationWaveform);
    fReemissionTimeIntegral = BuildTimeIntegral(reemissionWaveform);
    fOwnTimeIntegral = true;
  }
  else {
    // Use the default
    fScintillationTimeIntegral =
      MyPhysicsTable::GetDefault()->GetEntry(materialIndex)->fScintillationTimeIntegral;
    fReemissionTimeIntegral = fScintillationTimeIntegral;
    fOwnTimeIntegral = false;
  }

  // Retrieve vector of scintillation "modifications"
  // for the material from the material's optical
  // properties table ("SCINTMOD")    
  propname.str("");
  propname << "SCINTMOD" << name;
  G4MaterialPropertyVector* scintModVector = 
    matProps->GetProperty(propname.str().c_str());

  // Parse the entries in SCINTMOD:
  //
  //   0 - ResolutionScale
  //   1 - BirksConstant
  //   2 - Ref_dE_dx
  //
  // Use default if no particle-specific value is given
  if (!scintModVector) {
    scintModVector = matProps->GetProperty("SCINTMOD");
  }
  if (scintModVector) {
    for (unsigned int i=0; i<scintModVector->GetVectorLength(); i++) {
      G4double key = scintModVector->Energy(i);
      G4double value = (*scintModVector)[i];

      if (key == 0) {
        fResolutionScale = value;
      }
      else if (key == 1) {
        fBirksConstant = value;
      }
      else if (key == 2) {
        fRef_dEdx = value;
      }
      else {
        G4cerr << "GLG4Scint::MyPhysicsTable::Entry::Build: "
               << "Warning, unknown key " << key
               << "in SCINTMOD" << name << G4endl;
      }
    }
  }

  // Retrieve any energy-dependent quenching factor array
  propname.str("");
  propname << "QF" << name;
  fQuenchingArray =
    matProps->GetProperty(propname.str().c_str());
}


void GLG4Scint::SetNewValue(G4UIcommand* command, G4String newValues) {
  G4String commandName= command -> GetCommandName();
  if (commandName == "on") {
    SetDoScintillation(true);
  }
  else if (commandName == "off") {
    SetDoScintillation(false);
  }
  else if (commandName == "reemission") {
    char* endptr;
    G4int i = strtol((const char*)newValues, &endptr, 0);
    if (*endptr != '\0') { // non-numerical argument
      if (!(i = strcmp((const char*)newValues, "on"))) {
        SetDoReemission(true);
      }
      else if (!(i = strcmp((const char*)newValues, "off"))) {
        SetDoReemission(false);
      }
      else {
        G4cerr << "Command /glg4scint/reemission given unknown parameter "
               << '\"' << newValues << '\"' << G4endl
               << "  old value unchanged: "
               << ( GetDoReemission() ? "on" : "off" ) << G4endl;
      }
    }
    else {
      SetDoReemission(i != 0);
    }
  }
  else if (commandName == "maxTracksPerStep") {
    G4int i = strtol((const char*)newValues, NULL, 0);
    if (i > 0) {
      SetMaxTracksPerStep(i);
    }
    else {
      G4cerr << "Value must be greater than 0, old value unchanged" << G4endl;
    }
  }
  else if (commandName == "meanPhotonsPerSecondary") {
    G4double d = strtod((const char*)newValues, NULL);
    if (d >= 1.0) {
      SetMeanPhotonsPerSecondary(d);
    }
    else {
      G4cerr << "Value must be >= 1.0, old value unchanged" << G4endl;
    }
  }
  else if (commandName == "dump") {
    G4std::vector<GLG4Scint *>::iterator it = sScintProcessVector.begin();
    for (; it != sScintProcessVector.end(); it++) {
      (*it)->DumpInfo();
    }
  }
  else if (commandName == "setQF") {
    G4double d = strtod((const char *)newValues, NULL);
    if (d <= 1.0) {
      SetQuenchingFactor(d);
      SetUseUserQF(true);
    }
    else {
      G4cerr << "The quenching factor is <= 1.0, old value unchanged" << G4endl;
    }
  }
  else {
    G4cerr << "No GLG4Scint command named " << commandName << G4endl;
  }
}


G4String GLG4Scint::GetCurrentValue(G4UIcommand* command) {
  G4String commandName = command->GetCommandName();
  if (commandName == "on" || commandName == "off") {
    return GetDoScintillation() ? "on" : "off";
  }
  else if (commandName == "reemission") {
    return GetDoReemission() ? "1" : "0";
  }
  else if (commandName == "maxTracksPerStep") {
    char outbuff[64];
    snprintf(outbuff, 64, "%d", GetMaxTracksPerStep());
    return G4String(outbuff);
  }
  else if (commandName == "meanPhotonsPerSecondary") {
    char outbuff[64];
    snprintf(outbuff, 64, "%g", GetMeanPhotonsPerSecondary());
    return G4String(outbuff);
  }
  else if (commandName == "dump") {
    return "/glg4scint/dump not supported";
  }
  else if(commandName=="setQF"){
    char outbuff[64];
    snprintf(outbuff, 64, "%g", GetQuenchingFactor());
    return G4String(outbuff);
  }
  else {
    return (commandName + " is not a valid GLG4Scint command");
  }
}

