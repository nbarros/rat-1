#ifndef __GLG4Scint__
#define __GLG4Scint__

#include <globals.hh>
#include <local_g4compat.hh>
#include <templates.hh>
#include <vector>
#include <G4ThreeVector.hh>
#include <G4ParticleMomentum.hh>
#include <G4Step.hh>
#include <G4OpticalPhoton.hh>
#include <G4DynamicParticle.hh>
#include <G4Material.hh>
#include <G4PhysicsTable.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4PhysicsOrderedFreeVector.hh>
#include <G4ParticleChange.hh>
#include <G4UImessenger.hh>
#include <G4VProcess.hh>
#include <G4hParametrisedLossModel.hh>
#include <TF1.h>
#include <TGraph.h>
#include <TSpline.h>
#include <RAT/GLG4DummyProcess.hh>

class G4UIcommand;
class G4UIdirectory;

/**
 * @class GLG4Scint
 *
 * Declares GLG4Scint class and helpers.
 * 
 * This file is part of the GenericLAND software library.
 *
 * @author Glenn Horton-Smith (Tohoku) 28-Jan-1999
 *
 * GLG4Scint is an extremely modified version of the G4Scintillation
 * process, so much so that it's not even a G4Process anymore!
 * Features include arbitrary scintillation light time profile and
 * spectra, Birks' law, particle-dependent specification of all
 * parameters, and reemission of optical photons killed by other processes.
 *
 *   - Has a GenericPostPostStepDoIt() function (note two "Post"s)
 *     instead of a PostStepDoIt() function. GenericPostPostStepDoIt()
 *     should be called by user in UserSteppingAction. This guarantees
 *     that GLG4Scint will absolutely be the last process considered, and
 *     will definitely see the energy loss by charged particles accurately.
 *
 *   - Modified to allow specification of absolute yield spectra,
 *     resolution scale, Birk's-law coefficient, and digitized waveform,
 *     customized for medium and (optionally) particle type.
 *
 *   - No longer calls G4MaterialPropertiesTable::GetProperty() in
 *     [Post]PostStepDoit() -- all needed data can be found quickly in
 *     the internal physics table.
 *
 *   - Uses poisson random distribution for number of photons if
 *     mean number of photons <= 12.
 *
 *   - The total scintillation yield is now found implicitly from
 *     the integral of the scintillation spectrum, which must now be
 *     in units of photons per photon energy.
 *
 *   - A scintillation yield can be defined and -if found- used instead of
 *     the implicit integral of the scintillation spectrum. This allows
 *     having scintillators with the same spectrum, but different light
 *     yields. (by Dario Motta)
 *
 *   - The materials property tables used are
 *       SCINTILLATION  == scintillation spectrum
 *       SCINTWAVEFORM  == scintillation waveform or time constants
 *       SCINTMOD       == resolution scale, Birk's constant, reference dE/dx
 *       REEMISSION     == reemission spectrum (scintillation for optical photons)
 *       REEMITWAVEFORM == reemission waveform or time constants
 *
 *   - SCINTILLATION is required in each scintillating medium.
 *     (Okay to omit if you don't want the medium to scintillate.)
 *
 *   - If SCINTWAVEFORM is missing, uses exponential waveform with default
 *     scintillation time. If SCINTWAVEFORM contains negative "momenta"
 *     then each "momentum" is the decay time and its corresponding value
 *     is the relative strength of that exponential decay.
 *     Otherwise, the "photon nnergy" of each element is a time, and the
 *     Value of each element is the relative strength.
 *
 *   - Default values of resolution scale (=1.0), Birk's constant (=0.0)
 *     and reference dE/dx (=0.0) are used if all or part of SCINTMOD is
 *     is missing.  SCINTMOD "photon energy" values should be set to the
 *     index number (0.0, 1.0, 2.0, with no units).
 *
 *   - Birk's law (see 1998 Particle Data Booklet eq. 25.1) is implemented
 *     as:
 *
 *       yield(dE/dx) =
 *           yield_ref * dE/dx * (1 + kb*(dE/dx)_ref) / (1 + kb*(dE/dx)).
 *
 *     I.e. the scintillation spectrum given in SCINTILLATION is
 *     measured for particles with dE/dx = (dE/dx)_ref. The usual
 *     formula is recovered if (dE/dx)_ref = 0.0 (the default).
 *     This is useful if you have an empirically-measured spectrum for
 *     some densely-ionizing particle (like an alpha).
 *
 *   - The constructor accepts an additional string argument, tablename,
 *     which allows selection of alternate property tables. For example,
 *     tablename = "neutron" might be used to allow specification of a
 *     different waveform for scintillation due to neutron energy deposition.
 *     The code then searches for tables with names of the form
 *
 *        "SCINTILLATIONneutron"
 *
 *     If it finds such a table, that table is used in preference to
 *     the default (un-suffixed) table when stepping particles of that type.
 *
 *   - The process generates at most sMaxTracksPerStep secondaries per step.
 *     If more "real" photons are needed, it increases the weight of the
 *     tracked opticalphotons. opticalphotons are thus macro-particles in
 *     the high-scintillation case. The code preserves an integer number
 *     of real photons per macro-particle.
 */
class GLG4Scint : public G4UImessenger {
public:
  class MyPhysicsTable;
  
  GLG4Scint(const G4String& tableName="", G4double lowerMassLimit=0.0);

  virtual ~GLG4Scint();

  // This routine is called for each step of any particle
  // in a scintillator.  For accurate energy deposition, must be called
  // from user-supplied UserSteppingAction, which also must stack
  // any particles created.  A pseudo-Poisson-distributed number of
  // photons is generated according to the scintillation yield formula,
  // distributed evenly along the track segment and uniformly into 4pi.
  G4VParticleChange* PostPostStepDoIt(const G4Track& aTrack,
                                      const G4Step& aStep);

  G4double GetLowerMassLimit() const {
    return fLowerMassLimit;
  }

  void DumpInfo() const;
  
  MyPhysicsTable* GetMyPhysicsTable() const {
    return fPhysicsTable;
  }

  // Set/get a command in the G4UImessenger
  void SetNewValue(G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue(G4UIcommand* command);

  /* Static interface */
  static G4VParticleChange* GenericPostPostStepDoIt(const G4Step *pStep);

  // For energy deposition diagnosis
  static void ResetTotEdep() {
    sTotEdep = sTotEdepQuenched = sTotEdepTime = 0.0;
    sScintCentroidSum *= 0.0;
  }

  static G4double GetTotEdep() { return sTotEdep; }

  static G4double GetTotEdepQuenched() { return sTotEdepQuenched; }

  static G4double GetTotEdepTime() { return sTotEdepTime; }

  static G4ThreeVector GetScintCentroid() {
    return sScintCentroidSum * (1.0 / sTotEdepQuenched);
  }

  static G4double GetQuenchingFactor() { return sQuenchingFactor; }
  static void SetQuenchingFactor(G4double qf=1.0) { sQuenchingFactor = qf; }

  static G4double GetMeanPhotonsPerSecondary() { return sMeanPhotonsPerSecondary; }
  static void SetMeanPhotonsPerSecondary(G4double mpps) { sMeanPhotonsPerSecondary = mpps; }

  static G4int GetMaxTracksPerStep() { return sMaxTracksPerStep; }
  static void SetMaxTracksPerStep(G4int mtps) { sMaxTracksPerStep = mtps; }

  static G4bool GetDoScintillation() { return sDoScintillation; }
  static void SetDoScintillation(G4bool flag) { sDoScintillation = flag; }

  static G4bool GetDoReemission() { return sDoReemission; }
  static void SetDoReemission(G4bool flag) { sDoReemission = flag; }

  static G4bool GetUseUserQF() { return sUseUserQF; }
  static void SetUseUserQF(G4bool flag) { sUseUserQF = flag; }

protected:
  // Pointer to the physics table which this instance of GLG4Scint will use.
  // You may create a separate instance of GLG4Scint for each particle, if
  // you like.
  MyPhysicsTable* fPhysicsTable;

  // Lower mass limit
  G4double fLowerMassLimit;

  // Return value of PostPostStepDoIt
  G4ParticleChange aParticleChange;
  
  // Vector of all existing GLG4Scint objects.
  // They register themselves when created, and remove themselves when deleted.
  static G4std::vector<GLG4Scint*> sScintProcessVector;

  // Top level of scintillation command
  static G4UIdirectory* sUIDirectory;

  // Maximum number of secondary tracks per step
  static G4int sMaxTracksPerStep;

  // Mean number of true photons per secondary track
  static G4double sMeanPhotonsPerSecondary;

  // Enable/disable the production of scintillation photons
  static G4bool sDoScintillation;

  // Enable/disable the production of reemission photons
  static G4bool sDoReemission;

  // Total "real" (unquenched) energy deposited
  static G4double sTotEdep;

  // Total quenched energy deposited
  static G4double sTotEdepQuenched;

  static G4double sTotEdepTime;

  static G4ThreeVector sScintCentroidSum;

  // Bogus processes used to tag photons created in GLG4Scint
  static GLG4DummyProcess sScintProcess;
  static GLG4DummyProcess sReemissionProcess;

  // Quenching factor
  static G4double sQuenchingFactor;

  // Enable/disable use of a user-provided constant quenching factor
  static G4bool sUseUserQF;
};


/**
 * @class GLG4Scint::MyPhysicsTable
 */
class GLG4Scint::MyPhysicsTable {
public:
  class Entry;

  static MyPhysicsTable* FindOrBuild(const G4String& name);

  static const MyPhysicsTable* GetDefault() { return sHead; }

  void IncUsedBy(void) { ++fUserCount; }

  void DecUsedBy(void) {
    if (--fUserCount <= 0) {
      delete this;
    }
  }

  const Entry* GetEntry(int i) const;

  void Dump() const;

  G4String GetName() const { return fName; }
  void SetName(G4String name) { fName = name; }

private:
  friend class GLG4Scint;

  MyPhysicsTable();

  ~MyPhysicsTable();

  void Build(const G4String& name);

  MyPhysicsTable* fNext;
  G4int fUserCount;
  std::vector<Entry> fData;
  G4String fName;

  static MyPhysicsTable* sHead;
};


/**
 * @class GLG4Scint::MyPhysicsTable::Entry
 */
class GLG4Scint::MyPhysicsTable::Entry {
public:
  Entry() { Initialize(); }

  ~Entry() { Destroy(); }

  // Initialize this object to default values
  void Initialize();

  // Delete pointer members if we own them
  void Destroy();

  // Initialize this object using material properties from the provided
  // table
  void Build(const G4String& name, int materialIndex,
             G4MaterialPropertiesTable* matProps);

  // Parse a material property vector with scintillation timing into
  // an ordered vector, supporting both time-series data and decay
  // time constants.
  static G4PhysicsOrderedFreeVector*
  BuildTimeIntegral(G4MaterialPropertyVector* waveformData);

  G4MaterialPropertyVector* fQuenchingArray;
  std::vector<G4MaterialPropertyVector*> fReemissionProbVector;
  std::vector<G4PhysicsOrderedFreeVector*> fReemissionSpectrumIntegral;
  std::vector<G4PhysicsOrderedFreeVector*> fReemissionTimeIntegral;
  G4PhysicsOrderedFreeVector* fScintillationSpectrumIntegral;
  G4PhysicsOrderedFreeVector* fScintillationTimeIntegral;
  G4double fResolutionScale;
  G4double fBirksConstant;
  G4double fRef_dEdx;
  G4double fLightYield;
  bool fOwnSpectrumIntegral;
  bool fOwnTimeIntegral;
};


// Inline method definitions

inline void GLG4Scint::DumpInfo() const {
  G4cout << "GLG4Scint[" << *(fPhysicsTable->GetName()) << "] {" << G4endl
         << "  fLowerMassLimit=" << fLowerMassLimit << G4endl;
#ifdef G4DEBUG
  if (fPhysicsTable)
    fPhysicsTable->Dump();
  }
#endif
  G4cout << "}" << G4endl;
}


inline const GLG4Scint::MyPhysicsTable::Entry*
GLG4Scint::MyPhysicsTable::GetEntry(int i) const {
  return &(fData.at(i));
}

#endif  // __GLG4Scint__

