#include "G4ios.hh"
#include "GLG4OpAttenuation.hh"

#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4GeometryTolerance.hh"
#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Step.hh"

using namespace std;

/////////////////
// Hidden static variables and functions
/////////////////

// values in Cos2ThetaTable are used in equation to invert the equation
// for generating an angular distribution for scalar scattering:
//  dP/d(cos\theta) = (3/4) sin^2\theta
//   P_{cumulative} = (2+3*cos\theta-cos^3\theta)/4
//      cos(\theta) = 4*(P_{cumulative}-1/2) / (3-cos^2\theta)
//                                                ^^^^^^^^^^^ table value used
#define N_COSTHETA_ENTRIES 129
static double Cos2ThetaTable[N_COSTHETA_ENTRIES]; 
static int TableInitialized = 0;

GLG4DummyProcess GLG4OpAttenuation::fgAttenuation("Attenuation");
GLG4DummyProcess GLG4OpAttenuation::fgScattering("Scattering");


static void InitializeTable() {
  double angTolerance =
    G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  double cos2th = 0;

  for (int i=0; i<N_COSTHETA_ENTRIES-1; i++) {
    double x= i / (double)(N_COSTHETA_ENTRIES - 1);
    double old_cos2th;
    // Find root by iterating to convergence
    do {
      old_cos2th = cos2th;
      double costh = 2.0 * x / (3.0 - cos2th);
      cos2th = costh * costh;
    }
    while (fabs(old_cos2th - cos2th) > angTolerance);

    Cos2ThetaTable[i] = cos2th;
  }

  Cos2ThetaTable[N_COSTHETA_ENTRIES-1] = 1.0;
  TableInitialized = 1;
}

/////////////////
// Constructors and Destructor
/////////////////

GLG4OpAttenuation::GLG4OpAttenuation(const G4String& processName)
    : G4OpAbsorption(processName) {
  if (!TableInitialized)
    InitializeTable();
}


GLG4OpAttenuation::~GLG4OpAttenuation() {}

////////////
// Methods
////////////

G4VParticleChange*
GLG4OpAttenuation::PostStepDoIt(const G4Track& track, const G4Step& step) {
  aParticleChange.Initialize(track);

  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  G4double photonMomentum = particle->GetTotalMomentum();

  const G4Material* material = track.GetMaterial();
  G4MaterialPropertiesTable* materialPropertyTable =
    material->GetMaterialPropertiesTable();


  G4double opScatFrac = 0;

  if (materialPropertyTable) {
    std::stringstream scatteringPropertyName("OPSCATFRAC");

    // For multi-component, decide which one attenuates
    if (materialPropertyTable->ConstPropertyExists("NCOMPONENTS")) {
      G4double shortest = DBL_MAX;
      G4int compIndex = 0;
      G4int ncomp = (G4int) materialPropertyTable->GetConstProperty("NCOMPONENTS");
      for (G4int i=0; i<ncomp; i++) {
        G4double attLength = DBL_MAX;

        std::stringstream propname("");
        propname << "ABSLENGTH" << i;
        G4MaterialPropertyVector* attVector =
          materialPropertyTable->GetProperty(propname.str().c_str());

        if (attVector) {
          attLength = attVector->Value(photonMomentum);
        }

        G4double distance = -attLength * log(CLHEP::RandFlat::shoot());
        if (distance < shortest) {
          shortest = distance;
          compIndex = i;
        }
      }
      scatteringPropertyName << compIndex;

      // FIXME debugging output
      //G4cout << ">>> Attenuated by component " << compIndex << G4endl;
    }

    G4MaterialPropertyVector* opScatVector =
      materialPropertyTable->GetProperty(scatteringPropertyName.str().c_str());

    if (opScatVector) {
      opScatFrac = opScatVector->Value(photonMomentum);
    }
  }

  if (opScatFrac > 0 && G4UniformRand() < opScatFrac) {
    // Photon scattered coherently -- use scalar scattering (Rayleigh):
    // for fully polarized light, angular distribution \prop sin^2 \theta
    // where theta is angle between initial polarization and final
    // momentum vectors.  All light in Geant4 is fully polarized...
    G4double urand = G4UniformRand() - 0.5;
    G4double Cos2Theta0 =
      Cos2ThetaTable[(int)(fabs(urand)*2.0*(N_COSTHETA_ENTRIES-1)+0.5)];

    G4double CosTheta= 4.0 * urand / (3.0 - Cos2Theta0);

#ifdef G4DEBUG
    if (fabs(CosTheta)>1.0) {
      cerr << "GLG4OpAttenution: Warning, CosTheta=" << CosTheta
           << " urand=" << urand << endl;
      CosTheta = CosTheta > 0.0 ? 1.0 : -1.0;
    }
#endif

    G4double SinTheta = sqrt(1.0 - CosTheta * CosTheta);
    G4double Phi = (2.0 * G4UniformRand() - 1.0) * M_PI;
    G4ThreeVector e2(
      particle->GetMomentumDirection().cross(particle->GetPolarization()));

    G4ThreeVector NewMomentum =
      (CosTheta * particle->GetPolarization() +
       (SinTheta*cos(Phi)) * particle->GetMomentumDirection() +
       (SinTheta*sin(Phi)) * e2).unit();

    // Polarization is normal to new momentum and in same plane as
    // old new momentum and old polarization
    G4ThreeVector NewPolarization =
      (particle->GetPolarization() - CosTheta*NewMomentum).unit();

    aParticleChange.ProposeMomentumDirection(NewMomentum);
    aParticleChange.ProposePolarization(NewPolarization);

    postStepPoint->SetProcessDefinedStep(&fgScattering);
  }
  else {
    // Photon absorbed (may be re-radiated... but that is GLG4Scint's job)
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    postStepPoint->SetProcessDefinedStep(&fgAttenuation);

    if (verboseLevel > 0) {
      G4cout << "** Photon absorbed! **" << G4endl;
    }
  }

  return G4VDiscreteProcess::PostStepDoIt(track, step);
}

