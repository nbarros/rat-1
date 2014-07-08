//  "Attenuation" (absorption or scattering) of optical photons
//
//   GenericLAND Simulation
//
//   Original: Glenn Horton-Smith, Dec 2001
// 

#ifndef __GLG4OpAttenuation__
#define __GLG4OpAttenuation__

#include "G4OpAbsorption.hh"
#include "RAT/GLG4DummyProcess.hh"

class GLG4OpAttenuation : public G4OpAbsorption {
public:
  GLG4OpAttenuation(const G4String& processName="Attenuation");

  virtual ~GLG4OpAttenuation() {}

  // This is the method implementing attenuation of optical 
  // photons.  Fraction of photons scattered or absorbed is
  // determined by the MaterialProperyVector "OPSCATFRAC".
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

private:
  static GLG4DummyProcess fgAttenuation;
  static GLG4DummyProcess fgScattering;
};

#endif  // __GLG4OpAttenuation__

