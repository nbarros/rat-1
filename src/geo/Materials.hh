#ifndef __RAT_Materials__
#define __RAT_Materials__

#include <RAT/DB.hh>
#include <RAT/Log.hh>
#include <map>
#include <string>

class G4Material;
class G4PhysicsOrderedFreeVector;
class G4OpticalSurface;
typedef G4PhysicsOrderedFreeVector G4MaterialPropertyVector;

namespace RAT {

class Materials {
public:
  // Load all materials into memory
  static void LoadMaterials();

  // GEANT4 has no global store of optical surface information,
  // so we need to keep it here
  static void LoadOpticalSurfaces();

  static void ConstructMaterials();
  static void ReadPropertyTable();

  static std::map<std::string, G4OpticalSurface*> optical_surface;

private:
  inline static bool BuildMaterial(std::string name, DBLinkPtr ptr);

  // Load all entries from the OPTICS tables
  static void LoadOptics();

  // Create a MaterialPropertiesTable from an OPTICS table
  static void BuildMaterialPropertiesTable(G4Material* material, DBLinkPtr table);

  // Load a single material property vector from an OPTICS table by name
  // This handles energy vs. wavelength basis, etc.
  static G4MaterialPropertyVector* LoadProperty(DBLinkPtr table, std::string name);
};

}  // namespace RAT

#endif  // __RAT_Materials__

