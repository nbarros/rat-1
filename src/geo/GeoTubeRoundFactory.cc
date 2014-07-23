#include <RAT/GeoTubeRoundFactory.hh>

#include <G4UnionSolid.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>

using namespace std;
namespace RAT {


GeoTubeRoundFactory::~GeoTubeRoundFactory() {

}

G4VSolid* GeoTubeRoundFactory::ConstructSolid(DBLinkPtr table) {
	string volumeName = table->GetIndex();
	// Parameters for the tube
  G4double tub_r_max = table->GetD("tube_r_max") * mm;
  G4double tub_size_z = table->GetD("tube_size_z") * mm;

  // Optional parameters for the tube
  G4double tub_r_min = 0.0;
  try { tub_r_min = table->GetD("tube_r_min") * mm; }
  catch (DBNotFoundError &e) { };
  G4double tub_phi_start = 0.0;
  try { tub_phi_start = table->GetD("tube_phi_start") * deg; }
  catch (DBNotFoundError &e) { };
  G4double tub_phi_delta = twopi;
  try { tub_phi_delta = table->GetD("tube_phi_delta") * deg; }
  catch (DBNotFoundError &e) { };

  G4VSolid* baseTube = new G4Tubs(volumeName + "_tub",
  				     tub_r_min,tub_r_max , tub_size_z,
  				     tub_phi_start,tub_phi_delta);

  // Now get the sphere part
  // The optional parameters defaults are chosen to match a seamless
  // connection with the tube part

  // If not specified, it uses the same radius as the tube
  G4double sph_r_max = tub_r_max;
  try{ sph_r_max = table->GetD("sphere_r_max") * mm; }
  catch (DBNotFoundError &e) { };
  G4double sph_r_min = tub_r_min;
  try { sph_r_min = table->GetD("sphere_r_min") * mm; }
  catch (DBNotFoundError &e) { };

  // Offset moves the center of the sphere along the axis of the tube
  G4double sph_c_offset = 0.0;
  try { sph_c_offset = table->GetD("sphere_center_offset") * mm; }
  catch (DBNotFoundError &e) { };

  G4double sph_phi_start = 0.0;
  try { sph_phi_start = table->GetD("sphere_phi_start") * deg; }
  catch (DBNotFoundError &e) { };
  G4double sph_phi_delta = pi;
  try { sph_phi_start = table->GetD("sphere_phi_delta") * deg; }
  catch (DBNotFoundError &e) { };
  G4double sph_theta_start = 0.0;
  try { sph_theta_start = table->GetD("sphere_theta_start") * deg; }
  catch (DBNotFoundError &e) { };
  G4double sph_theta_delta = pi;
  try { sph_theta_delta = table->GetD("sphere_theta_delta") * deg; }
  catch (DBNotFoundError &e) { };

  // Construct the temporary sphere
  G4VSolid* baseSphere = new G4Sphere(volumeName + "_sph",
  																		sph_r_min, sph_r_max,
  																		sph_phi_start, sph_phi_delta,
  																		sph_theta_start,sph_theta_delta);

  /// Create the rotations. By default Geant4 builds the tube standing along Z
  // put the inter on the x-axis, then rotate
  G4RotationMatrix* rotate = new G4RotationMatrix();
  //rotate->rotateZ(90 * deg);
  //rotate->rotateY(90 * deg);
  rotate->rotateX(-90 * deg);
  // apply the offset along z to the top of the tube
  G4ThreeVector translate(0.0,0.0,tub_size_z);


  // Finally build the union
  G4VSolid* finalSolid = new G4UnionSolid(volumeName,baseTube,baseSphere,rotate,translate);

  return finalSolid;
}

} /* namespace RAT */
