Component Optics
----------------
Defining Composite Materials
````````````````````````````
RAT has a ''multi-component'' optical model which allows you to construct
composite materials which combine the optics of other materials. To use this
feature, add an entry to the `MATERIALS` DB table with the material/elemental
composition, and add an entry to the `OPTICS` table with this format::

    {
    name: "OPTICS",
    index: "my_composite_material",
    COMPONENTS: ["water", "scintillator"],
    COMPONENT_FRACTIONS: [0.25, 0.75],

    // Define the refractive index explicitly for the composite
    RINDEX_option: "wavelength",
    RINDEX_value1: [200d, 800d],
    RINDEX_value2: [1.5d, 1.5d],

    PROPERTY_LIST: ["COMPONENTS", "COMPONENT_FRACTIONS", "RINDEX"],
    }

The following optical properties, if present, will be pulled in from all
components::

    OPSCATFRAC, ABSLENGTH, SCINTILLATION, SCINTWAVEFORM, SCINTMOD,
    REEMISSION, REEMISSION_PROB, REEMITWAVEFORM

Only one component may have primary scintillation properties (`SCINTILLATION`,
`SCINTWAVEFORM`, `SCINTMOD`) defined. Other properties are copied into the
Geant4 material properties for the composite material with clever names like::

    OPSCATFRAC0      Scattering fraction for component 0 (water)
    OPSCATFRAC1      Scattering fraction for component 1 (scintillator)
    REEMITWAVEFORM1  Reemission timing for component 1 (scintillator)

Finally, the attenuation lengths are combined into one overall length, which
is used by Geant4 when determining whether the composite material attenuates.
If it does, we go back to component-specific properties to decide how to
proceed.

*Implementation notes:* See `src/geo/Materials.cc`.

Absorption and Scattering
`````````````````````````
If Geant4 decides, based on the overall attenuation length of the composite
material, that we ''attenuate'' an optical photon, a random number is thrown
based on the component fractions to decide which component was the culprit.

We then look up the scattering fraction (`OPSCATFRAC`) for that component and
choose whether to Rayleigh scatter the photon or to absorb it. If it is
absorbed in a scintillator, there is a chance that it will be reemitted
in the next section.

*Implementation notes:* See `src/physics/GLG4OpAttenuation.cc`. Note that the
index of the absorbing component is recorded in the `RAT::EventInfo`
associated with the `G4Event`, in order to communicate this to `GLG4Scint`.

Scintillation and Reemission
````````````````````````````
For the details of the scintillation physics, see :doc:`processes`. The only
difference with composite materials is that we use the properties of the
absorbing component to determine the reemission spectrum and time profile
of any reemitted photons.

*Implementation notes:* See `src/physics/GLG4Scint.cc`. The heavy lifting is
done on initialization in `GLG4Scint::MyPhysicsTable::Entry::Build`: here we
load the reemission spectra and time profiles, and just look up the right one
by index when dealing with an absorbed photon.

