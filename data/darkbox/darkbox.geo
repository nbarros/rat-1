{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "box",
size: [1000.0, 1000.0, 500.0], // mm, half-length
material: "air",
invisible: 1,
}

{
name: "GEO",
index: "glass",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world",
type: "sphere",
r_max: 240.0,
position: [0.0, 0.0, 0.0],
material: "glass",
color: [0.4, 0.4, 0.6, 0.05],
}

{
name: "GEO",
index: "target",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "glass",
type: "sphere",
r_max: 235.0,
position: [0.0, 0.0, 0.0],
material: "te_0p3_labppo_scintillator_bisMSB_Dec2013",
//material: "myscintillator",
color: [1.0, 0.9, 0.7, 0.35],
}

{ 
name: "GEO", 
index: "pmts", 
valid_begin: [0, 0], 
valid_end: [0, 0], 
mother: "world", 
type: "pmtarray", 
pmt_type: "r580", 
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner", 
pos_table: "PMTINFO", 
orientation: "point",
orient_point: [0.0, 0.0, 0.0], 
} 

{
name: "GEO",
index: "collimator",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world",
type: "tube",
r_max: 9.0,
r_min: 5.0,
size_z: 27.5,
position: [295.0, 0.0, 0.0],
orientation: [-1.0,0.0,0.0],
material: "pvc",
color: [0.0, 0.0, 0.0, 0.55],
}


{
name: "GEO",
index: "LED",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world",
type: "tuberound",
tube_r_max: 5.0,
tube_r_min: 0.0,
tube_size_z: 1.5,
position: [321.5, 0.0, 0.0],
orientation: [-1.0,0.0,0.0],
material: "acrylic_polycast",
color: [1.0, 0.0, 0.0, 0.0],
//color: [0.4, 0.4, 0.6, 0.55],
}
