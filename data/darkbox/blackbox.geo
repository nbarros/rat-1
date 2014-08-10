{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "box",
size: [1000.0, 1000.0, 1000.0], // mm, half-length
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
r_max: 6000.0,
position: [0.0, 0.0, 0.0],
material: "scintillator",
color: [0.4, 0.4, 0.6, 0.05],
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

