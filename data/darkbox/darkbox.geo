{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "box",
size: [2000.0, 2000.0, 2000.0], // mm, half-length
material: "air",
color: [0.4, 0.4, 0.6, 0.05],
invisible: 0,
}

{
name: "GEO",
index: "acrylic_box",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world",
type: "box",
size: [50.0, 50.0, 50.0],
position: [0.0, 0.0, 0.0],
material: "acrylic_uvt_good",
color: [0.4, 0.4, 0.6, 0.10],
invisible:0
}

{
name: "GEO",
index: "target",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "acrylic_box",
type: "box",
size: [45.0, 45.0, 45.0],
position: [0.0, 0.0, 0.0],
material: "mystery_scintillator",
color: [0.0, 1.0, 1.0,0.5],
invisible:0
}

{ 
name: "GEO", 
index: "pmts", 
valid_begin: [0, 0], 
valid_end: [0, 0], 
mother: "world", 
type: "pmtarray", 
pmt_type: "r11780_hqe", 
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner", 
pos_table: "PMTINFO",
orientation: "manual",
invisible:0
} 
