SetFactory("OpenCASCADE");

lc = 5e3;
lc_fault = 100; // Fine mesh in the fault zone //Change this to adjust mesh size

Fault_length = 22e3; //fault length
Fault_width = 15e3; //fault width
Fault_thickness = 1000; //fault thickness

HighDamage_thickness = 200;

// Nucleation in X,Z local coordinates
X_nucl = -7e3;
R_nucl = 2.5e3;

Xmax = 30e3;
Xmin = -Xmax;

Ymax = 0;
Ymin = -30e3;

//
Width_nucl = Ymin/2;

Zmin =  -30e3;
Zmax =   30e3;

Box(1) = {Xmin, 0, Zmin, 2*Xmax, Ymin, 2*Zmax};

// Create a damage zone
Box(2) = {-Fault_length/2, -Fault_width, -Fault_thickness/2, Fault_length, Fault_width, Fault_thickness};

// Create a nucleation patch
//Box(3) = {X_nucl-R_nucl/2, Width_nucl-R_nucl/2, -thickness_nucl/2, R_nucl, R_nucl, thickness_nucl};
Box(3) = {X_nucl-R_nucl/2, -Fault_width/2-R_nucl/2, -HighDamage_thickness/2, R_nucl, R_nucl, HighDamage_thickness};

// Create a cdbm allowable region
Box(4) = {-Fault_length/2, -Fault_width, -HighDamage_thickness/2, Fault_length, Fault_width, HighDamage_thickness};

// Boolean operation to fragment all volumes
BooleanFragments{ Volume{1,2,3,4}; Delete; }{}

// Define mesh sizes using a progression field for smooth transition

// Field 1: Mesh size inside the fault zone
Field[1] = Box;
Field[1].VIn = lc_fault;
Field[1].VOut = lc;
Field[1].XMin = -Fault_length/2;
Field[1].XMax = Fault_length/2;
Field[1].YMin = Ymin/2-Fault_width/2;
Field[1].YMax = 100 * Ymax;
Field[1].ZMin = -Fault_thickness/2;
Field[1].ZMax = Fault_thickness/2;

Background Field = 1;

// Mark all volumes as physical volumes
volumes[] = Volume{:};
For i In {0:#volumes[]-1}
    Physical Volume(Sprintf("Volume_%g", i+1)) = {volumes[i]};
EndFor

// Print the number of volumes created
Printf("Number of volumes created: %g", #volumes[]);

