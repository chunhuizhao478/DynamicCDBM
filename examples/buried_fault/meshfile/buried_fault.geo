SetFactory("OpenCASCADE");

lc = 5e3;
lc_fault = 50; // Fine mesh in the fault zone

Fault_length = 5e3;
Fault_width = 5e3;
Fault_thickness = 1000;
Fault_dip = Pi/2;

// Nucleation in X,Z local coordinates
R_nucl = 400;
thickness_nucl = 200;

Xmax = 10e3;
Xmin = -Xmax;

Ymax = 0;
Ymin = -Xmax;

//
Width_nucl = Ymin/2;

Zmin =  -Xmax +  0.5 * Fault_width  *Cos(Fault_dip);
Zmax =   Xmax + 0.5 * Fault_width  *Cos(Fault_dip);

Box(1) = {Xmin, 0, Zmin, 2*Xmax, Ymin, 2*Zmax};

// Create a damage zone
Box(2) = {-Fault_length/2, Ymin/2-Fault_width/2, -Fault_thickness/2, Fault_length, Fault_width, Fault_thickness};

// Create a nucleation patch
//Box(3) = {X_nucl-R_nucl/2, Width_nucl-R_nucl/2, -thickness_nucl/2, R_nucl, R_nucl, thickness_nucl};
Box(3) = {-R_nucl/2, Width_nucl-R_nucl/2, -thickness_nucl/2, R_nucl, R_nucl, thickness_nucl};

// Create a cdbm allowable region
Box(4) = {-Fault_length/2, Ymin/2-Fault_width/2, -thickness_nucl/2, Fault_length, Fault_width, thickness_nucl};

// Boolean operation to fragment all volumes
BooleanFragments{ Volume{1,2,3,4}; Delete; }{}

// Define mesh sizes using a progression field for smooth transition

// Field 1: Mesh size inside the fault zone
Field[1] = Box;
Field[1].VIn = lc_fault;
Field[1].VOut = lc/4;
Field[1].XMin = -Fault_length/2;
Field[1].XMax = Fault_length/2;
Field[1].YMin = Ymin/2-Fault_width/2;
Field[1].YMax = Ymin/2+Fault_width/2;
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

