/**
 * Coarsened variant of v3_100m_freesurface.geo
 * - Keeps 100 m elements near the fault and nucleation patch
 * - Coarsens aggressively away from the fault to reduce total elements
 * - Confines 100 m to a thin 3D slab around the fault plane (y ~ 0)
 * - Disables implicit point/curvature refinements to let fields control sizing
 */

// Base far-field size (coarser than original 1e4)
lc = 2e4;
// Minimum size near fault and nucleation
lc_fault = 100;

Fault_length = 30e3;
Fault_width = 15e3;
Fault_dip = 90*Pi/180.;

// Nucleation in X,Z local coordinates
X_nucl = 0e3;
Width_nucl = 0.5*Fault_width;
R_nucl = 1.5e3;
lc_nucl = 100;

Xmax = 40e3;
Xmin = -Xmax;
Ymin = -Xmax +  0.5 * Fault_width  *Cos(Fault_dip);
Ymax =  Xmax + 0.5 * Fault_width  *Cos(Fault_dip);
Zmin = -200e3;

// Move the fault to the center in depth
move_distance = -2e3;

// Create the Volume
Point(1) = {Xmin, Ymin, 0, lc};
Point(2) = {Xmin, Ymax, 0, lc};
Point(3) = {Xmax, Ymax, 0, lc};
Point(4) = {Xmax, Ymin, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};
Extrude {0,0, Zmin} { Surface{1}; }

// Create the fault - CENTERED IN DOMAIN
Point(100) = {-0.5*Fault_length, Fault_width *Cos(Fault_dip), -Fault_width  *Sin(Fault_dip) + move_distance, lc};
Point(101) = {-0.5*Fault_length, 0, move_distance, lc};
Point(102) = {0.5*Fault_length, 0,  move_distance, lc};
Point(103) = {0.5*Fault_length, Fault_width  *Cos(Fault_dip), -Fault_width  *Sin(Fault_dip) + move_distance, lc};
Line(100) = {100, 101};
Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 100};

// Create nucleation patch (rectangular ring)
Point(201) = {X_nucl + R_nucl , (Width_nucl + R_nucl) * Cos(Fault_dip), -(Width_nucl+R_nucl)  *Sin(Fault_dip), lc_nucl};
Point(202) = {X_nucl + R_nucl , (Width_nucl - R_nucl) * Cos(Fault_dip), -(Width_nucl-R_nucl)  *Sin(Fault_dip), lc_nucl};
Point(203) = {X_nucl - R_nucl , (Width_nucl - R_nucl) * Cos(Fault_dip), -(Width_nucl-R_nucl)  *Sin(Fault_dip), lc_nucl};
Point(204) = {X_nucl - R_nucl , (Width_nucl + R_nucl) * Cos(Fault_dip), -(Width_nucl+R_nucl)  *Sin(Fault_dip), lc_nucl};
Line(200) = {201, 202};
Line(201) = {202, 203};
Line(202) = {203, 204};
Line(203) = {204, 201};
Curve Loop(204) = {200,201,202,203};
Plane Surface(200) = {204};

Curve Loop(105) = {100,101,102,103};
Plane Surface(100) = {105, 204};

// There is a bug in "Attractor", we need to define a Ruled surface in FaceList
Line Loop(106) = {100,101,102,103};
Ruled Surface(101) = {106};
Ruled Surface(201) = {204};

Surface{100,200} In Volume{1};

// SIZING FIELDS
// Distance to the fault surface
Field[1] = Distance;
Field[1].FacesList = {101};

// Faster growth with distance: keep 100 m at fault, ramp quickly to km-scale
Field[2] = MathEval;
Field[2].F = Sprintf("%g + 0.2*F1 + (F1/1.8e3)^2", lc_fault);

// Nucleation refinement (tight halo, still 100 m at core)
Field[3] = Distance;
Field[3].FacesList = {201};

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc_nucl;
Field[4].LcMax = lc_fault;
Field[4].DistMin = 0.75*R_nucl;
Field[4].DistMax = 1.5*R_nucl;

Field[5] = Restrict;
Field[5].IField = 4;
Field[5].FacesList = {100,200};

// Near-fault transition band (clamp to far-field size beyond ~2 km from fault)
Field[6] = Threshold;
Field[6].IField = 1;
Field[6].LcMin = lc_fault;
Field[6].LcMax = lc;
Field[6].DistMin = 1.0e3;   // within 1 km -> keep lc_fault
Field[6].DistMax = 2.0e3;   // beyond 2 km -> use lc (far-field)

// Combine all controls
Field[7] = Min;
Field[7].FieldsList = {2,5,6};

Background Field = 7;

// Physical groups
Physical Surface(101) = {1};
Physical Surface(103) = {100,200};
Physical Surface(105) = {14,18,22,26,27};

Physical Volume(1) = {1};

// Let fields fully control sizes and keep file format compatibility
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.MshFileVersion = 2.2;
