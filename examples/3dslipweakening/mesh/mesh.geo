
// GMSH script for creating a structured hex mesh with refinement near the xy plane
// Box dimensions: x=(-60000,60000), y=(-60000,60000), z=(-60000,60000)
// Refined area: x=(-10000,10000), y=(0,-10000) in the xy plane
// Refinement: 100m elements near xy plane, growing to 10000m at box boundaries

// Set mesh verbosity level
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Define characteristic lengths for mesh refinement
lc_fine = 100;     // Fine mesh size near xy plane (100m)
lc_coarse = 10000; // Coarse mesh size at box boundaries (10000m)

// Box boundaries
x_min = -60000; x_max = 60000;
y_min = -60000; y_max = 60000;
z_min = -60000; z_max = 60000;

// Refinement zone boundaries
rx_min = -10000; rx_max = 10000;
ry_min = -10000; ry_max = 0;
rz_min = 0; rz_max = 0;  // xy plane is at z=0

// Create all points for the outer box
// Bottom face (z = z_min)
p1 = newp; Point(p1) = {x_min, y_min, z_min, lc_coarse};
p2 = newp; Point(p2) = {x_max, y_min, z_min, lc_coarse};
p3 = newp; Point(p3) = {x_max, y_max, z_min, lc_coarse};
p4 = newp; Point(p4) = {x_min, y_max, z_min, lc_coarse};

// Top face (z = z_max)
p5 = newp; Point(p5) = {x_min, y_min, z_max, lc_coarse};
p6 = newp; Point(p6) = {x_max, y_min, z_max, lc_coarse};
p7 = newp; Point(p7) = {x_max, y_max, z_max, lc_coarse};
p8 = newp; Point(p8) = {x_min, y_max, z_max, lc_coarse};

// Create intermediate points for refinement in z-direction
// Bottom refined region (z = 0, xy-plane)
p9 = newp; Point(p9) = {x_min, y_min, 0, lc_coarse};
p10 = newp; Point(p10) = {x_max, y_min, 0, lc_coarse};
p11 = newp; Point(p11) = {x_max, y_max, 0, lc_coarse};
p12 = newp; Point(p12) = {x_min, y_max, 0, lc_coarse};

// Points at the refined region in xy-plane
p13 = newp; Point(p13) = {rx_min, ry_min, 0, lc_fine};
p14 = newp; Point(p14) = {rx_max, ry_min, 0, lc_fine};
p15 = newp; Point(p15) = {rx_max, ry_max, 0, lc_fine};
p16 = newp; Point(p16) = {rx_min, ry_max, 0, lc_fine};

// Create lines for bottom face (z = z_min)
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};

// Create lines for top face (z = z_max)
l5 = newl; Line(l5) = {p5, p6};
l6 = newl; Line(l6) = {p6, p7};
l7 = newl; Line(l7) = {p7, p8};
l8 = newl; Line(l8) = {p8, p5};

// Create vertical lines connecting bottom and top faces
l9 = newl; Line(l9) = {p1, p9};
l10 = newl; Line(l10) = {p2, p10};
l11 = newl; Line(l11) = {p3, p11};
l12 = newl; Line(l12) = {p4, p12};

// Create lines for middle xy-plane (z = 0)
l13 = newl; Line(l13) = {p9, p10};
l14 = newl; Line(l14) = {p10, p11};
l15 = newl; Line(l15) = {p11, p12};
l16 = newl; Line(l16) = {p12, p9};

// Create vertical lines connecting middle to top faces
l17 = newl; Line(l17) = {p9, p5};
l18 = newl; Line(l18) = {p10, p6};
l19 = newl; Line(l19) = {p11, p7};
l20 = newl; Line(l20) = {p12, p8};

// Create lines for the refined region in xy-plane
l21 = newl; Line(l21) = {p13, p14};
l22 = newl; Line(l22) = {p14, p15};
l23 = newl; Line(l23) = {p15, p16};
l24 = newl; Line(l24) = {p16, p13};

// Connect refined region to outer boundary
l25 = newl; Line(l25) = {p9, p13};
l26 = newl; Line(l26) = {p10, p14};
l27 = newl; Line(l27) = {p11, p15};
l28 = newl; Line(l28) = {p12, p16};

// Create line loops for surfaces
// Bottom volume (z < 0)
ll1 = newll; Line Loop(ll1) = {l1, l2, l3, l4}; // Bottom face
ll2 = newll; Line Loop(ll2) = {l13, l14, l15, l16}; // Middle face at z=0
ll3 = newll; Line Loop(ll3) = {l1, l10, -l13, -l9}; // Front face
ll4 = newll; Line Loop(ll4) = {l2, l11, -l14, -l10}; // Right face
ll5 = newll; Line Loop(ll5) = {l3, l12, -l15, -l11}; // Back face
ll6 = newll; Line Loop(ll6) = {l4, l9, -l16, -l12}; // Left face

// Create surfaces for bottom volume
s1 = news; Plane Surface(s1) = {ll1};
s2 = news; Plane Surface(s2) = {ll2};
s3 = news; Plane Surface(s3) = {ll3};
s4 = news; Plane Surface(s4) = {ll4};
s5 = news; Plane Surface(s5) = {ll5};
s6 = news; Plane Surface(s6) = {ll6};

// Top volume (z > 0)
ll7 = newll; Line Loop(ll7) = {l5, l6, l7, l8}; // Top face
ll8 = newll; Line Loop(ll8) = {l5, -l18, -l13, l17}; // Front face
ll9 = newll; Line Loop(ll9) = {l6, -l19, -l14, l18}; // Right face
ll10 = newll; Line Loop(ll10) = {l7, -l20, -l15, l19}; // Back face
ll11 = newll; Line Loop(ll11) = {l8, -l17, -l16, l20}; // Left face

// Create surfaces for top volume
s7 = news; Plane Surface(s7) = {ll7};
s8 = news; Plane Surface(s8) = {ll8};
s9 = news; Plane Surface(s9) = {ll9};
s10 = news; Plane Surface(s10) = {ll10};
s11 = news; Plane Surface(s11) = {ll11};

// Create refined area in xy-plane
ll12 = newll; Line Loop(ll12) = {l21, l22, l23, l24}; // Refined area
s12 = news; Plane Surface(s12) = {ll12};

// Create sub-regions in xy-plane
ll13 = newll; Line Loop(ll13) = {l25, l21, -l26, -l13};
ll14 = newll; Line Loop(ll14) = {l26, l22, -l27, -l14};
ll15 = newll; Line Loop(ll15) = {l27, l23, -l28, -l15};
ll16 = newll; Line Loop(ll16) = {l28, l24, -l25, -l16};

s13 = news; Plane Surface(s13) = {ll13};
s14 = news; Plane Surface(s14) = {ll14};
s15 = news; Plane Surface(s15) = {ll15};
s16 = news; Plane Surface(s16) = {ll16};

// Create volume for bottom part (z < 0)
sl1 = newsl; Surface Loop(sl1) = {s1, s2, s3, s4, s5, s6};
v1 = newv; Volume(v1) = {sl1};

// Create volume for top part (z > 0)
sl2 = newsl; Surface Loop(sl2) = {s2, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16};
v2 = newv; Volume(v2) = {sl2};

// Set up transfinite mesh
// Number of nodes in each direction - adjust these values to control mesh density
nx = 15; // Number of nodes in x-direction
ny = 15; // Number of nodes in y-direction
nz_bottom = 20; // Number of nodes in z-direction below xy-plane
nz_top = 20; // Number of nodes in z-direction above xy-plane
nx_fine = 10; // Number of nodes in x-direction in refined area
ny_fine = 10; // Number of nodes in y-direction in refined area

// Set transfinite curves - outer box edges with progression
// Bottom face edges
Transfinite Curve {l1, l3} = nx Using Progression 1.1;
Transfinite Curve {l2, l4} = ny Using Progression 1.1;

// Top face edges
Transfinite Curve {l5, l7} = nx Using Progression 1.1;
Transfinite Curve {l6, l8} = ny Using Progression 1.1;

// Middle plane edges
Transfinite Curve {l13, l15} = nx Using Progression 1.1;
Transfinite Curve {l14, l16} = ny Using Progression 1.1;

// Vertical edges with progression towards the xy-plane
Transfinite Curve {l9, l10, l11, l12} = nz_bottom Using Progression 0.9; // Higher density near z=0
Transfinite Curve {l17, l18, l19, l20} = nz_top Using Progression 1.1; // Higher density near z=0

// Refined region edges
Transfinite Curve {l21, l23} = nx_fine;
Transfinite Curve {l22, l24} = ny_fine;

// Connecting edges from boundary to refined region
Transfinite Curve {l25, l27} = (nx - nx_fine) / 2 + 1 Using Progression 0.8;
Transfinite Curve {l26, l28} = (ny - ny_fine) / 2 + 1 Using Progression 0.8;

// Set transfinite surfaces
Transfinite Surface {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16};

// Set transfinite volumes
Transfinite Volume {v1, v2};

// Set recombine to create structured hex elements
Recombine Surface {s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16};

// Generate the mesh
Mesh.Algorithm = 6; // Frontal-Delaunay for quads
Mesh.Algorithm3D = 1; // Delaunay for tetrahedra
Mesh.SubdivisionAlgorithm = 0; // No subdivision
Mesh.Optimize = 1; // Optimize the mesh
Mesh.OptimizeNetgen = 1; // Optimize with Netgen algorithm
