b_0 = 0.04;
h_0 = 0.02;
Rq = 0.05;
L = 20*h_0;
R = L*Rq;
C_b = 0.5;
C_h = 0.5;
N_e = 40;
b_1 = b_0*(1-C_b);
h_1 = h_0*(1-C_h);
Point(1) = {R, b_0/2, h_0/2};
Point(2) = {R, b_0/2,-h_0/2};
Point(3) = {R,-b_0/2,-h_0/2};
Point(4) = {R,-b_0/2, h_0/2};
//
Point(5) = {R+L, b_1/2, h_1/2};
Point(6) = {R+L, b_1/2,-h_1/2};
Point(7) = {R+L,-b_1/2,-h_1/2};
Point(8) = {R+L,-b_1/2, h_1/2};//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 4};
//+
Line(10) = {7, 3};
//+
Line(11) = {2, 6};
//+
Line(12) = {5, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -4, -10, 8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, 2, 11, -6};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, 8, 5, 6};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {11, 7, 10, -3};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {9, 1, -12, -5};
//+
Plane Surface(6) = {6};
//+
Transfinite Curve {9, 10, 12, 11} = 41 Using Progression 1;
//+
Transfinite Curve {1, 3, 7, 5} = 6 Using Progression 1;
//+
Transfinite Curve {4, 2, 8, 6} = 6 Using Progression 1;
//+
Transfinite Surface {6} = {4, 1, 5, 8};
//+
Transfinite Surface {2} = {4, 8, 7, 3};
//+
Transfinite Surface {5} = {2, 6, 7, 3};
//+
Transfinite Surface {3} = {6, 5, 1, 2};
//+
Transfinite Surface {4} = {7, 6, 5, 8};
//+
Transfinite Surface {1} = {4, 1, 2, 3};
//+
Recombine Surface {2, 6, 3, 5, 4, 1};
//+
Surface Loop(1) = {6, 2, 1, 3, 5, 4};
//+
Volume(1) = {1};
//+
Transfinite Volume{1} = {7, 6, 2, 3, 8, 5, 1, 4};
//+
Physical Volume("BEAM", 13) = {1};
