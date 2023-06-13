// Gmsh project created on Fri Aug 16 13:09:28 2019
SetFactory("OpenCASCADE");

remove = 0;

// Points to be used
Point(1) = {0.05-remove, 0, 0, 1.0};
Point(2) = {0.05-remove, 31.4, 0, 1.0};
Point(3) = {1.75-remove, 31.4, 0, 1.0};
Point(4) = {1.75-remove, 0, 0, 1.0};
Point(5) = {2.05-remove, 0, 0, 1.0};
Point(6) = {2.05-remove, 31.4, 0, 1.0};
Point(7) = {0.35, 31.4, 0, 1.0};
Point(8) = {0.35, 0, 0, 1.0};


//+
Line(1) = {2, 7};
//+
Line(2) = {7, 8};
//+
Line(3) = {8, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {7, 3};
//+
Line(6) = {3, 4};
//+
Line(7) = {4, 8};
//+
Line(8) = {4, 5};
//+
Line(9) = {5, 6};
//+
Line(10) = {6, 3};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Line Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(1) = {2};
//+
Line Loop(3) = {5, 6, 7, -2};
//+
Plane Surface(2) = {3};
//+
Line Loop(4) = {10, 6, 8, 9};
//+
Plane Surface(3) = {4};
//+
Transfinite Line {1, 3, 10, 8} = 100 Using Progression 1;
//+
Transfinite Line {5, 7} = 50 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Recombine Surface {1};
//+
Extrude {0, 0, -10} {
  Surface{1}; Surface{2}; Surface{3}; 
  Layers{1};
  Recombine;
}


//+
Physical Volume("internalField") = {1, 2, 3};
//+
Physical Surface("otherfaces") = {8, 4, 9, 12, 16, 13, 2, 1, 3, 14, 11, 6};
//+
Physical Surface("dielectric") = {7, 15};
