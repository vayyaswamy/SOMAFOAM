// Gmsh project created on Thu Feb 14 15:39:37 2019
SetFactory("OpenCASCADE");

// Parameters to change if desired - This is for uniform mesh so there is that
lhs1 = 0;
lhs2 = 1;
gapsize = 1;

// Mesh 
mesh_size = 41;

// Points of dielectric (left side)
Point(1) = {-0.95+0, 0, 0, 1.0};
Point(2) = {-0.95+0, 31.4, 0, 1.0};
Point(3) = {0.05+0, 31.4, 0, 1.0};
Point(4) = {0.05+0, 0, 0, 1.0};

// Points of dielectric (left side)
Point(5) = {lhs1+2.05, 0, 0, 1.0};
Point(6) = {lhs1+2.05, 31.4, 0, 1.0};
Point(7) = {lhs2+2.05, 31.4, 0, 1.0};
Point(8) = {lhs2+2.05, 0, 0, 1.0};


//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 1};
//+
Line(4) = {1, 2};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 5};
//+
Line(8) = {5, 6};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};
//+
Transfinite Line {1, 3, 5, 7} = 10 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Recombine Surface {1, 2};
//+
Extrude {0, 0, -10} {
  Surface{1}; Surface{2};
  	Layers{1};	
	Recombine;
}

//+
Physical Volume("internal") = {1, 2};
//+
Physical Surface("electrode") = {6};
//+
Physical Surface("ground") = {9};
//+
Physical Surface("plasma") = {11, 4};
//+
Physical Surface("otherFaces") = {1, 2, 12, 7, 5, 10, 3, 8};
