//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.40, 0, 0, 1.0};
//+
Point(3) = {1.60, 0, 0, 1.0};
//+
Point(4) = {2, 0, 0, 1.0};
//+
Point(5) = {0, 20, 0, 1.0};
//+
Point(6) = {0.40, 20, 0, 1.0};
//+
Point(7) = {1.60, 20, 0, 1.0};
//+
Point(8) = {2, 20, 0, 1.0};
//+
Line(1) = {5, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 3};
//+
Line(7) = {2, 3};
//+
Line(8) = {3, 4};
//+
Line(9) = {4, 8};
//+
Line(10) = {8, 7};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {7, -6, -5, -3};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {10, 6, 8, 9};
//+
Plane Surface(3) = {3};
//+
Transfinite Line {4, 2, 8, 10} = 100 Using Progression 1;
//+
Transfinite Line {5, 7} = 25 Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Recombine Surface {1, 2, 3};
//+
Extrude {0, 0, 400} {
  Surface{1}; Surface{2}; Surface{3};
  Layers{1};
  Recombine; 
}
//+
Physical Surface("otherFaces") = {54, 32, 76, 1, 2, 3, 49, 31, 63, 23, 41, 71};
//+
Physical Surface("electrode") = {19};
//+
Physical Surface("ground") = {75};
//+
Physical Volume("internalfield") = {2, 1, 3};
