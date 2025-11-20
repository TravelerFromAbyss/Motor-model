size = 2.5;

maxXY = 54;
minXY = 10;

// Outside points
Point(1) = {maxXY, maxXY, 0, size};
Point(2) = {-maxXY, maxXY, 0, size};
Point(3) = {-maxXY, -maxXY, 0, size};
Point(4) = {maxXY, -maxXY, 0, size};

//Outside lines
Line(1) = {1, 2};
Line(2) = {4, 1};
Line(3) = {3, 4};
Line(4) = {2, 3};

// Inside points
Point(5) = {minXY, minXY, 0, size / 3};
Point(6) = {-minXY, minXY, 0, size / 3};
Point(7) = {-minXY, -minXY, 0, size / 3};
Point(8) = {minXY, -minXY, 0, size / 3};

//Inside lines
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Physical Curve("Boundary", 1) = {1, 2, 3, 4};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
//Physical Curve("Outer", 2) = {1, 2, 3, 4, 5, 6, 7, 8};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};
Physical Curve("Inner", 3) = {5, 6, 7, 8};

Physical Surface("Steel", 1) = {1};
Physical Surface("Wire", 2) = {2};

