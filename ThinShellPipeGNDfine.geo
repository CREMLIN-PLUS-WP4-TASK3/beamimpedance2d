a=0.01;
b=0.04;
h=0.0403;
h2=0.1;

GNDAngle=Pi/10;
Printf("Angles'%g' and '%g'", Cos(GNDAngle), Sin(GNDAngle));

lbeam = a;
lvac = 0.05;
lpipe =0.0001;
lGND=0.05;



Point(1) = {0, 0, 0, lbeam};

Point(2) = {a, 0, 0, lbeam};
Point(3) = {0,a, 0, lbeam};








Line(1) = {1,2};
Line(2) = {1,3};





Circle(9) = {2,1,3};


Symmetry { 1,0,0,0 } { Duplicata{Line{9}; Line{1};} }

Symmetry { 0,1,0,0 } { Duplicata{Line{9}; Line{10}; Line{2};} }


Line Loop(30)={9,-10,13,-12};





Point(11) = {b,0, 0, lpipe};
Point(12) = {0,b, 0, lpipe};
Point(13) = {-b,0, 0, lpipe};
Point(14) = {0,-b, 0, lpipe};


Circle(21) = {11,1,12};
Circle(22) = {12,1,13};
Circle(23) = {13,1,14};
Circle(24) = {14,1,11};

Line Loop(31)={21,22,23,24};



//The pipe
Point(21) = {h,0, 0, lpipe};
Point(22) = {0,h, 0, lpipe};
Point(23) = {-h,0, 0, lpipe};
Point(24) = {0,-h, 0, lpipe};

Point(25) = {h*Cos(GNDAngle),h*Sin(GNDAngle), 0, lpipe};

Circle(30) = {21,1,25};

Circle(31) = {25,1,22};
Circle(32) = {22,1,23};
Circle(33) = {23,1,24};
Circle(34) = {24,1,21};



Line Loop(32)={30,31,32,33,34};



//The outside vacuum
Point(31) = {h2,0, 0, lGND};
Point(32) = {0,h2, 0, lvac};
Point(33) = {-h2,0, 0, lvac};
Point(34) = {0,-h2, 0, lvac};

Point(35) = {h2*Cos(GNDAngle),h2*Sin(GNDAngle), 0, lGND};


Circle(40) = {31,1,35};

Circle(41) = {35,1,32};
Circle(42) = {32,1,33};
Circle(43) = {33,1,34};
Circle(44) = {34,1,31};

Line Loop(33)={40,41,42,43,44};


//For GND
Line(100)={21,31};
Line(101)={25,35};
Line Loop(34)={100,40,-101,-30};



//Physical Surfaces

Plane Surface(1) = {30};


Plane Surface(2) = {31,30};


Plane Surface(3) = {32,31,30}; //Pipe
Plane Surface(5) = {34}; //GND





Plane Surface(4) = {33,32,31,30,34};
//Physical Surface(4) = {4};   //OutsideVac


Physical Surface(1) = {1};   //beam
Physical Surface(2) = {2,4};   //vacuum + Outsidevac
Physical Surface(3) = {3};   //Pipe
Physical Surface(6) = {5};   //GND



//For Mesh
t=a*0.5;
Point(115)={t,0,0,lbeam};
Point(116)={t,t,0,lbeam};
Point(117)={0,t,0,lbeam};
Point(118)={-t,t,0,lbeam};
Point(119)={-t,0,0,lbeam};
Point(120)={-t,-t,0,lbeam};
Point(121)={0,-t,0,lbeam};
Point(122)={t,-t,0,lbeam};


aa=1.5*a;
Point(215)={aa,0,0,lbeam};
Point(216)={a,a,0,lbeam};
Point(217)={0,aa,0,lbeam};
Point(218)={-a,a,0,lbeam};
Point(219)={-aa,0,0,lbeam};
Point(220)={-a,-a,0,lbeam};
Point(221)={0,-aa,0,lbeam};
Point(222)={a,-a,0,lbeam};


Point {1}  In Surface{1};   //Force mesh on these points

Point {115}  In Surface{1};
Point {116}  In Surface{1};
Point {117}  In Surface{1};
Point {118}  In Surface{1};
Point {119}  In Surface{1};
Point {120}  In Surface{1};
Point {121}  In Surface{1};
Point {122}  In Surface{1};

Point {215}  In Surface{2};
Point {216}  In Surface{2};
Point {217}  In Surface{2};
Point {218}  In Surface{2};
Point {219}  In Surface{2};
Point {220}  In Surface{2};
Point {221}  In Surface{2};
Point {222}  In Surface{2};

