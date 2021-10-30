a=0.01;
b=0.04;



lbeam = a/10;
lfine = a/5;
lm=a/5;



Point(1) = {0, 0, 0, lbeam};

Point(2) = {a, 0, 0, lbeam};
Point(3) = {0,a, 0, lbeam};








Line(1) = {1,2};
Line(2) = {1,3};





Circle(9) = {2,1,3};


Symmetry { 1,0,0,0 } { Duplicata{Line{9}; Line{1};} }

Symmetry { 0,1,0,0 } { Duplicata{Line{9}; Line{10}; Line{2};} }


Line Loop(30)={9,-10,13,-12};
Plane Surface(1) = {30};
//Physical Line(1) ={1};
Physical Surface(1) = {1};   //beam




Point(11) = {b,0, 0, lfine};
Point(12) = {0,b, 0, lfine};
Point(13) = {-b,0, 0, lfine};
Point(14) = {0,-b, 0, lfine};


Circle(21) = {11,1,12};
Circle(22) = {12,1,13};
Circle(23) = {13,1,14};
Circle(24) = {14,1,11};

Line Loop(31)={21,22,23,24};
Plane Surface(2) = {31,30};
Physical Surface(2) = {2};   //vacuum





//For Mesh
tx=a*0.4;
ty=a*0.3;
Point(115)={tx,0,0,lbeam};
Point(116)={tx,ty,0,lbeam};
Point(117)={0,ty,0,lbeam};
Point(118)={-tx,ty,0,lbeam};
Point(119)={-tx,0,0,lbeam};
Point(120)={-tx,-ty,0,lbeam};
Point(121)={0,-ty,0,lbeam};
Point(122)={tx,-ty,0,lbeam};


aa=3.0*a;
aaa=aa/1.4;
Point(215)={aa,0,0,lm};
Point(216)={aaa,aaa,0,lm};
Point(217)={0,aa,0,lm};
Point(218)={-aaa,aaa,0,lm};
Point(219)={-aa,0,0,lm};
Point(220)={-aaa,-aaa,0,lm};
Point(221)={0,-aa,0,lm};
Point(222)={aaa,-aaa,0,lm};


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
