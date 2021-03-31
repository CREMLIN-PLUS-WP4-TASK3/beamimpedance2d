//beam radius
a=0.0005;

// inner right beam screen radius
rx=0.018;
//outer right beam screen radius
rx2=0.019;

ry=0.01273;
rx3=0.00973;

ry2=0.014;
rx4=0.0102;

xhole=0.01493;
yhole=0.00157;

xhole2=0.01584;
yhole2=0.0025;

h=0.0215;

lvac = 0.002;
lpipe =6e-4;
lfine = 1e-4;
lveryfine = 9e-5;

eps=0.000005;


// origin
Point(1) = {0, 0, 0, lfine};

// beam

Point(2) = {a-eps, 0, 0, lveryfine};
Point(3) = {a,0, 0, lveryfine};
Point(4) = {a+eps,0, 0, 1.0*lveryfine};


Point(12) = {0,a-eps, 0, lveryfine};
Point(13) = {0,a,0,  lveryfine};
Point(14) = {0,a+eps, 0, 1.0*lveryfine};


Point(22) = {-(a-eps), 0, 0, lveryfine};
Point(23) = {-a,0, 0, lveryfine};
Point(24) = {-(a+eps),0, 0, 1.0*lveryfine};


Point(32) = {0,-(a-eps), 0, lveryfine};
Point(33) = {0,-a,0,  lveryfine};
Point(34) = {0,-(a+eps), 0, 1.0*lveryfine};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

Line(11) = {1,12};
Line(12) = {12,13};
Line(13) = {13,14};

Line(21) = {1,22};
Line(22) = {22,23};
Line(23) = {23,24};

Line(31) = {1,32};
Line(32) = {32,33};
Line(33) = {33,34};

Circle(51) = {2,1,12};
Circle(52) = {3,1,13};
Circle(53) = {4,1,14};

Circle(61) = {12,1,22};
Circle(62) = {13,1,23};
Circle(63) = {14,1,24};

Circle(71) = {22,1,32};
Circle(72) = {23,1,33};
Circle(73) = {24,1,34};

Circle(81) = {32,1,2};
Circle(82) = {33,1,3};
Circle(83) = {34,1,4};


Line Loop(101)={1,51,-11};
Line Loop(102)={2,52,-12,-51};
Line Loop(103)={3,53,-13,-52};

Line Loop(111)={11,61,-21};
Line Loop(112)={12,62,-22,-61};
Line Loop(113)={13,63,-23,-62};

Line Loop(121)={21,71,-31};
Line Loop(122)={22,72,-32,-71};
Line Loop(123)={23,73,-33,-72};

Line Loop(131)={31,81,-1};
Line Loop(132)={32,82,-2,-81};
Line Loop(133)={33,83,-3,-82};




// pipe
Point(41) = {rx, 0, 0, lpipe};
//Point(42) = {-rx+ry-rx3, 0, 0, lpipe};

Point(51) = {rx2, 0, 0, lpipe};
//Point(52) = {-rx2+ry-rx3, 0, 0, lpipe};

Point(43) = {ry, ry, 0, lpipe};
Point(44) = {ry, -ry, 0, lpipe};
Point(45) = {-rx3, ry, 0, lpipe};
Point(46) = {-rx3, -ry, 0, lpipe};

Point(53) = {ry2, ry2, 0, lpipe};
Point(54) = {ry2, -ry2, 0, lpipe};
Point(55) = {-rx4, ry2, 0, lpipe};
Point(56) = {-rx4, -ry2, 0, lpipe};

Point(57) = {(2*ry2*ry2-rx2*rx2)/(2*(ry2-rx2)), 0, 0, lpipe};
Point(58) = {(rx4*rx4+ry2*ry2-(ry-rx2-rx3)*(ry-rx2-rx3))/(2*(rx2-ry-rx4+rx3)), 0, 0, lpipe};

Point(47) = {(rx3*rx3+ry*ry-(ry-rx-rx3)*(ry-rx-rx3))/(2*(rx-ry)),0,0,lpipe};

//Point(61) = {-xhole,yhole,0,lpipe};
//Point(62) = {-xhole,-yhole,0,lpipe};
//Point(63) = {-xhole2,yhole2,0,lpipe};
//Point(64) = {-xhole2,-yhole2,0,lpipe};

Line(101) = {43,45};
Line(102) = {44,46};

Circle(103) = {43,1,44};
Circle(104) = {45,47,46};
//Circle(105) = {62,47,46};
//Circle(106) = {61,47,62};


Line(111) = {53,55};
Line(112) = {54,56};

Circle(113) = {53,57,54};
Circle(114) = {55,58,56};
//Circle(115) = {64,58,56};
//Circle(116) = {63,58,64};

//Line(116) = {61,63}; //hole
//Line(117) = {62,64};

//Line Loop(301) = {-112,-113,111,114,-116,-104,-101,103,102,-105,117,115};
Line Loop(301) = {-102,-103,101,104};
Line Loop(303) = {-112,-113,111,114};

// outside vacuum
Point(71) = {h,0, 0, lvac};
Point(72) = {0,h, 0, lvac};
Point(73) = {-h,0, 0, lvac};
Point(74) = {0,-h, 0, lvac};

Circle(201) = {71,1,72};
Circle(202) = {72,1,73};
Circle(203) = {73,1,74};
Circle(204) = {74,1,71};

Line Loop(302)={201,202,203,204};

// surfaces
// beam and inner vac
Plane Surface(1) = {101};
Plane Surface(2) = {102};
Plane Surface(3) = {103};

Plane Surface(4) = {111};
Plane Surface(5) = {112};
Plane Surface(6) = {113};

Plane Surface(7) = {121};
Plane Surface(8) = {122};
Plane Surface(9) = {123};

Plane Surface(10) = {131};
Plane Surface(11) = {132};
Plane Surface(12) = {133};

// pipe
Plane Surface(13) = {303,301};

//vac
Plane Surface(14) = {301,103,113,123,133};
Plane Surface(15) = {302,303};


//Physical surfaces

Physical Surface(1) = {2,3,5,6,8,9,11,12};   //beam
Physical Surface(2) = {1,4,7,10,14,15};   //vacuum
Physical Surface(4) = {13};   //pipe


Color Green{ Surface{ 14:15 }; }




