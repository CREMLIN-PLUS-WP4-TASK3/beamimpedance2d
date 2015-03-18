a=0.001;

rx=0.018;
rx2=0.019;

ry=0.01273;
rx3=0.00973;

ry2=0.014;
rx4=0.0102;

xhole=0.01493;
yhole=0.00157;

xhole2=0.01584;
yhole2=0.0025;

//h=0.0215;
h=0.03;

lvac = 1e-3;
lpipe = 1e-3;
lcooling = 1e-3;
lfine = 5e-4;
lveryfine = 1e-4;
lcoating = 1e-3;
lsing=4e-4;

eps=0.000005;
d1=0.00008;


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




// pipe & coating
Point(41) = {rx, 0, 0, lcoating};
Point(1041)={rx-d1, 0, 0, lcoating};

Point(51) = {rx2, 0, 0, lpipe};

Point(43) = {ry, ry, 0, lcoating};
Point(1043) = {Sqrt((Sqrt(2)*ry-d1)^2-(ry-d1)^2), ry-d1, 0 ,lcoating};
Point(44) = {ry, -ry, 0, lcoating};
Point(1044) = {Sqrt((Sqrt(2)*ry-d1)^2-(ry-d1)^2), -ry+d1, 0, lcoating};
Point(45) = {-rx3, ry, 0, lcoating};
Point(1045) = {-Sqrt((Sqrt(2)*(-rx3)+d1)^2-(-rx3+d1)^2), ry-d1, 0 ,lcoating};
Point(46) = {-rx3, -ry, 0, lcoating};
Point(1046) = {-Sqrt((Sqrt(2)*(-rx3)+d1)^2-(-rx3+d1)^2), -ry+d1, 0, lcoating};

Point(53) = {ry2, ry2, 0, lpipe};
Point(54) = {ry2, -ry2, 0, lpipe};
Point(55) = {-rx4, ry2, 0, lpipe};
Point(56) = {-rx4, -ry2, 0, lpipe};

Point(57) = {(2*ry2*ry2-rx2*rx2)/(2*(ry2-rx2)), 0, 0, lpipe};
Point(58) = {(rx4*rx4+ry2*ry2-(ry-rx2-rx3)*(ry-rx2-rx3))/(2*(rx2-ry-rx4+rx3)), 0, 0, lpipe};

Point(47) = {(rx3*rx3+ry*ry-(ry-rx-rx3)*(ry-rx-rx3))/(2*(rx-ry)),0,0,lpipe};

Point(61) = {-xhole,yhole,0,lsing};
Point(1061) = {-xhole+d1*Cos(ArcTan(yhole/xhole)),yhole-d1*Sin(ArcTan(yhole/xhole)),0,lsing};
Point(62) = {-xhole,-yhole,0,lsing};
Point(1062) = {-xhole+d1*Cos(ArcTan(yhole/xhole)),-yhole+d1*Sin(ArcTan(yhole/xhole)),0,lsing};
Point(63) = {-xhole2,yhole2,0,lsing};
Point(64) = {-xhole2,-yhole2,0,lsing};

Line(101) = {43,45};
Line(102) = {44,46};

Line(10101) = {1043,1045};
Line(10102) = {1044,1046};

Circle(103) = {43,1,44};
Circle(104) = {45,47,61};
Circle(105) = {62,47,46};

Circle(10103) = {1043,1,1044};
Circle(10104) = {1045,47,1061};
Circle(10105) = {1062,47,1046};

Point(111) = {-0.003958,ry2,0,lcooling};
Point(121) = {-0.003958,-ry2,0,lcooling};

Line(111) = {53,111};
Line(118) = {111,55};
Line(112) = {54,121};
Line(119) = {121,56};

Circle(113) = {53,57,54};
Circle(114) = {55,58,63};
Circle(115) = {64,58,56};

Line(116) = {61,63};
Line(117) = {62,64};

Line(10116) = {1061,61};
Line(10117) = {1062,62};

// pipe loop
Line Loop(301) = {-119,-112,-113,111,118,114,-116,-104,-101,103,102,-105,117,115};
// coating loop
Line Loop(10301) = {-102,-103,101,104,-10116,-10104,-10101,10103,10102,-10105,10117,105};

// outside vacuum
Point(71) = {h,0, 0, lvac};
Point(72) = {0,h, 0, lvac};
Point(73) = {-h,0, 0, lvac};
Point(74) = {0,-h, 0, lvac};

Circle(201) = {71,1,72};
Circle(202) = {72,1,73};
Circle(203) = {73,1,74};
Circle(204) = {74,1,71};

Line Loop(302) = {201,202,203,204};

// cooling
Point(101) = {-0.017,0,0,lsing};
Point(102) = {-0.0187,0.0033,0,lcooling};
Point(103) = {-0.005826,0.01807,0,lcooling};
Point(104) = {-0.005125,0.01933,0,lcooling};
Point(105) = {-0.005125,-0.01933,0,lcooling};
Point(106) = {-0.005826,-0.01807,0,lcooling};
Point(107) = {-0.0187,-0.0033,0,lcooling};

Point(110) = {-0.003958,0.016463,0,lcooling};
// Point(111) see above
Point(112) = {-0.003958,2*0.016463-ry2,0,lcooling};

Point(120) = {-0.003958,-0.016463,0,lcooling};
// Point(121) see above
Point(122) = {-0.003958,-2*0.016463+ry2,0,lcooling};

Line(300) = {101,102};
Circle(301) = {102,1,103};
Circle(302) = {104,1,105};
Circle(303) = {106,1,107};
Line(304) = {107,101};

Circle(305) = {103,110,111};
Circle(306) = {111,110,112};
Line(307) = {112,104};

Circle(308) = {122,120,121};
Circle(309) = {121,120,106};
Line(310) = {105,122};

Line Loop(303) = {300,301,305,306,307,302,310,308,309,303,304};

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
Plane Surface(13) = {301};

// coating
Plane Surface(1013) = {10301};

// vacuum
Line Loop(304)={-10102,-10103,10101,10104,10116,116,-114,-118,-305,-301,-300,-304,-303,-309,119,-115,-117,-10117,10105};
Line Loop(305)={-112,-113,111,306,307,302,310,308};
Plane Surface(14)={304,103,113,123,133};
Plane Surface(15)={302,305};

// cooling
Plane Surface(16)={303};

//Physical surfaces

Physical Surface(1) = {2,3,5,6,8,9,11,12};   //beam
Physical Surface(2) = {1,4,7,10,14,15};      //vacuum
Physical Surface(5) = {13};                  //pipe
Physical Surface(3) = {16};                  //Reflector
Physical Surface(4) = {1013};		//Coating


Color Grey50{ Surface{ 13 }; }
Color Purple{ Surface{ 16 }; }
Color Green{ Surface{ 14:15 }; }



