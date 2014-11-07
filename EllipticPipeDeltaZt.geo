lveryfine=8e-4;  //4e-4
lfine = 5e-4;
lrough = 4e-4;

a=0.01;
bx=0.073;
by=0.035;
eps=0.000005;

//r_h=0.073;
//r_v=0.035;

Point(1) = {0, 0, 0, lfine};

Point(2) = {a-eps, 0, 0, lveryfine};
Point(3) = {a,0, 0, lveryfine};
Point(4) = {a+eps,0, 0, 1.0*lveryfine};
Point(5) = {bx,0, 0, lrough};				//Outside

Point(12) = {0,a-eps, 0, lveryfine};
Point(13) = {0,a,0,  lveryfine};
Point(14) = {0,a+eps, 0, 1.0*lveryfine};
Point(15) = {0,by, 0, lrough};				//Outside

Point(22) = {-(a-eps), 0, 0, lveryfine};
Point(23) = {-a,0, 0, lveryfine};
Point(24) = {-(a+eps),0, 0, 1.0*lveryfine};
Point(25) = {-bx,0, 0, lrough};				//Outside

Point(32) = {0,-(a-eps), 0, lveryfine};
Point(33) = {0,-a,0,  lveryfine};
Point(34) = {0,-(a+eps), 0, 1.0*lveryfine};
Point(35) = {0,-by, 0, lrough};				//Outside



Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};  

Line(11) = {1,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};

Line(21) = {1,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,25};

Line(31) = {1,32};
Line(32) = {32,33};
Line(33) = {33,34};
Line(34) = {34,35};




Circle(51) = {2,1,12};
Circle(52) = {3,1,13};
Circle(53) = {4,1,14};
//Circle(54) = {5,1,15};   //Outside
Ellipse(54) = {5, 1, 2, 15};

Circle(61) = {12,1,22};
Circle(62) = {13,1,23};
Circle(63) = {14,1,24};
//Circle(64) = {15,1,25};  //Outside
Ellipse(64) = {15, 1, 2, 25};

Circle(71) = {22,1,32};
Circle(72) = {23,1,33};
Circle(73) = {24,1,34};
//Circle(74) = {25,1,35};  //Outside
Ellipse(74) = {25, 1, 2, 35};

Circle(81) = {32,1,2};
Circle(82) = {33,1,3};
Circle(83) = {34,1,4};
//Circle(84) = {35,1,5};  //Outside
Ellipse(84) = {35, 1, 2, 5};

//Line Loop sectors
Line Loop(101)={1,51,-11};
Line Loop(102)={2,52,-12,-51};
Line Loop(103)={3,53,-13,-52};
Line Loop(104)={4,54,-14,-53};

Line Loop(111)={11,61,-21};
Line Loop(112)={12,62,-22,-61};
Line Loop(113)={13,63,-23,-62};
Line Loop(114)={14,64,-24,-63};

Line Loop(121)={21,71,-31};
Line Loop(122)={22,72,-32,-71};
Line Loop(123)={23,73,-33,-72};
Line Loop(124)={24,74,-34,-73};

Line Loop(131)={31,81,-1};
Line Loop(132)={32,82,-2,-81};
Line Loop(133)={33,83,-3,-82};
Line Loop(134)={34,84,-4,-83};


//Line Loop circles
Line Loop(151)={51,61,71,81};
Line Loop(152)={52,62,72,82};
Line Loop(153)={53,63,73,83};
Line Loop(154)={54,64,74,84};



/*
Plane Surface(1) = {151};    //a -eps
Plane Surface(2) = {152,101,102,103,104};  //a 
Plane Surface(3) = {153,101,102,103,104,111,112,113,114}; //a+eps
Plane Surface(4) = {154,101,102,103,104,111,112,113,114,121,122,123,124};
doesn't work
*/


Plane Surface(1) = {101};
Plane Surface(2) = {102};
Plane Surface(3) = {103};
Plane Surface(4) = {104};
Plane Surface(5) = {111};
Plane Surface(6) = {112};
Plane Surface(7) = {113};
Plane Surface(8) = {114};
Plane Surface(9) = {121};
Plane Surface(10) = {122};
Plane Surface(11) = {123};
Plane Surface(12) = {124};
Plane Surface(13) = {131};
Plane Surface(14) = {132};
Plane Surface(15) = {133};
Plane Surface(16) = {134};

//Line {1} In Surface {30};


Physical Surface(1) = {1,2,3,4,5,6,7,8,9,10,11,12};   //beam   (+/-eps)
//Physical Surface(2) = {2};
//Physical Surface(3) = {3};
//Physical Surface(4) = {4};
Physical Surface(2) = {13,14,15,16};   //vacuum

//Physical Line(6)={1};





