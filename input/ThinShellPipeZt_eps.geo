lveryfine=8e-4;
lfine = 1e-3;

//lrough = 5e-4;
lpipe=2e-4;
lvac=1e-2;


a=0.01;
b=0.04;
h=0.041;
h2=0.1;
eps=0.000005;



am=0.025;
lm=5e-3;






Point(1) = {0, 0, 0, lfine};

Point(2) = {a-eps, 0, 0, lveryfine};
Point(3) = {a,0, 0, lveryfine};
Point(4) = {a+eps,0, 0, 1.0*lveryfine};
Point(5) = {am,0, 0, lm};
Point(6) = {b,0, 0, lpipe};

Point(12) = {0,a-eps, 0, lveryfine};
Point(13) = {0,a,0,  lveryfine};
Point(14) = {0,a+eps, 0, 1.0*lveryfine};
Point(15) = {0,am, 0, lm};
Point(16) = {0,b, 0, lpipe};

Point(22) = {-(a-eps), 0, 0, lveryfine};
Point(23) = {-a,0, 0, lveryfine};
Point(24) = {-(a+eps),0, 0, 1.0*lveryfine};
Point(25) = {-am,0, 0, lm};
Point(26) = {-b,0, 0, lpipe};

Point(32) = {0,-(a-eps), 0, lveryfine};
Point(33) = {0,-a,0,  lveryfine};
Point(34) = {0,-(a+eps), 0, 1.0*lveryfine};
Point(35) = {0,-am, 0, lm};
Point(36) = {0,-b, 0, lpipe};


Point(221) = {h,0, 0, lpipe};
Point(222) = {0,h, 0, lpipe};
Point(223) = {-h,0, 0, lpipe};
Point(224) = {0,-h, 0, lpipe};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,221};

Line(11) = {1,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,222};

Line(21) = {1,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,25};
Line(25) = {25,26};
Line(26) = {26,223};

Line(31) = {1,32};
Line(32) = {32,33};
Line(33) = {33,34};
Line(34) = {34,35};
Line(35) = {35,36};
Line(36) = {36,224};



Circle(51) = {2,1,12};
Circle(52) = {3,1,13};
Circle(53) = {4,1,14};
Circle(54) = {5,1,15};
Circle(55) = {6,1,16};


Circle(61) = {12,1,22};
Circle(62) = {13,1,23};
Circle(63) = {14,1,24};
Circle(64) = {15,1,25};
Circle(65) = {16,1,26};

Circle(71) = {22,1,32};
Circle(72) = {23,1,33};
Circle(73) = {24,1,34};
Circle(74) = {25,1,35};
Circle(75) = {26,1,36};


Circle(81) = {32,1,2};
Circle(82) = {33,1,3};
Circle(83) = {34,1,4};
Circle(84) = {35,1,5};
Circle(85) = {36,1,6};


//mesh
//Circle(91)={6,1,16};
//Circle(92)={16,1,26};
//Circle(93)={26,1,36};
//Circle(94)={36,1,6};




//Line Loop sectors
Line Loop(101)={1,51,-11};
Line Loop(102)={2,52,-12,-51};
Line Loop(103)={3,53,-13,-52};
Line Loop(104)={4,54,-14,-53};
Line Loop(105)={5,55,-15,-54};

Line Loop(111)={11,61,-21};
Line Loop(112)={12,62,-22,-61};
Line Loop(113)={13,63,-23,-62};
Line Loop(114)={14,64,-24,-63};
Line Loop(115)={15,65,-25,-64};

Line Loop(121)={21,71,-31};
Line Loop(122)={22,72,-32,-71};
Line Loop(123)={23,73,-33,-72};
Line Loop(124)={24,74,-34,-73};
Line Loop(125)={25,75,-35,-74};

Line Loop(131)={31,81,-1};
Line Loop(132)={32,82,-2,-81};
Line Loop(133)={33,83,-3,-82};
Line Loop(134)={34,84,-4,-83};
Line Loop(135)={35,85,-5,-84};


//Line Loop circles
Line Loop(151)={51,61,71,81};
Line Loop(152)={52,62,72,82};
Line Loop(153)={53,63,73,83};
Line Loop(154)={54,64,74,84};
Line Loop(155)={55,65,75,85};





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
Plane Surface(41)={105};

Plane Surface(5) = {111};
Plane Surface(6) = {112};
Plane Surface(7) = {113};
Plane Surface(8) = {114};
Plane Surface(81)={115};

Plane Surface(9) = {121};
Plane Surface(10) = {122};
Plane Surface(11) = {123};
Plane Surface(12) = {124};
Plane Surface(121)={125};

Plane Surface(13) = {131};
Plane Surface(14) = {132};
Plane Surface(15) = {133};
Plane Surface(16) = {134};
Plane Surface(161)={135};

//Line {1} In Surface {30};






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
//The pipe



Circle(131) = {221,1,222};
Circle(132) = {222,1,223};
Circle(133) = {223,1,224};
Circle(134) = {224,1,221};

Line Loop(32)={131,132,133,134};
Plane Surface(17) = {32,105,115,125,135};   //Pipe


//The outside vacuum
Point(331) = {h2,0, 0, lvac};
Point(332) = {0,h2, 0, lvac};
Point(333) = {-h2,0, 0, lvac};
Point(334) = {0,-h2, 0, lvac};


Circle(41) = {331,1,332};
Circle(42) = {332,1,333};
Circle(43) = {333,1,334};
Circle(44) = {334,1,331};

Line Loop(33)={41,42,43,44};
Plane Surface(18) = {33,32,105,104};


//Physical Surface(4) = {18};   //OutsideVac

Physical Surface(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13};   //beam   (+/-eps)
Physical Surface(2) = {14,15,16,18,41,81,121,161};   //vacuum
Physical Surface(3) = {17};   //Pipe






