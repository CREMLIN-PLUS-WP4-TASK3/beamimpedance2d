a=0.001;
halfgap=3*a;
halfheight=0.03;
xm=0.05;
ym=0.05;
thick=xm-halfgap;





lveryfine=5e-5;
lfine = 5e-4;
lrough = 5e-3;



///////////////////////////////////////////////////////////
//The beam
eps=0.000005;



Point(1) = {0, 0, 0, lfine};

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
/////////////////////////////////////////////////////////

//////////////////////////////////////////
//bounding box
Point(101)={xm,ym,0,lrough};
Point(102)={-xm,ym,0,lrough};
Point(103)={-xm,halfheight,0,lrough};
Point(104)={-xm,-halfheight,0,lrough};
Point(105)={-xm,-ym,0,lrough};
Point(106)={xm,-ym,0,lrough};
Point(107)={xm,-halfheight,0,lrough};
Point(108)={xm,halfheight,0,lrough};


Line(101)={101,102};
Line(102)={102,103};
//Line(103)={103,104};
Line(104)={104,105};
Line(105)={105,106};
Line(106)={106,107};
//Line(107)={107,108};
Line(108)={108,101};
/////////////////////////////////////////////////


/////////////////////////////////////////////////////////
//Collimator
Point(201)={-halfgap,halfheight,0,lfine};
Point(202)={-halfgap,-halfheight,0,lfine};
Point(203)={halfgap,-halfheight,0,lfine};
Point(204)={halfgap,halfheight,0,lfine};

Line(111)={103,201};
Line(112)={201,202};
Line(113)={202,104};

Line(114)={108,204};
Line(115)={204,203};
Line(116)={203,107};

////////////////////////////////////////////////////////


//Line Loop sectors (beam)
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


// Line loop vacuum
Line Loop(201)={101,102,111,112,113,104,105,106,-116,-115,-114,108};

// Line loops collimator
//Line Loop(301)={103,-113,-112,-111};
//Line Loop(302)={107,114,115,116};


//beam and inner vac
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


//vac
Plane Surface(13)={201,103,113,123,133};  //the beam is a hole in the vacuum :-)

//collimator jaws
//Plane Surface(14)={301};
//Plane Surface(15)={302};


Physical Surface(1) = {2,3,5,6,8,9,11,12};   //beam   (+/-eps)
Physical Surface(2) = {1,4,7,10,13};   //vacuum
//Physical Surface(3) = {14,15}; //jaws

