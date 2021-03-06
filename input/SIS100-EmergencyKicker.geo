a=0.005;
halfgap=3*a;
halfheight=0.03;
xm=0.15;
ym=0.15;

FT=0.06;      // Ferrite Thickness
KHH=0.11;    // Kicker Half Height
KHW=0.1275;  // Kicker Half Width
HG=0.0014;    // Half Gap

r1=0.0178;
r2=0.0305;
r3=0.033;



lveryfine=1e-3;
lfine = 2e-2;
lrough = 1e-2;
lCopper = 1e-3;
lFerrite = 1e-2;
lOutsideferrite=1e-2;

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
Line(103)={103,104};
Line(104)={104,105};
Line(105)={105,106};
Line(106)={106,107};
Line(107)={107,108};
Line(108)={108,101};
/////////////////////////////////////////////////


/////////////////////////////////////////////////////////
//Upper Yoke
Point(201)={KHW,HG,0,lCopper};
Point(202)={KHW,KHH,0,lOutsideferrite};
Point(203)={-KHW,KHH,0,lOutsideferrite};
Point(204)={-KHW,HG,0,lCopper};

Point(205)={-KHW+FT,HG,0,lCopper};
Point(206)={-KHW+FT,KHH-FT,0,lFerrite};
Point(207)={KHW-FT,KHH-FT,0,lFerrite};
Point(208)={KHW-FT,HG,0,lCopper};

Line(201)={201,202};
Line(202)={202,203};
Line(203)={203,204};
Line(204)={204,205};
Line(205)={205,206};
Line(206)={206,207};
Line(207)={207,208};
Line(208)={208,201};

//Lower Yoke
Point(301)={KHW,-HG,0,lCopper};
Point(302)={KHW,-KHH,0,lOutsideferrite};
Point(303)={-KHW,-KHH,0,lOutsideferrite};
Point(304)={-KHW,-HG,0,lCopper};

Point(305)={-KHW+FT,-HG,0,lCopper};
Point(306)={-KHW+FT,-KHH+FT,0,lFerrite};
Point(307)={KHW-FT,-KHH+FT,0,lFerrite};
Point(308)={KHW-FT,-HG,0,lCopper};

Line(301)={301,302};
Line(302)={302,303};
Line(303)={303,304};
Line(304)={304,305};
Line(305)={305,306};
Line(306)={306,307};
Line(307)={307,308};
Line(308)={308,301};


//Right Gap  -->Lines oriented down
Line(401)={201,301};
Line(402)={208,308};

//Left Gap  -->Lines oriented down
Line(411)={204,304};
Line(412)={205,305};



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


// Line loop outside vacuum
Line Loop(201)={101,102,103,104,105,106,107,108};

//Line Loop outside Kicker
//Line Loop(202)={201,202.203, 411, -303,-302,-301, -401};



// Line loops upper yoke
Line Loop(301)={201,202,203,204,205,206,207,208};   

//Line loop lower yoke
Line Loop(302)={-301,-302,-303,-304,-305,-306,-307,-308};      


//Line loop right gap
Line Loop(303)={-401,-208,402,308};   

//Line loop left gap
Line Loop(304)={411,304,-412,-204};   

//Line loop inner vacuum
Line Loop(500)={305,306,307,-402, -207,-206,-205, 412};



///////////////////////////////////////////////////////////////////////



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
Plane Surface(13)={500,103,113,123,133};  //inside -->the beam is a hole in the vacuum :-)
Plane Surface(14)={201,301,302,303,304};    //outside   -->the kicker is a hole in the vacuum :-)


//Upper yoke
Plane Surface(15)={301};

//Lower yoke
Plane Surface(16)={302};

//Right gap
Plane Surface(17)={303};

//Left gap
Plane Surface(18)={304};



Physical Surface(1) = {2,3,5,6,8,9,11,12};   //beam   (+/-eps)
Physical Surface(2) = {1,4,7,10,13,14};   //vacuum
Physical Surface(5) = {15,16}; //Ferrite
Physical Surface(4) = {17,18}; //Copper


Color Grey50{ Surface{ 13:14 }; }
//Color Purple{ Surface{ 2:3:5:6:8:9:11:12 }; }
Color Blue{ Surface{ 15:16 }; }
Color Green{ Surface{ 17:18 }; }


