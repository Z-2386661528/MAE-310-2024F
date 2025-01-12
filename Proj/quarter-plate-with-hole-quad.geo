//定义了两个变量，方便后续画表格
R = 0.5;
L = 1.0;

//找到所有需要的点，按一定顺序编排，1-4在正方形四个端点
Point(1) = {L, -L, 0};
Point(2) = {L, L, 0};
Point(3) = {-L, L, 0};
Point(4) = {-L, -L, 0};

//567分别是圆上三个点，
Point(5) = {-L + R, -L, 0};
Point(6) = {-L, -L + R, 0};
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};

//使用圆弧连线，4是圆心
Circle(1) = {5, 4, 7};
Circle(2) = {7, 4, 6};

//继续连线，绘制整个图
Line(3) = {6, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};
Line(7) = {2, 7};

//画出第一个闭合图形，由4，7，2，3构成，位置在左上，并通过这条曲线定义第一个平面
Curve Loop(1) = {-7, -4, -3, -2};
Plane Surface(1) = {1};

//由7， -1， -6， -5构成，位置在右下，并通过这条曲线定义第二个平面

Curve Loop(2) = {7, -1, -6, -5};
Plane Surface(2) = {2};

//有限元，每条线分段为2
Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 80;

//对平面12表面按无边界的方法平均画网格
Transfinite Surface{1};
Transfinite Surface{2};

//对网格优化
Recombine Surface{1};
Recombine Surface{2};

//1：线性生成单元，8：网格生成算法
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;

//+
Physical Curve("right", 8) = {5};
//+
Physical Curve("up", 9) = {4};
//+
Physical Curve("left", 10) = {3};
//+
Physical Curve("down", 11) = {6};
//+
Physical Curve("circle", 12) = {1, 2};
//+
Physical Surface("surf", 13) = {2, 1};
