function out = init_(n)
%%给出初始的f,g
%采用了跟matlab矩阵下标一致的写法，即（1，1）从左上开始，分别为纵坐标和横坐标
d = 1/n;
%% F
x = (0:d:1)';
y = d/2:d:1;
f = -(2*cos(2*pi*x)-1)*sin(2*pi*y)*4*pi^2+(x.^2)*ones(1,n);
f(1:n+1,1) = f(1:n+1,1)-(1-cos(2*pi*x))*2*pi/d;
f(1:n+1,n) = f(1:n+1,n)+(1-cos(2*pi*x))*2*pi/d;
out.f = f;

%% G
x = (d/2:d:1)';
y = 0:d:1;
g = 4*pi^2*sin(2*pi*x)*(2*cos(2*pi*y)-1);
g(1,1:n+1) = g(1,1:n+1)+(1-cos(2*pi*y))*2*pi/d;
g(n,1:n+1) = g(n,1:n+1)-(1-cos(2*pi*y))*2*pi/d;
out.g=g;
