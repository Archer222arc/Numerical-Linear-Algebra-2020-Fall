function out = uzawa(n,u,v,p,alpha,tol)
% uzawa，输入n,u,v,alpha和tol，默认为zeros,0.3,1e-8
tic;
if nargin < 6
    u = zeros(n+1,n);
    u(1,:) = u(1,:)*0;
    u(n+1,:) = u(n+1,:)*0;
    v = zeros(n,n+1);
    v(:,1) = v(:,1)*0;
    v(:,n+1) = v(:,n+1)*0;
    tol = 1e-8;
    alpha = 1;
    p = zeros(n);
end
%%初始化
h = 1/n;
init = init_(n);
f = init.f;
g = init.g;
d = zeros(n);

r0 = residual_(n,u,v,p,f,g,zeros(n));
rf0 = r0.rf;
rg0 = r0.rg;
rd0 = r0.rd;
res0 = sqrt(norm(rf0(:))^2+norm(rg0(:))^2)+norm(rd0(:))^2;
res0 = sqrt(res0);
res = res0;
cnt = 0;
%% uzawa迭代
while res/res0 > tol
    cnt = cnt+1;

    tmp = Conjugate_gradient(n,u,v,p,f,g,1e-8);
    u = tmp.u;
    v = tmp.v;
    rd = (u(2:n+1,:)-u(1:n,:)+v(:,2:n+1)-v(:,1:n))/h;
    p = p-alpha*rd;
    res = residual_(n,u,v,p,f,g,d);
    res = norm(res.rf(:))^2+norm(res.rg(:))^2+norm(res.rd(:))^2;
    res = sqrt(res);
%     fprintf("迭代次数为%d相对误差为%e\n",cnt,res/res0);
end

%% 输出
time = toc;
x = (1-cos(2*pi*(0:h:1)))'*sin(2*pi*(h/2:h:1));
y = -sin(2*pi*(h/2:h:1))'*(1-cos(2*pi*(0:h:1)));
rx = x-u;
ry = y-v;
e = h*sqrt(norm(rx(:))^2+norm(ry(:))^2);
fprintf("迭代次数为%d,最终的相对误差为%e,e0=%e,花的时间为%f\n",cnt,res/res0,e,time);
out.res = res;
out.u = u;
out.v = v;
out.p = p;
out.e = e;
out.time = time;