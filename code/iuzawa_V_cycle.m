function out = iuzawa_V_cycle(n,opt)
%大小为n的DCG多重循环网格,提升到n/L
if nargin < 2
    opt.tol = 1e-8;
    opt.n1 = 1;
    opt.n2 = 1; 
    opt.L = n/2;
    opt.tau = 1e-2;
    opt.alpha = 1;
end
% if ~isfield("L",opt);   opt.L = n/4;    end
% if ~isfield("n1",opt);  opt.n1 = 2;     end
% if ~isfield("n2",opt);  opt.n2 = 2;     end
% if ~isfield("tol",opt); opt.tol = 1e-8; end
%% 初始化
tic;

h = 1/n;
tol = opt.tol;
init = init_(n);
f = init.f;
g = init.g;
u = zeros(n+1,n);
u([1,n+1],:) = u([1,n+1],:)*0;
v = zeros(n,n+1);
v(:,[1,n+1]) = v(:,[1,n+1])*0;
p = zeros(n);


r0 = residual_(n,u,v,p,f,g,zeros(n));
rf0 = r0.rf;
rg0 = r0.rg;
rd0 = r0.rd;
res0 = sqrt(norm(rf0(:))^2+norm(rg0(:))^2)+norm(rd0(:))^2;
res0 = sqrt(res0);

%% 开始迭代
res = res0;
cnt = 0;
rd = zeros(n);
while res/res0 > tol && cnt < 40000
    cnt = cnt+1;
    out = V_cycle_iuzawa(n,p,f,g,opt,rd,u,v);
    u = out.u;
    v = out.v;
    p = out.p;
    res = out.res;
    res = sqrt(norm(res.rf(:))^2+norm(res.rg(:))^2+norm(res.rd(:))^2);

    fprintf("第%d次多重循环网格迭代,误差为%e\n",cnt,res/res0);
end
time = toc;
x = (1-cos(2*pi*(0:h:1)))'*sin(2*pi*(h/2:h:1));
y = -sin(2*pi*(h/2:h:1))'*(1-cos(2*pi*(0:h:1)));
rx = x-u;
ry = y-v;
e = h*sqrt(norm(rx(:))^2+norm(ry(:))^2);
fprintf("迭代次数为%d,最终的相对误差为%e,e0=%e,花的时间为%f\n",cnt,res/res0,e,time);
% disp("最终迭代次数为%d,误差为%e\n, 花的",cnt,res/res0);
out.u = u;
out.v = v;
out.p = p;
out.e = e;
out.ite = cnt;
out.res = res;
out.time = toc;