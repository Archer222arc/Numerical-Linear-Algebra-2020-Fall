function out = iuzawa(n,u,v,p,alpha,tol,f,g,d,tau)
% iuzawa，输入n,u,v,alpha和tol，默认为zeros,1,1e-8
if nargin < 6
    u = zeros(n+1,n);
    u(1,:) = u(1,:)*0;
    u(n+1,:) = u(n+1,:)*0;
    v = zeros(n,n+1);
    v(:,1) = v(:,1)*0;
    v(:,n+1) = v(:,n+1)*0;
    tol = 1e-2;
    alpha = 1;
    p = zeros(n);
end
%%初始化
h = 1/n;
if nargin < 9
    init = init_(n);
    f = init.f;
    g = init.g;
    d = zeros(n);
end


r0 = residual_(n,u,v,p,f,g,zeros(n));
rf0 = r0.rf;
rg0 = r0.rg;
rd0 = r0.rd;
res0 = sqrt(norm(rf0(:))^2+norm(rg0(:))^2)+norm(rd0(:))^2;
res0 = sqrt(res0);
res = res0;
cnt = 0;
%% iuzawa迭代
while res/res0 > tol
    cnt = cnt+1;

    tmp = Conjugate_gradient(n,u,v,p,f,g,tau);
    u = tmp.u;
    v = tmp.v;
    rd = (u(2:n+1,:)-u(1:n,:)+v(:,2:n+1)-v(:,1:n))/h-d;
    p = p-alpha*rd;
    res = residual_(n,u,v,p,f,g,d);
    res = norm(res.rf(:))^2+norm(res.rg(:))^2+norm(rd(:))^2;
    res = sqrt(res);
end

%% 输出
out.res = res;
out.u = u;
out.v = v;
out.p = p;
