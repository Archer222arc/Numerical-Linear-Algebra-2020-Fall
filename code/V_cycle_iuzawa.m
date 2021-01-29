function out = V_cycle_iuzawa(n,p,f,g,opt,d,u,v)
%V-cycle
if ~isfield(opt,'n1'), opt.n1 = 1;   end %粗网格
if ~isfield(opt,'n2'), opt.n2 = 1;   end %细网格

fprintf("在对大小 %d 的网格进行迭代\n",n);
h = 1/n;
if nargin < 8
    u = zeros(n+1,n);
    v = zeros(n,n+1);
end
if opt.L == 1
    tmp = iuzawa(n,u,v,p,opt.alpha,opt.tol,f,g,d,opt.tau);
    out.u = tmp.u;
    out.v = tmp.v;
    out.p = tmp.p;
    out.res = tmp.res;
    return;
end

n1 = opt.n1;    n2 = opt.n2;
opt.L = opt.L/2;

%% 在细网格迭代
for i = 1:n1
    tmp = Conjugate_gradient(n,u,v,p,f,g,opt.tau);
    u = tmp.u;
    v = tmp.v;
    rd = (u(2:n+1,:)-u(1:n,:)+v(:,2:n+1)-v(:,1:n))/h-d;
    p = p-opt.alpha*rd;
%     res = residual_(n,u,v,p,f,g,d);
%     rd = tmp.res.rd;
end

if opt.L == 1
    out.u = u;
    out.v = v;
    out.p = p;
    return
end


res = residual_(n,u,v,p,f,g,d);
rf = res.rf;
rg = res.rg;
rd = res.rd;

%% 限制



tmp = down(n,rf,rg,rd);
% tmp = V_cycle(n/2,zeros(n/2),tmp.rf,tmp.rg,opt,tmp.rd);
tmp = V_cycle_iuzawa(n/2,zeros(n/2),tmp.rf,tmp.rg,opt,tmp.rd);


%% 提升回细网格
tmp = lift(n/2,tmp.u,tmp.v,tmp.p);

u = u+tmp.u;
v = v+tmp.v;
p = p+tmp.p;


for i = 1:n2
    tmp = Conjugate_gradient(n,u,v,p,f,g,opt.tau);
    u = tmp.u;
    v = tmp.v;
    rd = (u(2:n+1,:)-u(1:n,:)+v(:,2:n+1)-v(:,1:n))/h-d;
    p = p-opt.alpha*rd;
%     rd = tmp.res.rd;
end
res = residual_(n,u,v,p,f,g,d);



%% 返回

out.u = u;
out.v = v;
out.p = p;
out.res = res;