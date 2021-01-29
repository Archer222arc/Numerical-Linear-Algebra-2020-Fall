function out = V_cycle(n,p,f,g,opt,d,u,v)
%V-cycle
if ~isfield(opt,'n1'), opt.n1 = 1;   end %粗网格
if ~isfield(opt,'n2'), opt.n2 = 1;   end %细网格

fprintf("在对大小 %d 的网格进行迭代\n",n);

if opt.L == 1
    u = zeros(n+1,n);
    v = zeros(n,n+1);
    res = residual_(n,u,v,p,f,g,d);
    r = norm(res.rf(:))^2+norm(res.rg(:))^2+norm(res.rd(:))^2;
    r = sqrt(r);
    while r > 1e-8
        tmp = DGS(n,u,v,p,f,g,d);
        u = tmp.u;
        v = tmp.v;
        p = tmp.p;
        res = residual_(n,u,v,p,f,g,d);
        r = norm(res.rf(:))^2+norm(res.rg(:))^2+norm(res.rd(:))^2;
        r = sqrt(r);
    end
    out.u = u;
    out.v = v;
    out.p = p;
    out.res = res;
    return;
end

n1 = opt.n1;    n2 = opt.n2;
opt.L = opt.L/2;
if nargin < 8
    u = zeros(n+1,n);
    v = zeros(n,n+1);
end
%% 在细网格迭代
for i = 1:n1
    tmp = DGS(n,u,v,p,f,g,d);
    u = tmp.u;
    v = tmp.v;
    p = tmp.p;
%     rd = tmp.res.rd;
end

res = residual_(n,u,v,p,f,g,d);
rf = res.rf;
rg = res.rg;
rd = res.rd;

%% 限制



tmp = down(n,rf,rg,rd);
% tmp = V_cycle(n/2,zeros(n/2),tmp.rf,tmp.rg,opt,tmp.rd);
tmp = V_cycle(n/2,zeros(n/2),tmp.rf,tmp.rg,opt,tmp.rd);


%% 提升回细网格
tmp = lift(n/2,tmp.u,tmp.v,tmp.p);

u = u+tmp.u;
v = v+tmp.v;
p = p+tmp.p;


for i = 1:n2
    tmp = DGS(n,u,v,p,f,g,d);
    u = tmp.u;
    v = tmp.v;
    p = tmp.p;
%     rd = tmp.res.rd;
end

res = tmp.res;


%% 返回

out.u = u;
out.v = v;
out.p = p;
out.res = res;