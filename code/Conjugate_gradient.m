function out = Conjugate_gradient(n,u,v,p,f,g,tol)
% 共轭梯度法求解方程
if nargin < 7
    tol = 1e-8;
end
%% 初始化
h = 1/n;
f(2:n,:) = f(2:n,:)*h^2-(p(2:n,:)-p(1:n-1,:))*h;
g(:,2:n) = g(:,2:n)*h^2-(p(:,2:n)-p(:,1:n-1))*h;

% r0 = norm(f(:))^2+norm(g(:))^2;
% r0 = sqrt(r0);
% rho = r0;
% rf = zeros(n+1,n);
% rg = zeros(n,n+1);
ite = 0;


pu = zeros(n+1,n);
pv = zeros(n,n+1);
pu(2:n,1) = f(2:n,1)+u(1:n-1,1)+u(3:n+1,1)+u(2:n,2)-u(2:n,1)*3;
pu(2:n,2:n-1) = f(2:n,2:n-1)+u(1:n-1,2:n-1)+u(3:n+1,2:n-1)+u(2:n,1:n-2)+u(2:n,3:n)-u(2:n,2:n-1)*4;
pu(2:n,n) = f(2:n,n)+u(1:n-1,n)+u(3:n+1,n)+u(2:n,n-1)-u(2:n,n)*3;
pv(1,2:n) = g(1,2:n)+v(1,1:n-1)+v(1,3:n+1)+v(2,2:n)-v(1,2:n)*3;
pv(2:n-1,2:n) = g(2:n-1,2:n)+v(2:n-1,1:n-1)+v(2:n-1,3:n+1)+v(1:n-2,2:n)+v(3:n,2:n)-v(2:n-1,2:n)*4;
pv(n,2:n) = g(n,2:n)+v(n,1:n-1)+v(n,3:n+1)+v(n-1,2:n)-v(n,2:n)*3;
% u = zeros(n+1,n);
% v = zeros(n,n+1);
ru = pu;
rv = pv;
rho_u0 = norm(pu(:))^2;
rho_v0 = norm(pv(:))^2;
rho_u = rho_u0;
rho_v = rho_v0;

wu = zeros(n+1,n);
wv = zeros(n,n+1);
maxn = 3;
%% 共轭梯度法
while sqrt(rho_u) > tol*sqrt(rho_u0) && ite < maxn %&& sqrt(rho_u) > tol
    ite = ite+1;
    if ite ~= 1
        betau = rho_u/rho_u1;
        pu = ru+betau*pu;        
    end
    wu(2:n,1) = 3*pu(2:n,1)-pu(1:n-1,1)-pu(3:n+1,1)-pu(2:n,2);
    wu(2:n,2:n-1) = 4*pu(2:n,2:n-1)-pu(1:n-1,2:n-1)-pu(3:n+1,2:n-1)-pu(2:n,1:n-2)-pu(2:n,3:n);
    wu(2:n,n) = 3*pu(2:n,n)-pu(1:n-1,n)-pu(3:n+1,n)-pu(2:n,n-1);
    alpha = rho_u/(pu(:)'*wu(:));
    if alpha > maxn
        break;
    end
    u = u+alpha*pu;
    ru = ru-alpha*wu;
    rho_u1 = rho_u;
    rho_u = norm(ru(:))^2;

end
ite = 0;
while sqrt(rho_v) > tol*sqrt(rho_v0)  && ite < maxn %&& sqrt(rho_v) > tol
    ite = ite+1;
    if ite ~= 1
        betav = rho_v/rho_v1; 
        pv = rv+betav*pv;
    end
    wv(1,2:n) = 3*pv(1,2:n)-pv(1,1:n-1)-pv(1,3:n+1)-pv(2,2:n);
    wv(2:n-1,2:n) = 4*pv(2:n-1,2:n)-pv(2:n-1,1:n-1)-pv(2:n-1,3:n+1)-pv(1:n-2,2:n)-pv(3:n,2:n);
    wv(n,2:n) = 3*pv(n,2:n)-pv(n,1:n-1)-pv(n,3:n+1)-pv(n-1,2:n);
    alpha = rho_v/(pv(:)'*wv(:));
    if alpha > maxn
        break;
    end
    v = v+alpha*pv;
    rv = rv-alpha*wv;
    rho_v1 = rho_v;
    rho_v = norm(rv(:))^2;
end
%% 输出
out.u = u;
out.v = v;