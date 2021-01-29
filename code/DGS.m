function out = DGS(n,u,v,p,f,g,d)
%单层DCG
if nargin < 2
    u = zeros(n+1,n);
    v = zeros(n,n+1);
    p = zeros(n,n);
    init = init_(n);
    f = init.f;
    g = init.g;
end

h = 1/n;

%% 更新 u 
for i = 2:n
    u(i,1) = (h^2*f(i,1)-h*(p(i,1)-p(i-1,1))+u(1+i,1)+u(i-1,1)+u(i,2))/3;
end

for j = 2:n-1
    for i = 2:n
        u(i,j) = (h^2*f(i,j)-h*(p(i,j)-p(i-1,j))+u(i-1,j)+u(i,j+1)+u(i+1,j)+u(i,j-1))/4;
    end
end

for i = 2:n
    u(i,n) = (h^2*f(i,n)-h*(p(i,n)-p(i-1,n))+u(i+1,n)+u(i-1,n)+u(i,n-1))/3;
end
%% 更新 v
for i = 2:n
    v(1,i) = (h^2*g(1,i)-h*(p(1,i)-p(1,i-1))+v(1,1+i)+v(1,i-1)+v(2,i))/3;
end
for j = 2:n-1
    for i = 2:n
        v(j,i) = (h^2*g(j,i)-h*(p(j,i)-p(j,i-1))+v(j+1,i)+v(j-1,i)+v(j,i+1)+v(j,i-1))/4;
    end
end
for i = 2:n
    v(n,i) = (h^2*g(n,i)-h*(p(n,i)-p(n,i-1))+v(n,i-1)+v(n,i+1)+v(n-1,i))/3;
end

%% 计算内部残量的同时更新内部速度和单元压力
% res = zeros(n);
for i = 2:n-1
    for j = 2:n-1
        r = (v(i,j)-v(i,j+1))/h-(u(i+1,j)-u(i,j))/h+d(i,j);
%         res(i,j) = r;
        
        delta = r*h/4;
        v(i,j+1) = v(i,j+1)+delta;
        v(i,j) = v(i,j)-delta;
        u(i+1,j) = u(i+1,j)+delta;
        u(i,j) = u(i,j)-delta;
        
        p(i,j) = p(i,j)+r;
        p(i+1,j) = p(i+1,j)-r/4;
        p(i-1,j) = p(i-1,j)-r/4;
        p(i,j-1) = p(i,j-1)-r/4;
        p(i,j+1) = p(i,j+1)-r/4;
        
    end
end

% for i = 2:n-1
%     for j = 2:n-1
%         delta = res(i,j)*h/4;
%         v(i,j+1) = v(i,j+1)+delta;
%         v(i,j) = v(i,j)-delta;
%         u(i+1,j) = u(i+1,j)+delta;
%         u(i,j) = u(i,j)-delta;
%     end
% end
%% 计算边界残量并更新速度和压力

for i = 2:n-1
    r = (u(1,i)-u(2,i))/h+(v(1,i)-v(1,i+1))/h+d(1,i);
%     r = (u(1,i)-u(2,i))/h+(v(1,i)-v(1,i+1))/h;
%     res(1,i) = r;
    
    delta = r*h/3;
    u(2,i) = u(2,i)+delta;
    v(1,i) = v(1,i)-delta;
    v(1,i+1) = v(1,i+1)+delta;
    
    p(1,i) = p(1,i)+r*4/3;
%     p(n,i) = p(n,i)+r;
    p(1,i-1) = p(1,i-1)-r/3;
    p(1,i+1) = p(1,i+1)-r/3;
    p(2,i) = p(2,i)-r/3;
end

for i = 2:n-1
    r = (u(n,i)-u(n+1,i))/h+(v(n,i)-v(n,i+1))/h+d(n,i);
%     res(n,i) = r;
    delta = r*h/3;
    u(n,i) = u(n,i)-delta;
    v(n,i) = v(n,i)-delta;
    v(n,i+1) = v(n,i+1)+delta;
    
    p(n,i) = p(n,i)+r*4/3;
%     p(n,i) = p(n,i)+r;
    p(n,i-1) = p(n,i-1)-r/3;
    p(n,i+1) = p(n,i+1)-r/3;
    p(n-1,i) = p(n-1,i)-r/3;
end

for i = 2:n-1
    r = (u(i,1)-u(i+1,1))/h+(v(i,1)-v(i,2))/h+d(i,1);
%     res(i,1) = r;
    delta = r*h/3;
    u(i,1) = u(i,1)-delta;
    u(i+1,1) = u(i+1,1)+delta;
    v(i,2) = v(i,2)+delta;
%     
    p(i,1) = p(i,1)+r*4/3;
%     p(i,1) = p(i,1)+r;
    p(i+1,1) = p(i+1,1)-r/3;
    p(i-1,1) = p(i-1,1)-r/3;
    p(i,2) = p(i,2)-r/3;
end

for i = 2:n-1
    r = (u(i,n)-u(i+1,n))/h+(v(i,n)-v(i,n+1))/h+d(i,n);
%     r = (u(i,n)-u(i+1,n))/h+(v(i,n)-v(i,n+1))/h;
%     res(i,n) = r;
    
    delta = r*h/3;
    u(i,n) = u(i,n)-delta;
    u(i+1,n) = u(i+1,n)+delta;
    v(i,n) = v(i,n)-delta;
    
    p(i,n) = p(i,n)+r*4/3;
%     p(i,n) = p(i,n)+r;
    p(i+1,n) = p(i+1,n)-r/3;
    p(i-1,n) = p(i-1,n)-r/3;
    p(i,n-1) = p(i,n-1)-r/3;
end




%% 更新顶点残量、速度和压力

r = (u(1,1)-u(2,1))/h+(v(1,1)-v(1,2))/h+d(1,1);
% r = (u(1,1)-u(2,1))/h+(v(1,1)-v(1,2))/h;
% res(1,1) = r;
delta = r*h/2;
v(1,2) = v(1,2)+delta;
u(2,1) = u(2,1)+delta;

% p(1,1) = p(1,1)+r*2;
p(1,1) = p(1,1)+r;
p(1,2) = p(1,2)-r/2;
p(2,1) = p(2,1)-r/2;


r = (u(n,1)-u(n+1,1))/h+(v(n,1)-v(n,2))/h+d(n,1);
% r = (u(n,1)-u(n+1,1))/h+(v(n,1)-v(n,2))/h;
% res(n,1) = r;
delta = r*h/2;
u(n,1) = u(n,1)-delta;
v(n,2) = v(n,2)+delta;

% p(n,1) = p(n,1)+r*2;
p(n,1) = p(n,1)+r;
p(n,2) = p(n,2)-r/2;
p(n-1,1) = p(n-1,1)-r/2;


r = (v(n,n)-v(n,n+1))/h+(u(n,n)-u(n+1,n))/h+d(n,n);
% r = (v(n,n)-v(n,n+1))/h+(u(n,n)-u(n+1,n))/h;
% res(n,n) = r;
% delta = r*h/2;
v(n,n) = v(n,n)-delta;
u(n,n) = u(n,n)-delta;

% p(n,n) = p(n,n)+r*2;
p(n,n) = p(n,n)+r;
p(n,n-1) = p(n,n-1)-r/2;
p(n-1,n) = p(n-1,n)-r/2;


r = (v(1,n)-v(1,n+1))/h+(u(1,n)-u(2,n))/h+d(1,n);
% r = (v(1,n)-v(1,n+1))/h+(u(1,n)-u(2,n))/h;
% res(1,n) = r;
delta = r*h/2;
u(2,n) = u(2,n)+delta;
v(1,n) = v(1,n)-delta;
% p(1,n) = p(1,n)+r*2;
p(1,n) = p(1,n)+r;
p(1,n-1) = p(1,n-1)-r/2;
p(2,n) = p(2,n)-r/2;
% 
% 
% %% 更新 u,v
% u(3:n,2:n-1) = u(3:n,2:n-1)+res(2:n-1,2:n-1)*h/4;
% u(2:n-1,2:n-1) = u(2:n-1,2:n-1)-res(2:n-1,2:n-1)*h/4;
% v(2:n-1,3:n) = v(2:n-1,3:n)+res(2:n-1,2:n-1)*h/4;
% v(2:n-1,2:n-1) = v(2:n-1,2:n-1)-res(2:n-1,2:n-1)*h/4;
% 
% u(3:n,[1,n]) = u(3:n,[1,n])+res(2:n-1,[1,n])*h/3;
% u(2:n-1,[1,n]) = u(2:n-1,[1,n])-res(2:n-1,[1,n])*h/3;
% u(2,2:n-1) = u(2,2:n-1)+res(1,2:n-1)*h/3;
% u(n,2:n-1) = u(n,2:n-1)-res(n,2:n-1)*h/3;
% v([1,n],3:n) = v([1,n],3:n)+res([1,n],2:n-1)*h/3;
% v([1,n],2:n-1) = v([1,n],2:n-1)-res([1,n],2:n-1)*h/3;
% v(2:n-1,2) = v(2:n-1,2)+res(2:n-1,1)*h/3;
% v(2:n-1,n) = v(2:n-1,n)-res(2:n-1,n)*h/3;
% 
% u(2,[1,n]) = u(2,[1,n])+res(1,[1,n])*h/2;
% u(n,[1,n]) = u(n,[1,n])-res(n,[1,n])*h/2;
% v([1,n],2) = v([1,n],2)+res([1,n],1)*h/2;
% v([1,n],n) = v([1,n],n)-res([1,n],n)*h/2;

% %% 使用迭代法更新残量
% 
% rd = d-(u(2:n+1,1:n)-u(1:n,1:n))/h-(v(1:n,2:n+1)-v(1:n,1:n))/h;
% rd = rd*h^2;
% r = 1; r0 = norm(rd(:));
% 
% delta = zeros(n);
% while r/r0 > 1e-8
% %     cnt = cnt+1;
%     tmp = delta;
%     delta(1,1) = (rd(1,1)+delta(1,2)+delta(2,1))/2;
%     for i = 2:n-1
%         delta(1,i) = (rd(1,i)+delta(1,i-1)+delta(1,i+1)+delta(2,i))/3;
%     end
%     delta(1,n) = (rd(1,n)+delta(1,n-1)+delta(2,n))/2;
%     for i = 2:n-1
%         for j = 2:n-1
%             delta(i,j) = (rd(i,j)+delta(i-1,j)+delta(i+1,j)+delta(i,j-1)+delta(i,j+1))/4;
%         end
%     end
%     delta(n,1) = (rd(n,1)+delta(n,2)+delta(n-1,1))/2;
%     for i = 2:n-1
%         delta(n,i) = (rd(n,i)+delta(n,i-1)+delta(n,i+1)+delta(n-1,i))/3;
%     end
%     delta(n,n) = (rd(n,n)+delta(n,n-1)+delta(n-1,n))/2;
%     r = norm(delta(:)-tmp(:));
% end
% 
% u(2:n,:) = u(2:n,:)+(delta(2:n,:)-delta(1:n-1,:))/h;
% v(:,2:n) = v(:,2:n)+(delta(:,2:n)-delta(:,1:n-1))/h;
% p(2:n-1,2:n-1) = p(2:n-1,2:n-1)-(4*delta(2:n-1,2:n-1)-delta(1:n-2,2:n-1)-delta(2:n-1,1:n-2)-delta(3:n,2:n-1)-delta(2:n-1,3:n))/h^2;
% p(1,1) = p(1,1)-(2*delta(1,1)-delta(1,2)-delta(2,1))/h^2;
% p(1,n) = p(1,1)-(2*delta(1,n)-delta(1,n-1)-delta(2,n))/h^2;
% p(n,1) = p(n,1)-(2*delta(n,1)-delta(n-1,1)-delta(n,2))/h^2;
% p(n,n) = p(n,n)-(2*delta(n,n)-delta(n,n-1)-delta(n-1,n))/h^2;
% p(1,2:n-1) = p(1,2:n-1)-(3*delta(1,2:n-1)-delta(2,2:n-1)-delta(1,1:n-2)-delta(1,3:n))/h^2;
% p(n,2:n-1) = p(n,2:n-1)-(3*delta(n,2:n-1)-delta(n-1,2:n-1)-delta(n,1:n-2)-delta(n,3:n))/h^2;
% p(2:n-1,1) = p(2:n-1,1)-(3*delta(2:n-1,1)-delta(2:n-1,2)-delta(1:n-2,1)-delta(3:n,1))/h^2;
% p(2:n-1,n) = p(2:n-1,n)-(3*delta(2:n-1,n)-delta(2:n-1,n-1)-delta(1:n-2,n)-delta(3:n,n))/h^2;

%% 返回
res = residual_(n,u,v,p,f,g,d);

out.p = p;
out.u = u;
out.v = v;
out.res = res;
