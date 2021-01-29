function out = residual_(n,u,v,p,f,g,d)
%返回残量
rf = zeros(n+1,n);
rg = zeros(n,n+1);
h = 1/n;
for i = 2:n
    rf(i,1) = f(i,1)-(p(i,1)-p(i-1,1))/h+(u(i+1,1)+u(i-1,1)+u(i,2)-3*u(i,1))/h^2;
end
for j = 2:n-1
    for i = 2:n
        rf(i,j) = f(i,j)-(p(i,j)-p(i-1,j))/h+(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-4*u(i,j))/h^2;
    end
end
for i = 2:n
    rf(i,n) = f(i,n)-(p(i,n)-p(i-1,n))/h+(u(i+1,n)+u(i-1,n)+u(i,n-1)-3*u(i,n))/h^2;
end


for i = 2:n
    rg(1,i) = g(1,i)-(p(1,i)-p(1,i-1))/h+(v(1,1+i)+v(1,i-1)+v(2,i)-3*v(1,i))/h^2;
end
for j = 2:n-1
    for i = 2:n
        rg(j,i) = g(j,i)-(p(j,i)-p(j,i-1))/h+(v(j+1,i)+v(j-1,i)+v(j,i+1)+v(j,i-1)-4*v(j,i))/h^2;
    end
end
for i = 2:n
    rg(n,i) = g(n,i)-(p(n,i)-p(n,i-1))/h+(v(n,i-1)+v(n,i+1)+v(n-1,i)-3*v(n,i))/h^2;
end

%% 散度残量
rd = -(u(2:n+1,:)-u(1:n,:)+v(:,2:n+1)-v(:,1:n))/h+d;

out.rf = rf;
out.rg = rg;
out.rd = rd;