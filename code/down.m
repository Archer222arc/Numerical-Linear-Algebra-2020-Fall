function out = down(n,rf,rg,rd)
%把n限制到n/2的网格上
if nargin < 4
    disp("错误输入，请输入n,u,v,p")
    return 
end

%% 初始化
dg = zeros(n/2,n/2+1);
df = zeros(n/2+1,n/2);
m = n/2;

%% 限制
dg(1:m,2:m) = rg(1:2:2*m-1,2:2:n-2)/8+rg(1:2:2*m-1,4:2:n)/8+rg(2:2:2*m,2:2:n-2)/8+...
    rg(2:2:2*m,4:2:n)/8+rg(1:2:2*m-1,3:2:n-1)/4+rg(2:2:2*m,3:2:n-1)/4;
df(2:m,1:m) = rf(2:2:n-1,1:2:n-1)/8+rf(4:2:n,1:2:n-1)/8+rf(2:2:n-2,2:2:n)/8+...
    rf(4:2:n,2:2:n)/8+rf(3:2:n-1,1:2:n-1)/4+rf(3:2:n-1,2:2:n)/4;
dd = (rd(1:2:n,1:2:n)+rd(2:2:n,2:2:n)+rd(1:2:n,2:2:n)+rd(2:2:n,1:2:n))/4;

%% 输出
out.rf = df;
out.rg = dg;
out.rd = dd;