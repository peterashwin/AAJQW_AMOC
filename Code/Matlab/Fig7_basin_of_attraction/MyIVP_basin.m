function [tend,xend] = MyIVP_basin(f,x0,tspan,h,Xu,Xl)
% Solving initial-value problems for ODEs of form x*(t)=f(t,x(t))
%   1D time input and n-D x input to n-D x output
%   Trapezoidal 2nd order

N = round((tspan(2)-tspan(1))/h);
M = size(x0,1);
m = size(x0,2);
t = NaN(N+1,1);
t(1) = tspan(1);

if m > 1
    xt = NaN(M,m,N+1);
    xt(:,:,1) = x0;
else
    xt(:,1) = x0;
end
n=1;
dis=1;
while and(n<=N,dis>1e-3)
    t0 = t(n);
    if m > 1
        x0 = xt(:,:,n);
    else
        x0 = xt(:,n);
    end
    k = f(t0,x0);
    x_temp = x0 + h.*k;
    t0 = tspan(1)+ n*h;
    t(n+1) = t0;
    if m > 1
        xt(:,:,n+1) = x0 + h./2.*(k+f(t0,x_temp));
    else
        xt(:,n+1) = x0 + h./2.*(k+f(t0,x_temp));
    end
    if m > 1
        xend = xt(:,:,n+1);
    else
        xend = xt(:,n+1);
    end
    dis= min(norm(xend-Xu),norm(xend-Xl));
    n=n+1;
    tend=t(n);
end

end