function plotCheb(f,a,b,nb,col)
%Plots the function represented by the Chebyshev coefficients in f, on 
%[a,b] discretized with step size nb.
%Normalization: f(t)=f_0 + 2\sum_{k\geq 2} f_k T_k(2*x/(b-a)+(b+a)/(b-a)).
if nargin<5
    col='b';
end
if nargin<4
    nb=500;
end
if nargin<2
    a=-1;
    b=1;
end

k=0:length(f)-1;
t=(a:(b-a)/nb:b)';
M=cos(acos(2*t/(b-a)-(b+a)/(b-a))*k);
f(2:end)=2*f(2:end);
ft=M*f;

plot(t,ft,col,'Linewidth',2)

