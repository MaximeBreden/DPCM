function f_coeffs=computeCheb(f,K,disp)

if nargin<3
    disp=0;
end

thetaK=pi*(1:2:2*K-1)'/(2*K);
xK=cos(thetaK);

f_eval=f(xK);
f_coeffs=dct(f_eval);
f_coeffs(1)=f_coeffs(1)/sqrt(K);
f_coeffs(2:end)=f_coeffs(2:end)/sqrt(2*K);

if disp
    figure
    plot(xK,f_eval,'*r')
    hold on
    plotCheb(f_coeffs,-1,1)
end