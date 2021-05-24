function Y=change_size(X,K)

%Change the number of Chebyshev modes used to K, padding with 0 if K needs
%to be increased, and truncating if K needs to be decreased.

K_sol=(length(X)-5)/5;

L=min(K_sol,K);
if K_sol<=K
    ind_X=1:5*K_sol+5;
else
    ind_X=[1:L,K_sol+(1:L),2*K_sol+(1:L),3*K_sol+(1:L),4*K_sol+(1:L),5*K_sol+(1:5)];
end

ind_Y=[1:L,K+(1:L),2*K+(1:L),3*K+(1:L),4*K+(1:L),5*K+(1:5)];

Y=zeros(5*K+5,1);
Y(ind_Y)=X(ind_X);