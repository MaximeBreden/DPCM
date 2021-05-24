function [coeffs_Y,coeffs_Z1,coeffs_Z2,result,A]=coeffs_for_proof(X,para,nu,rstar,eta,A)

%Here we compute the various bounds described in Section 4, block-wise,
%i.e. for each component Psi, E, C, etc, but without putting everything
%together yet, because we still have to chose the weight eta.

if strcmp(para.case,'Poten')
    K=(length(X)-5)/5;
    V=para.V;
else
    error('Invalid value of para.case')
end

Adag=DF(X,para);
if isempty(A)
    if exist('intval','file') && isintval(X(1))
        A=intval(inv(Adag.mid));
    else
        A=inv(Adag);
    end
end
result=1;

%% Z1          
coeffs_Z1=max(coeffs_Z1_finite(X,para,nu,A),coeffs_Z1_tail(X,para,nu,A));

if isempty(eta)
    if exist('intval','file') && isintval(X(1))
        [eta,~]=eigs(coeffs_Z1.mid',1);
        eta=intval(abs(eta)');
    else
        [eta,~]=eigs(coeffs_Z1',1);
        eta=abs(eta)';
    end
end

Z1_opt=max((eta*coeffs_Z1)./eta);
if Z1_opt>=1
    disp('The proof will fail no matter the choice of weights, Z1_opt is too large')
    disp(Z1_opt)
    coeffs_Y=NaN;
    coeffs_Z2=NaN;
    result=0;
    return
end


%% Z0
temp=eye(5*K+5)-A*Adag;
coeffs_Z0=norm_op_blocks_int(temp,nu); % BC?
coeffs_Z1=coeffs_Z1+coeffs_Z0;

%% Y
K_ext=2*K;
ind_red=[1:K,K_ext+(1:K),2*K_ext+(1:K),3*K_ext+(1:K),4*K_ext+(1:K),5*K_ext+(1:5)];

X_ext=zeros(5*K_ext+5,1);
if exist('intval','file') && isintval(X(1))
    X_ext=intval(X_ext);
end
X_ext(ind_red)=X;

FX_ext=F(X_ext,para);
FX_red=FX_ext(ind_red);

AFX_red=A*FX_red;
AFX_ext=FX_ext;
AFX_ext(ind_red)=AFX_red;

coeffs_Y=zeros(10,1);
if exist('intval','file') && isintval(nu)
    coeffs_Y=intval(coeffs_Y);
end
weights_ext=[1 2*nu.^(1:K_ext-1)];
coeffs_Y(1)=weights_ext*abs(AFX_ext(1:K_ext));
coeffs_Y(2)=weights_ext*abs(AFX_ext(K_ext+1:2*K_ext));
coeffs_Y(3)=weights_ext*abs(AFX_ext(2*K_ext+1:3*K_ext));
coeffs_Y(4)=weights_ext*abs(AFX_ext(3*K_ext+1:4*K_ext));
coeffs_Y(5)=weights_ext*abs(AFX_ext(4*K_ext+1:5*K_ext));
coeffs_Y(6:10)=abs(AFX_ext(5*K_ext+(1:5)));

%% Z2 

Psi=X(1:K);
E=X(K+1:2*K);
C=X(2*K+1:3*K);
N=X(3*K+1:4*K);
P=X(4*K+1:5*K);
JC=X(5*K+1);
JN=X(5*K+2);
JP=X(5*K+3);
delta=X(5*K+4);
l=X(5*K+5);

Cstar=rstar/eta(3);
Nstar=rstar/eta(4);
Pstar=rstar/eta(5);
if exist('intval','file') && isintval(nu)
    JC=JC+infsup(-rstar,rstar)/eta(6);
    delta=delta+infsup(-rstar,rstar)/eta(9);
    l=l+infsup(-rstar,rstar)/eta(10);
    Psistar_BC=infsup(-rstar,rstar)/eta(1);
    Cstar_BC=infsup(-rstar,rstar)/eta(3);
    Nstar_BC=infsup(-rstar,rstar)/eta(4);
    Pstar_BC=infsup(-rstar,rstar)/eta(5);
else
    Psistar_BC=0;
    Cstar_BC=0;
    Nstar_BC=0;
    Pstar_BC=0;
end

BC0=[1 2*(-1).^(1:K-1)];
BC1=[1 2*ones(1,K-1)];
if exist('intval','file') && isintval(nu)
    Int=diag(1./intval(2*(1:K)))*spdiags([-ones(K+1,1),ones(K+1,1)],[0,2],K,K+1);
    e0=intval([1;zeros(K,1)]);
else
    Int=diag(1./(2*(1:K)))*spdiags([-ones(K+1,1),ones(K+1,1)],[0,2],K,K+1);
    e0=[1;zeros(K,1)];
end
C0=[C;0];
N0=[N;0];
P0=[P;0];

coeffs_A_int=norm_op_blocks_int(A,nu);
coeffs_A_BC=norm_op_blocks_BC(A,nu);
%clear A
for i=1:5
    coeffs_A_int(i,i)=max(coeffs_A_int(i,i),1);
end

D2F_int=zeros(10,10,10);
D2F_BC=zeros(10,10,10);
if exist('intval','file') && isintval(nu)
    D2F_int=intval(D2F_int);
    D2F_BC=intval(D2F_BC);
end    
    
%FPsi
D2F_BC(1,2,10)=abs(2*para.alpha0/l^2); %(l,E')
D2F_BC(1,10,2)=abs(2*para.alpha0/l^2); %(l',E)
D2F_BC(1,10,10)=abs(4*para.alpha0/l^3*BC0*E); %(l,l')

%FE
D2F_BC(2,2,10)=abs(2*para.alpha1/l^2); %(l,E')
D2F_BC(2,10,2)=abs(2*para.alpha1/l^2); %(l',E)
D2F_BC(2,10,10)=abs(4*para.alpha1/l^3*BC0*E); %(l,l')

D2F_int(2,3,10)=abs(para.coupling*l/(2*para.lambda^2)*2); %(l,C')
D2F_int(2,10,3)=abs(para.coupling*l/(2*para.lambda^2)*2); %(l',C)
D2F_int(2,4,10)=abs(para.coupling*l/(2*para.lambda^2)*1); %(l,N')
D2F_int(2,10,4)=abs(para.coupling*l/(2*para.lambda^2)*1); %(l',N)
D2F_int(2,5,10)=abs(para.coupling*l/(2*para.lambda^2)*3); %(l,P')
D2F_int(2,10,5)=abs(para.coupling*l/(2*para.lambda^2)*3); %(l',P)
D2F_int(2,10,10)=para.coupling/(2*para.lambda^2)*(2*nu.^(1:K)*abs(Int*(3*P0-N0+2*C0+para.rho*e0))+3*Pstar+Nstar+2*Cstar); %(l,l')

%FC
D2F_BC(3,1,3)=abs(2*l*para.der21rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(Psi,C')
D2F_BC(3,3,1)=abs(2*l*para.der12rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(Psi',C)
D2F_BC(3,1,1)=abs(2*l*para.der22rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(Psi,Psi')
D2F_BC(3,3,3)=abs(2*l*para.der11rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(C,C')
D2F_BC(3,1,10)=abs(2*para.der2rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(Psi,l')
D2F_BC(3,10,1)=abs(2*para.der2rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(l',C)
D2F_BC(3,3,10)=abs(2*para.der1rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(C,l')
D2F_BC(3,10,3)=abs(2*para.der1rC0(BC0*[C Psi]+[Cstar_BC Psistar_BC])); %(C',l)

D2F_int(3,2,3)=abs(para.zC); %(E,C')
D2F_int(3,3,2)=abs(para.zC); %(E',C)
D2F_int(3,3,9)=abs(para.epsC/2*l); %(C,delta')
D2F_int(3,9,3)=abs(para.epsC/2*l); %(C',delta)
D2F_int(3,3,10)=abs(para.epsC/2*delta); %(C,l')
D2F_int(3,10,3)=abs(para.epsC/2*delta); %(C',l)
D2F_int(3,9,10)=2*nu.^(1:K)*abs(Int*(para.epsC/2*C0))+para.epsC/2*Cstar; %(delta,l')
D2F_int(3,10,9)=2*nu.^(1:K)*abs(Int*(para.epsC/2*C0))+para.epsC/2*Cstar; %(delta',l)

%FN
D2F_BC(4,1,4)=abs(2*l*para.der21rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(Psi,N')
D2F_BC(4,4,1)=abs(2*l*para.der12rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(Psi',N)
D2F_BC(4,1,1)=abs(2*l*para.der22rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(Psi,Psi')
D2F_BC(4,4,4)=abs(2*l*para.der11rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(N,N')
D2F_BC(4,1,10)=abs(2*para.der2rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(Psi,l')
D2F_BC(4,10,1)=abs(2*para.der2rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(l',N)
D2F_BC(4,4,10)=abs(2*para.der1rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(N,l')
D2F_BC(4,10,4)=abs(2*para.der1rN0(BC0*[N Psi]+[Nstar_BC Psistar_BC])); %(N',l)

D2F_int(4,2,4)=abs(para.zN); %(E,N')
D2F_int(4,4,2)=abs(para.zN); %(E',N)
D2F_int(4,4,9)=abs(para.epsN/2*l); %(N,delta')
D2F_int(4,9,4)=abs(para.epsN/2*l); %(N',delta)
D2F_int(4,4,10)=abs(para.epsN/2*delta); %(N,l')
D2F_int(4,10,4)=abs(para.epsN/2*delta); %(N',l)
D2F_int(4,9,10)=2*nu.^(1:K)*abs(Int*(para.epsN/2*N0))+para.epsN/2*Nstar; %(delta,l')
D2F_int(4,10,9)=2*nu.^(1:K)*abs(Int*(para.epsN/2*N0))+para.epsN/2*Nstar; %(delta',l)

%FP
D2F_BC(5,1,5)=abs(2*l*para.der21rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(Psi,P')
D2F_BC(5,5,1)=abs(2*l*para.der12rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(Psi',P)
D2F_BC(5,1,1)=abs(2*l*para.der22rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(Psi,Psi')
D2F_BC(5,5,5)=abs(2*l*para.der11rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(P,P')
D2F_BC(5,1,10)=abs(2*para.der2rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(Psi,l')
D2F_BC(5,10,1)=abs(2*para.der2rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(l',P)
D2F_BC(5,5,10)=abs(2*para.der1rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(P,l')
D2F_BC(5,10,5)=abs(2*para.der1rP0(BC0*[P Psi]+[Pstar_BC Psistar_BC])); %(P',l)

D2F_int(5,2,5)=abs(para.zP); %(E,P')
D2F_int(5,5,2)=abs(para.zP); %(E',P)
D2F_int(5,5,9)=abs(para.epsP/2*l); %(P,delta')
D2F_int(5,9,5)=abs(para.epsP/2*l); %(P',delta)
D2F_int(5,5,10)=abs(para.epsP/2*delta); %(P,l')
D2F_int(5,10,5)=abs(para.epsP/2*delta); %(P',l)
D2F_int(5,9,10)=2*nu.^(1:K)*abs(Int*(para.epsP/2*P0))+para.epsP/2*Pstar; %(delta,l')
D2F_int(5,10,9)=2*nu.^(1:K)*abs(Int*(para.epsP/2*P0))+para.epsP/2*Pstar; %(delta',l)

%FJC
D2F_BC(6,1,3)=abs(2*l*para.der21rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(Psi,C')
D2F_BC(6,3,1)=abs(2*l*para.der12rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(Psi',C)
D2F_BC(6,1,1)=abs(2*l*para.der22rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(Psi,Psi')
D2F_BC(6,3,3)=abs(2*l*para.der11rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(C,C')
D2F_BC(6,1,10)=abs(2*para.der2rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(Psi,l')
D2F_BC(6,10,1)=abs(2*para.der2rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(l',C)
D2F_BC(6,3,10)=abs(2*para.der1rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(C,l')
D2F_BC(6,10,3)=abs(2*para.der1rC1([BC1*[C Psi]+[Cstar_BC Psistar_BC] V])); %(C',l)

%FJN
D2F_BC(7,1,4)=abs(2*l*para.der21rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(Psi,N')
D2F_BC(7,4,1)=abs(2*l*para.der12rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(Psi',N)
D2F_BC(7,1,1)=abs(2*l*para.der22rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(Psi,Psi')
D2F_BC(7,4,4)=abs(2*l*para.der11rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(C,C')
D2F_BC(7,1,10)=abs(2*para.der2rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(Psi,l')
D2F_BC(7,10,1)=abs(2*para.der2rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(l',N)
D2F_BC(7,4,10)=abs(2*para.der1rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(N,l')
D2F_BC(7,10,4)=abs(2*para.der1rN1([BC1*[N Psi]+[Nstar_BC Psistar_BC] V])); %(N',l)

%FJP
D2F_BC(8,1,5)=abs(2*l*para.der21rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(Psi,P')
D2F_BC(8,5,1)=abs(2*l*para.der12rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(Psi',P)
D2F_BC(8,1,1)=abs(2*l*para.der22rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(Psi,Psi')
D2F_BC(8,5,5)=abs(2*l*para.der11rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(C,C')
D2F_BC(8,1,10)=abs(2*para.der2rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(Psi,l')
D2F_BC(8,10,1)=abs(2*para.der2rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(l',P)
D2F_BC(8,5,10)=abs(2*para.der1rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(P,l')
D2F_BC(8,10,5)=abs(2*para.der1rP1([BC1*[P Psi]+[Pstar_BC Psistar_BC] V])); %(P',l)

%Fdelta
D2F_BC(9,1,1)=abs(para.der2chi(BC0*Psi+Psistar_BC)); %(Psi,Psi')

%Fl
D2F_BC(10,6,9)=abs(para.kappa/(2*delta^2*para.epsC)); %(JC,delta')
D2F_BC(10,9,6)=abs(para.kappa/(2*delta^2*para.epsC)); %(JC',delta)
D2F_BC(10,10,10)=abs(para.kappa*JC/(delta^3*para.epsC)); %(delta,delta')

coeffs_Z2=zeros(10,10,10);
if exist('intval','file') && isintval(nu)
    coeffs_Z2=intval(coeffs_Z2);
end
for i1=1:10
    for i2=1:10
        coeffs_Z2(i1,:,:)=coeffs_Z2(i1,:,:)+coeffs_A_int(i1,i2)*D2F_int(i2,:,:)+coeffs_A_BC(i1,i2)*D2F_BC(i2,:,:);
    end
end


%% Auxiliary functions
function Mat_coeffs=norm_op_blocks_int(B,nu)
%!! Does not take into account the boundary terms !!
weights=[1 2*nu.^(1:K-1)];
Mat_coeffs=zeros(10,10);
if exist('intval','file') && isintval(nu)
    Mat_coeffs=intval(Mat_coeffs);
end
for m=1:5
    for n=1:5
        Mat_coeffs(m,n)=max((weights*abs(B(K*(m-1)+(1:K),K*(n-1)+(2:K))))./weights(2:K));
    end
end
for m=6:10
    for n=1:5
        Mat_coeffs(m,n)=max(abs(B(5*K+m-5,K*(n-1)+(2:K)))./weights(2:K));
    end
end
end

function Mat_coeffs=norm_op_blocks_BC(B,nu)
weights=[1 2*nu.^(1:K-1)];
Mat_coeffs=zeros(10,10);
if exist('intval','file') && isintval(nu)
    Mat_coeffs=intval(Mat_coeffs);
end
for m=1:5
    for n=1:5
        Mat_coeffs(m,n)=weights*abs(B(K*(m-1)+(1:K),K*(n-1)+1));
    end
    for n=6:10
        Mat_coeffs(m,n)=weights*abs(B(K*(m-1)+(1:K),5*K+n-5));
    end
end
for m=6:10
    for n=1:5
        Mat_coeffs(m,n)=abs(B(5*K+m-5,K*(n-1)+1));
    end
    for n=6:10
        Mat_coeffs(m,n)=abs(B(5*K+m-5,5*K+n-5));
    end
end
end

end
