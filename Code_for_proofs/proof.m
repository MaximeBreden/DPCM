function [result,rmin,rmax,eta,A]=proof(X,para,nu,rstar,criterion,eta0,eta,A)

%Here we try to validate a posteriori the approximate solution X, using the
%strategy described in Section 4 of the paper.

rmin=NaN;
rmax=NaN;

if nargin<6
    eta0=[];
end

if nargin<7
    eta=[];
end

if nargin<8
    A=[];
end

[coeffs_Y,coeffs_Z1,coeffs_Z2,result,A]=coeffs_for_proof(X,para,nu,rstar,eta,A);
if not(result)
    eta=NaN*ones(1,10);
    return
end

if nargin<5
    criterion='NoOpt';
end

if isempty(eta) %Computation of an approximately "optimal" set of weights eta
    if exist('intval','file') && isintval(nu)
        eta=intval(optimize_eta(coeffs_Y.sup,coeffs_Z1.sup,coeffs_Z2.sup,eta0,criterion));
    else
        eta=optimize_eta(coeffs_Y,coeffs_Z1,coeffs_Z2,eta0,criterion);
    end
    if min(eta)<=0
        error('Something went wrong during the optimization process')
    end
end

Y=eta*coeffs_Y;
Z1=max((eta*coeffs_Z1)./eta); 
if Z1>=1
    disp('The proof failed, Z1 is too large')
    disp(Z1)
    result=0;
    return
end
for i=1:10
    coeffs_Z2(i,:,:)=coeffs_Z2(i,:,:)*eta(i);
end
for j1=1:10
    coeffs_Z2(:,j1,:)=coeffs_Z2(:,j1,:)/eta(j1);
end
for j2=1:10
    coeffs_Z2(:,:,j2)=coeffs_Z2(:,:,j2)/eta(j2);
end
Z2=max(max(sum(coeffs_Z2,1)));

if exist('intval','file') && isintval(nu)
    disp(['Y = ',num2str(Y.sup)])
    disp(['Z1 = ',num2str(Z1.sup)])
    disp(['Z2 = ',num2str(Z2.sup)])
else
    disp(['Y = ',num2str(Y)])
    disp(['Z1 = ',num2str(Z1)])
    disp(['Z2 = ',num2str(Z2)])
end

Delta=(1-Z1)^2-2*Y*Z2;
if not(Delta>0)
    result=0;
    disp('Validation failed :(')
    if exist('intval','file') && isintval(nu)
        disp(['no positive roots, Delta = ',num2str(Delta.inf)])
    else
        disp(['no positive roots, Delta = ',num2str(Delta)])
    end
    return
end


rmin=((1-Z1)-sqrt(Delta))/Z2;
rmax=min((1-Z1)/Z2,rstar);
if rmax<=rmin
    result=0;
    disp('rstar too small')
    return
end
result=1;
if exist('intval','file') && isintval(nu)
    rmin=rmin.sup;
    rmax=rmax.inf;
end
disp(['[r_min,r_max] = ',num2str(rmin),',',num2str(rmax)])
fprintf('\n')
 