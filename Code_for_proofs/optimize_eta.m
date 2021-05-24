function eta=optimize_eta(coeffs_Y,coeffs_Z1,coeffs_Z2,eta0,criterion)

% Here we numerically choose a set of weigths eta, according to criterion.
% We can try to either make Z1 as small as possible, or make the
% discriminant of the radii polynomial P as large as possible, or make
% r_max-r_min as large as possible.

% eta should be a row
% Y=eta*coeffs_Y;
% Z1=max((eta*coeffs_Z1)./eta); 
% Z2=max(max(sum(permute(bsxfun(@rdivide,permute(bsxfun(@rdivide,bsxfun(@times,coeffs_Z2,eta'),eta),[1,3,2]),eta),[1,3,2]),1)));

if nargin<5
    criterion='Delta';
end

if isempty(eta0)
      eta0=ones(1,10);
    noguess=1;
else
    noguess=0;
end

switch criterion
    case 'Z1'
        disp('Optimizing eta to make Z1 as small as possible')
        [eta,~]=eigs(coeffs_Z1',1);
        eta=abs(eta)';
        eta=max(eta,10^-5);
        disp(['eta = ',num2str(eta)])
        return
    case 'Delta'
        disp('Optimizing eta to make Delta as large as possible')
        g=@(eta) 2*(eta*coeffs_Y)*(max(max(sum(permute(bsxfun(@rdivide,permute(bsxfun(@rdivide,bsxfun(@times,coeffs_Z2,eta'),eta),[1,3,2]),eta),[1,3,2]),1))))...
                 -(1-(max((eta*coeffs_Z1)./eta)))^2;
    case 'Largest_interval'
        disp('Optimizing eta to make the validation interval as large as possible')
        g=@(eta) 2*(eta*coeffs_Y)/(max(max(sum(permute(bsxfun(@rdivide,permute(bsxfun(@rdivide,bsxfun(@times,coeffs_Z2,eta'),eta),[1,3,2]),eta),[1,3,2]),1))))...
                 -(1-(max((eta*coeffs_Z1)./eta)))^2/(max(max(sum(permute(bsxfun(@rdivide,permute(bsxfun(@rdivide,bsxfun(@times,coeffs_Z2,eta'),eta),[1,3,2]),eta),[1,3,2]),1))))^2;
    otherwise
        warning('Invalid criterion for the optimization of eta')
        eta=eta0;
        disp(['eta = ',num2str(eta)])
        return
end

eta=patternsearch(g,eta0,[],[],[],[],10^-2*ones(1,10),ones(1,10),@constraint);
%eta=ga(g,6,[],[],[],[],10^-2*ones(1,10),ones(1,10),@constraint);
if noguess
    eta1=patternsearch(g,eta+10^-2*(2*rand(1,10)-1),[],[],[],[],10^-2*ones(1,10),ones(1,10),@constraint);
    if g(eta1)<g(eta)
        eta=eta1;
    end
end
disp(['eta = ',num2str(eta)])


function [c,ceq] = constraint(eta)
ceq = [];
c = [max((eta*coeffs_Z1)./eta)-0.95;
     2*(eta*coeffs_Y)*(max(max(sum(permute(bsxfun(@rdivide,permute(bsxfun(@rdivide,bsxfun(@times,coeffs_Z2,eta'),eta),[1,3,2]),eta),[1,3,2]),1))))-(1-(max((eta*coeffs_Z1)./eta)))^2];
%the nonlinear constraint imposes that Z1<=0.95 and Delta = (1-Z1)^2-2Y*Z2>=0
end

end