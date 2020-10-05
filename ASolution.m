function [sigma] = ASolution(Vecs, victimSpectra, meanSpec_targetClass)
%Analytical Solution: Returns the analytical perturbation for an adversarial
%attack
%   Ziyi Cao wrote the code and confirmed the math
%   Garth did the math behind this

% updated by Youlin on 2019-12-19: to aloow for multile number of classes by matrix operation 
sigma = Vecs*Vecs'*(victimSpectra-meanSpec_targetClass) ./ (sum(Vecs .^2,2) + 1);

%%%%%%%%%%%%%%%%%%%%%%older scripts%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To Do: Check that Vecs is the right dimensions and then transpose if
%neccessary.
% sizecheck = size(Vecs);
% if sizecheck(2) > sizecheck(1)
%     Vecs = Vecs';
%     sizecheck = size(Vecs);
% end

% if sizecheck(2) == 2
% %analytical solution for perturbation
% sigma=zeros(1340,1);
% denominator=Vecs*Vecs'*(victimSpectra-meanSpec_targetClass);
% sigma = denominator ./ (Vecs(:,1).^2 + Vecs(:,2).^2 + 1);%use matrix operator instead of for loop can increase speed
% end

% if sizecheck(2) == 3
% %analytical solution for perturbation
% sigma=zeros(1340,1);
% denominator=Vecs*Vecs'*(victimSpectra-meanSpec_targetClass);
% sigma = denominator ./ (Vecs(:,1).^2 + Vecs(:,2).^2 + Vecs(:,3).^2 + 1);%use matrix operator instead of for loop can increase speed
% end

end

