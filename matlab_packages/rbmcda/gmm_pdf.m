%GMM_PDF  Multivariate Mixture of Gaussians PDF
%
% Syntax:
%   [P,PJX] = GMM_PDF(X,M,S,PJ)
%
% In:
%   X  - Dx1 value or N values as DxN matrix
%   M  - DxNC   matrix of means.
%   S  - DxDxNC matrix of covariances
%   PJ - Probabilities of mixture components
%
% Out:
%   P  - Probability of each PX
%   PJX - Probabilities P(j|x) for each x[n]
%   
% Description:
%   Calculate values of PDF (Probability Density
%   Function) of multivariate mixture Gaussian
%   distribution
%
%     p(X) = sum p(j) N(X | M_j, S_j)
%             j
%
%   Function returns probability of X in PDF. If multiple
%   X's are given (as multiple columns), function returns
%   probabilities for each of them.
%
% See also:
%   GAUSS_PDF

% History:
%   23.06.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%
% $Id: gmm_pdf.m,v 1.2 2003/07/01 10:43:15 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [P,PJX] = gmm_pdf(X,M,S,PJ)

  P = zeros(1,size(X,2));
  SJ = sum(PJ);
  if SJ==0
    PJ = ones(1,length(PJ))/length(PJ);
  else
    PJ = PJ ./ SJ;
  end
  PJX = zeros(length(PJ),size(X,2));
  for j=1:length(PJ)
    PJX(j,:) = PJ(j)*gauss_pdf(X,M(:,j),S(:,:,j));
    P = P + PJX(j,:);
  end
  PJX = PJX ./ repmat(sum(PJX,1),length(PJ),1);

