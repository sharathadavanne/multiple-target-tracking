%GMM_RND  RND from Multivariate Mixture of Gaussians
%
% Syntax:
%   X = GMM_RND(M,S,PJ,N)
%
% Author:
%   Simo Särkkä, 2003
%
% In:
%   M  - DxNC   matrix of means.
%   S  - DxDxNC matrix of covariances
%   PJ - Probabilities of mixture components
%   N  - Number of samples (default 1)
%
% Out:
%   X  - DxN matrix of samples
%   
% Description:
%   Draw values of from multivariate mixture
%   of Gaussians distribution
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
% $Id: gmm_pdf.m,v 1.1 2003/06/23 18:09:05 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function X = gmm_rnd(M,S,PJ,N)

  if nargin<4
    N=[];
  end
  if isempty(N)
    N=1;
  end

  SJ = sum(PJ);
  if SJ==0
    PJ = ones(1,length(PJ))/length(PJ);
  else
    PJ = PJ ./ SJ;
  end
  PJ = cumsum(PJ);
  X  = zeros(size(M,1),N);
  for j=1:N
    c = min(find(PJ > rand));
    if isempty(c)
      c = 1;
    end
    X(:,j) = gauss_rnd(M(:,c),S(:,:,c),1);
  end

  
