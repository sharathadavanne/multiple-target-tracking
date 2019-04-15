%GAUSS_RND  Multivariate Gaussian random variables
%
% Syntax:
%   P = GAUSS_RND(M,S,N)
%
% Author:
%   Simo Särkkä, 2002
%
% In:
%   M - Dx1 mean of distibution or N values as DxN matrix.
%   S - DxD covariance matrix
%   N - Number of samples (optional, default 1)
%
% Out:
%   X - DxN matrix of samples.
%   
% Description:
%   Draw N samples from multivariate Gaussian distribution
% 
%     X ~ N(M,S)
%
% See also:
%   GAUSS_PDF

% History:
%   20.11.2002  The first official version.
%
% Copyright (C) 2002 Simo Särkkä
%
% $Id: gauss_rnd.m,v 1.1.1.1 2003/09/15 10:54:34 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.


function X = gauss_rnd(M,S,N)

  if nargin < 3
    N = 1;
  end
  
  L = chol(S)';
  X = repmat(M,1,N) + L*randn(size(M,1),N);
  