%GAMMA_PDF  Probability density function of a Gamma distribution
%
% Syntax:
%   p = gamma_pdf(x,alpha,beta,mu)
%
% In:
%   x     - Locations where to evaluate the PDF 
%   alpha - Parameter of the distribution
%   beta  - Parameter of the distribution
%   mu    - Mean of the distribution
%
%
% Out:
%   p     - Density at the given locations
%
% Description:
%   Probability density function of a Gamma distribution
%
% See also:
%   GAMMA_CDF

% History:
%    11.02.2008 First Official version
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: gamma_pdf.m, $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.


function p = gamma_pdf(x,gam,beta,mu)

  if nargin < 3
    beta = [];
  end
  if nargin < 4
    mu = [];
  end

  if isempty(beta)
    beta = 1;
  end
  if isempty(mu)
    mu = 0;
  end
  
  %
  % Convert to standard form
  %
  x = (x-mu)/beta;

  %
  % Compute probability
  %
  p = x.^(gam-1).*exp(-x)./gamma(gam);
  p = p / beta;

