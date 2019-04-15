%BIN_PDF  Probability density function of a Binomial distribution.
%
% Syntax:
%    p = bin_pdf(x,rho,n)
%
% In:
%   x   - Location where to evaluate the PDF.
%   rho - 'Probability' parameter 
%   n   - 'Sample size' parameter

% Out:
%   p   - Density at the given location.
%
% Description:
%   Probability density function for a Binomial distribution.

% History:
%   12.02.2008  Changed the name to bin_pdf.
%   07.02.2008  The first official version.
%
% Copyright (C)  2003 Simo Särkkä
%                2008 Jouni Hartikainen
%
% $Id: bin_pdf.m $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function p = bin_pdf(x,rho,n)
  if (x > n) | (n < 0)
    p = 0;
  else
    p = nchoosek(n,x) * rho^x * (1-rho)^(n-x);
  end
  
