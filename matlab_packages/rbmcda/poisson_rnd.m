%POISSON_RND  Draws random samples from a Poisson distribution.
%
% Syntax:
%   X = poisson_rnd(lambda,N)
%
% In:
%   lambda - Rate parameter of the distribution
%   N      - Number of samples (default 1)
%
% Out:
%   X      - 1xN vector of samples
%   
% Description:
%   Draw values of from a Poisson distribution.
%
% See also:
%   GAMMA_RND, GAMMA_PDF

% History:
%   13.02.2008  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: poisson_rnd.m, $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function x = poisson_rnd(lambda,N)
    if nargin < 2
        N = 1;
    end
    
    x = zeros(1,N);
    for i = 1:N
        s = - log (1 - rand) ./ lambda;
        while (s < 1)
            s = (s - log (1 - rand) ./ lambda);
            x(i) = x(i) + 1;
        end
    end