%EFF_PARTICLES  Estimate the number of effective particles
%
% Syntax:
%   N_eff = EFF_WEIGHTS(W)
%
% Author:
%   Simo Särkkä, 2002
%
% In:
%   W - Importance weights as a numeric array or a cell array
%       containing particle structures.
%
% Out:
%   N_eff - Estimated number of effective weights
%
% Description:
%   Estimate the effective number of particles from a given
%   set of importance weights.
%
% See also:
%   RESAMPLE

% History:
%   23.01.2008 JH  Added support for particle structures
%   18.02.2003 SS  The first official version
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: eff_weights.m,v 1.1.1.1 2003/09/15 10:54:33 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function N = eff_particles(W)

  % If the weights are given in a vector 
  if isnumeric(W)      
      W = W ./ sum(W);
      N = 1 / sum(W.^2);
  % If the particle structures are given
  elseif iscell(W) & isfield([W{:}],'W')
      tmp = [W{:}];
      WW = [tmp.W];
      WW = WW ./ sum(WW);
      N = 1 / sum(WW.^2);
  else
      error('Weights are not specified correctly');
  end
