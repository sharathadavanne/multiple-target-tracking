%NMCDA_INIT  NMCDA Data structure initialization
%
% Syntax:
%   S = nmcda_init(np,M0,P0,dt,[W])
%
% Author:
%   Simo Särkkä, 2003
%
% In:
%   np - Number of Monte Carlo samples
%   M0 - Target prior mean
%   P0 - Target prior covariance
%   dt - Time between measurements
%   W  - Prior importance weights (optional, default: uniform)
% 
% Out:
%   S - Struct array of particles
%
% Description:
%   Initialize data structure for Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%
%   Each element of array S represents one particle, which is a
%   data structure containing the following fields:
%
%      M  : Cell array of size 1xT of T target means
%      P  : Cell array of size 1xT of T target covariances
%      W  : Importance weight
%      M0 : Prior mean for targets
%      P0 : Prior covariance for targets
%      T  : Times of last associated measurements
%      B  : Birth indicator
%      D  : Death index, zero means none
%      dt : Elapsed time from last measurement
%
% See also:
%   KF_NMCDA_PREDICT, KF_NMCDA_UPDATE, EKF_NMCDA_PREDICT, EKF_NMCDA_UPDATE,
%   UKF_NMCDA_PREDICT, UKF_NMCDA_UPDATE

% History:
%   29.01.2008  Changed the name to nmcda_init as the same
%               can be used for KF, EKF and UKF based RBMCDA.
%   26.01.2008  Moved weights inside the particle
%               structures.
%   09.12.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: nmcda_init.m, $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function S = nmcda_init(np,M0,P0,dt,W)
  if nargin < 5
      W = ones(1,np)/np;
  end
  str = struct(...
      'M',{{}},... % Cell array 1xT of T target means
      'P',{{}},... % Cell array 1xT of T target covariances
      'W',0,...    % Importance weight
      'M0',M0,...  % Prior mean for targets
      'P0',P0,...  % Prior covariance for targets
      'T',[],...   % Times of last associated measurements
      'B',0,...    % Birth indicator
      'D',0,...    % Death index, zero means none
      'dt',dt ...  % Elapsed time from last measurement
  );
  
  S = cell(1,np);
  
  for i=1:np
    S{i} = str;
    S{i}.W = W(i);
  end

