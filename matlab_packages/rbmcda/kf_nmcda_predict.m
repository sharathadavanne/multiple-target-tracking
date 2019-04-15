%KF_NMCDA_PREDICT  KF/NMCDA Prediction step
%
% Syntax:
%   S = kf_nmcda_predict(S,A,Q,B,u)
%
% Author:
%   Simo Särkkä, 2003
%
% In:
%   S - 1xNP cell array of particle structures
%   A - Transition matrix of dynamic model  (optional, default identity)
%   Q - Process noise of dynamic model      (optional, default zero)
%   B - Input effect matrix                 (optional, default identity)
%   u - Constant input                      (optional, default empty)
% 
% Out:
%   S - 1xNP cell array of updated particle structures
%
% Description:
%   Perform prediction step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%
% See also:
%   KF_NMCDA_INIT, KF_NMCDA_UPDATE, KF_PREDICT

% History:
%   30.10.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%
% $Id: kf_nmcda_predict.m,v 1.2 2003/12/10 16:40:09 ssarkka Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function S = kf_nmcda_predict(S,A,Q,B,u)

  %
  % Check arguments
  %
  if nargin < 2
    A = [];
  end
  if nargin < 3
    Q = [];
  end
  if nargin < 4
    B = [];
  end
  if nargin < 5
    u = [];
  end

  %
  % Loop over particles
  %
  for i=1:length(S)

    %
    % Loop over targets
    %
    for j=1:length(S{i}.M)
      %
      % Kalman Filter prediction for the target
      %
      [S{i}.M{j},S{i}.P{j}] = kf_predict(S{i}.M{j},S{i}.P{j},A,Q,B,u);
    end
  end

