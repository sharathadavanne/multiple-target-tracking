%UKF_NMCDA_PREDICT  UKF/NMCDA Prediction step
%
% Syntax:
%   S = ukf_nmcda_predict(S,a,Q,param)
%
% Author:
%   Simo Särkkä, 2003
%
% In:
%   S     - 1xN cell array of particle structures
%   a     - Mean prediction function as inline
%           function or name of function in
%           form a(x,param).                (optional, default A(x,i)*X)
%   Q     - Process noise of discrete model     (optional, default zero)
%   param - Parameters of a
% 
% Out:
%   S     - 1xN cell array of predicted particle structures
%
% Description:
%   Perform prediction step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%
% See also:
%   UKF_NMCDA_INIT, UKF_NMCDA_UPDATE, UKF_PREDICT1

% History:
%   09.12.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
%
% $Id: ekf_nmcda_predict.m,v 1.1 2003/12/10 16:40:09 ssarkka Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function S = ukf_nmcda_predict(S,a,Q,param)
  %
  % Check arguments
  %
  if nargin < 2
      a = [];
  end
  if nargin < 3
      Q = [];
  end
  if nargin < 4
      param = [];
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
      % Extended Kalman Filter prediction for the target
      %
      [S{i}.M{j},S{i}.P{j}] = ...
	  ukf_predict1(S{i}.M{j},S{i}.P{j},a,Q,param);
    end
  end

