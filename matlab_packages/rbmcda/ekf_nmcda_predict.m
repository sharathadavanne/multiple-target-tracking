%EKF_NMCDA_PREDICT  EKF/NMCDA Prediction step
%
% Syntax:
%   S = ekf_nmcda_predict(S,[A,Q,a,AW,param])
%
% Author:
%   Simo Särkkä, 2003
%
% In:
%   S  - Struct array 1xNP of particles
%   A  - Derivative of a() with respect to state 
%        as matrix, or as inline function or name of
%        function in form A(x,param).     (optional, default eye)
%   Q  - Process noise of discrete model     (optional, default zero)
%   a  - Mean prediction function as inline
%        function or name of function in
%        form a(x,param).                (optional, default A(x,i)*X)
%   AW - Derivative of a() with respect to
%        noise q as matrix, or as inline
%        function or name of function in
%        form AW(x,param).                (optional, default eye)
%   param - Parameters of a
% 
% Out:
%   S - Predicted struct array of particles
%
% Description:
%   Perform prediction step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%
% See also:
%   EKF_NMCDA_INIT, EKF_NMCDA_UPDATE, EKF_PREDICT

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

function S = ekf_nmcda_predict(S,A,Q,h,AW,param)
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
    h = [];
  end
  if nargin < 5
    AW = [];
  end
  if nargin < 6
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
	  ekf_predict1(S{i}.M{j},S{i}.P{j},A,Q,h,AW,param);
    end
  end

