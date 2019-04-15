%EKF_NMCDA_UPDATE_DP  KF/NMCDA Update step with no target death processing
%
% Syntax:
%   [S,EV_STRS] = kf_nmcda_update_dp(S,Y,t,H,R,h,V,CP,CD,pb,param)
%
% Author:
%   Simo Särkkä, 2003
%
% In:
%   S     - Struct array 1xNP of particles
%   Y     - Dx1 measurement vector
%   t     - Time stamp of the measurement
%   H     - Measurement matrix.
%   R     - Measurement noise covariance.
%   h     - Mean measurement prediction as vector,
%           inline function or name of function in
%           form h(x [,P1,P2,...]). (optional, for default see EKF_UPDATE)
%   V     - Derivative of h() with respect to noise as matrix,
%           inline function or name of function in
%           form V(x [,P1,P2,...]). (optional, for default see EKF_UPDATE)
%   CP    - Prior probability of a measurement being due
%           to clutter.                              (optional, default zero)
%   CD    - Probability density of clutter measurements,
%           which could be for example 1/V, where V is
%           the volume of clutter measurement space. (optional, default 0.01)
%   pb    - Prior probability of birth               (optional, default 0.01)
%   param - Parameters of h. See, for instance, ekf_predict1 or ekf_update1
%           for more details.
%
%
% Out:
%         S - Predicted struct array of particles
%   EV_STRS - Comment strings of happened events
%
% Description:
%   Perform update step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%   Assumes that the target deaths are processed during the prediction
%   step.
%
% See also:
%   EKF_NMCDA_INIT, EKF_NMCDA_PREDICT_DP

% History:
%   25.01.2008  The first official version.
%
% Copyright (C) 2008 Simo Särkkä
%
% $Id: kf_nmcda_update_dp.m,v 1.1 2007-09-24 19:03:04 simosark Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [S,ev_strs] = ekf_nmcda_update_dp(S,Y,t,H,R,h,V,CP,CD,pb,param)
%
% Which argument are there
%
if nargin < 8
    CP = [];
end
if nargin < 9
    CD = [];
end
if nargin < 10
    pb = [];
end

%
% Apply defaults
%
if isempty(CP)
    CP = 0;
end
if isempty(CD)
    CD = 0.01;
end
if isempty(pb)
    pb = 0.01;
end

%
% Loop over particles
%
ev_strs = cell(size(S));
for i=1:length(S)

    %
    % Association priors to targets
    %
    TP0 = (1-CP)/length(S{i}.M);

    % Number of different events
    n_events = 2 + length(S{i}.M);
    % Possible events are:
    evt = cell(1,n_events);   % Event descriptions
    evp = zeros(1,n_events);  % Event priors
    evl = zeros(1,n_events);  % Event likelhoods
    str = cell(1,n_events);   % Structure after each event
    count = 0; % Event counter

    %
    % Association to clutter
    %
    count = count+1;
    evt{count} = 'Clutter';
    evp(count) = (1-pb)*CP;
    evl(count) = CD;
    str{count} = S{i};
    str{count}.B = 0;

    %
    % Loop over associations to targets
    %
    for j=1:length(S{i}.M)

        %
        % Compute update result and likelihood for
        % association to signal j
        %
        [M,P,K,IM,IS,lhood] = ...
            ekf_update1(S{i}.M{j},S{i}.P{j},Y,H,R,h,V,param);

        %
        % Assocation to target j
        %
        count = count+1;
        evt{count} = sprintf('Target %d',j);

        evp(count) = (1-pb)*TP0;
        evl(count) = lhood;
        str{count} = S{i};
        str{count}.M{j} = M;
        str{count}.P{j} = P;
        str{count}.T(j) = t;
        str{count}.B = 0;
    end

    %
    % Initialization of new target
    %
    [M,P,K,IM,IS,lhood] = ...
        ekf_update1(S{i}.M0,S{i}.P0,Y,H,R,h,V,param);

    %
    % Association to new target
    %
    count = count+1;
    j = length(S{i}.M)+1;
    evt{count} = sprintf('New target %d',j);
    evp(count) = pb;
    evl(count) = lhood;
    str{count} = S{i};
    str{count}.M{j} = M;
    str{count}.P{j} = P;
    str{count}.T(j) = t;
    str{count}.B = 1;

    %
    % Draw sample from importance distribution
    %
    evp = evp./sum(evp);
    imp = evp.*evl;
    imp = imp./sum(imp);

    ev = categ_rnd(imp);  % Event index
    S{i} = str{ev};       % Copy the updated structure
    ev_comment = evt{ev}; % Event description string

    S{i}.W = S{i}.W*evl(ev)*evp(ev)/imp(ev);
    ev_strs{i} = ev_comment;
end

S = normalize_weights(S);
