%KF_NMCDA_UPDATE_DP  KF/NMCDA Update step with no target death processing
%
% Syntax:
%   [S,EV_STRS] = kf_nmcda_update_dp(S,Y,t,H,R,CP,CD,pb)
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
%   CP    - Prior probability of a measurement being due
%           to clutter.                              (optional, default zero)
%   CD    - Probability density of clutter measurements,
%           which could be for example 1/V, where V is
%           the volume of clutter measurement space. (optional, default 0.01)
%   pb    - Prior probability of birth               (optional, default 0.01)
%
%
% Out:
%         S - Struct array of updated particles
%   EV_STRS - Comment strings of happened events
%
% Description:
%   Perform update step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%   Assumes that target deaths are processed during the prediction step.
%
%
% See also:
%   KF_NMCDA_INIT, KF_NMCDA_PREDICT_DP

% History:
%   26.01.2008  Weights are now inside particle structures
%   19.08.2004  Changed birth/death model to that of article
%   30.10.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: kf_nmcda_update_dp.m,v 1.1 2007-09-24 19:03:04 simosark Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [S,ev_strs] = kf_nmcda_update_dp(S,Y,t,H,R,CP,CD,pb)

%
% Which argument are there
%
if nargin < 6
    CP = [];
end
if nargin < 7
    CD = [];
end
if nargin < 8
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

    %
    % Possible events are:
    %
    %
    evt = {};  % Event descriptions
    evp = [];  % Event priors
    evl = [];  % Event likelhoods
    str = {};  % Structure after each event
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
    str{count}.W = S{i}.W;

    %
    % Loop over associations to targets
    %
    for j=1:length(S{i}.M)

        %
        % Compute update result and likelihood for
        % association to signal j
        %
        [M,P,K,IM,IS,lhood] = ...
            kf_update(S{i}.M{j},S{i}.P{j},Y,H,R);

        %
        % Assocation to target j
        %
        count = count+1;
        evt{count} = sprintf('Target %d',j);
        evb(count) = 0;
        evd(count) = 0;
        evp(count) = (1-pb)*TP0;
        evl(count) = lhood;
        str{count} = S{i};
        str{count}.M{j} = M;
        str{count}.P{j} = P;
        str{count}.T(j) = t;
        str{count}.B = 0;
        str{count}.W = S{i}.W;
    end

    %
    % Initialization of new target
    %
    [M,P,K,IM,IS,lhood] = ...
        kf_update(S{i}.M0,S{i}.P0,Y,H,R);

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
    str{count}.W = S{i}.W;

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
