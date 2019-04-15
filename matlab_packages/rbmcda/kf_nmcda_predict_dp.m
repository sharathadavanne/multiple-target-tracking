%KF_NMCDA_PREDICT_DP  KF/NMCDA Prediction step
%
% Syntax:
%   [S,EV_STRS] = kf_nmcda_predict_dp(S,A,Q,B,u,t,alpha,beta)
%
% In:
%   S - Struct array 1xNP of particles
%   A - Transition matrix of discrete model (optional, default identity)
%   Q - Process noise of discrete model     (optional, default zero)
%   B - Input effect matrix                 (optional, default identity)
%   U - Constant input                      (optional, default empty)
%   t - Time of the prediction step
%   alpha - Parameter alpha for the gamma
%           distribution model for time to death     (optional, default 1)
%   beta  - Parameter beta for the gamma
%           distribution model for time to death     (optional, default 10)
%
% Out:
%   S - Predicted struct array of particles
%   EV_STRS - Comment strings of happened events
%
% Description:
%   Perform prediction step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%   Processes the target deaths with the assumption that only one
%   target can die in each time step.
%
% See also:
%   KF_NMCDA_INIT, KF_NMCDA_UPDATE, KF_PREDICT

% History:
%   30.10.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%
% $Id: kf_nmcda_predict_dp.m,v 1.2 2007-09-25 00:29:20 simosark Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [S,ev_strs] = kf_nmcda_predict_dp(S,A,Q,B,u,t,alpha,beta)

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
if nargin < 7
    alpha = [];
end
if nargin < 8
    beta = [];
end

%
% Apply defaults
%
if isempty(alpha)
    alpha = 10;
end
if isempty(beta)
    beta = 1;
end

%
% Loop over particles
%
ev_strs = cell(size(S));
for i=1:length(S)
    dead = 0;
    %
    % Loop over targets
    %
    for j=1:length(S{i}.M)

        % No target has died yet
        if dead == 0
            %
            % Probability of death
            %
            dt0 = t - S{i}.T(j);
            dt1 = t - S{i}.T(j) + S{i}.dt;
            if dt0 == 0
                p_death = gamma_cdf(dt1,alpha,beta);
            else
                p_death = 1 - (1-gamma_cdf(dt1,alpha,beta)) / ...
                    (1-gamma_cdf(dt0,alpha,beta));
            end

            if (rand < p_death)
                %
                % Target dies
                %
                dead = j;
            end
        end
        if j ~= dead
            % Kalman Filter prediction for the target if alive
            [S{i}.M{j},S{i}.P{j}] = kf_predict(S{i}.M{j},S{i}.P{j},A,Q,B,u);
        end
    end
    % Remove the died target
    if dead ~= 0
        ind = find((1:length(S{i}.M))~=dead); % Remove index
        ev_strs{i} = sprintf('Target %d died ',dead);
        S{i}.M = S{i}.M(ind);
        S{i}.P = S{i}.P(ind);
        S{i}.T = S{i}.T(ind);
        S{i}.D = dead;
        S{i}.B = 0;
    else
        S{i}.D = 0;
    end

end


