%UKF_MCDA_SMOOTH  RTS Smoothing of particles in UKF-RBMCDA algorithm
%
% Syntax:
%   [S,SM] = UKF_MCDA_SMOOTH(S,a,Q,param,same_p)
%
% In:
%   S  - NxNP cell array containing NP particle structures for N time steps.
%   a  - Dynamical model as cell array
%        of size TxN, for each target in each particle,
%        or inline function or name of function in
%        form a(x,i,param), where i is the
%        index of the target x.                    (optional, default A(x,i)*X)
%   Q  - Process noise of discrete model as matrix
%        if common for all targets, or as cell array
%        of size Tx1 for all targets separately.   (optional, default zero)
%  param - Parameters of a.
%  same_p - 1 if the same parameters should be
%           used on every time step                (optional, default 1)
%
% Out:
%   S  - NxNP cell array containing the smoothed particles for each time step.
%   SM - 1xT cell array containing smoothed means for each target as a
%        matrix DxN.
% 
% Description:
%   Perform Extended Kalman Filter prediction step for each target
%   and each association hypothesis particle. The model is
%
%     x_i[k] = a_i(x_i[k-1], q , param),  q ~ N(0,Q_i)
%
%   for each target i. Dynamics a_i() for each target
%   are assumed to have known statistics.
%
% See also:
%   EKF_MCDA_UPDATE, EKF_PREDICT, LTI_DISC, KF_PREDICT

% History:
%    24.01.2008  JH  The first official version.
%
% Copyright (C) 2008 Jouni Hartikainen
%
% $Id:  $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.


function [SS,SM] = ukf_mcda_smooth(S,A,Q,a,WW,param,same_p)
    
    
    % Default values
    if nargin < 4
        a = [];
    end
    if nargin < 5
        WW = [];
    end
    if nargin < 6
        param = [];
    end
    if nargin < 7
        same_p = [];
    end
    
    % Number of particles
    NP = size(S,2);
    % Number of targets
    NT = size(S{1,1}.M,2);

    % Number of data points
    n = size(S,1);
    
    % State space dimensionality
    m = size(S{1,1}.M{1},1);
    
    % Struct for the smoothed particles
    str = struct(...
        'M',{{}},... % Cell array 1xT of T target means
        'P',{{}},... % Cell array 1xT of T target covariances
        'W',0 ...
    );

    % Initialize the cell array 
    SS = cell(n,NP);
    for i = 1:n
        for j = 1:NP
            SS{i,j} = str;
            SS{i,j}.W = S{n,j}.W;
        end
    end
    
    MM = zeros(m,n);
    PP = zeros(m,m,n);

    % Space for smoothed means of each target
    SM = cell(1,NT);
    for i = 1:NT
        SM{i} = zeros(m,n);
    end
    
    for k = 1:NT
        for i=1:NP
            % Form each trajectory for the smoother
            for j = 1:n
                PP(:,:,j) = S{j,i}.P{k};
                MM(:,j) = S{j,i}.M{k};
            end
            
            % Smooth each trajectory
            [SM_i, SP_i] = urts_smooth1(MM,PP,A,Q,a,WW,param,same_p);

            % Save smoothed particles
            for j = 1:n
                SS{j,i}.M{k} = SM_i(:,j);
                SS{j,i}.P{k} = SP_i(:,:,j);
            end

            % Smoothed mean 
            SM{k} = SM{k} + S{n,i}.W*SM_i;
        end
    end