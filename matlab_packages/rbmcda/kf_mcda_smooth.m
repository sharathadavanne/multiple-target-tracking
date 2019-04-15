%KF_MCDA_SMOOTH  RTS Smoothing of particles in KF-RBMCDA algorithm
%
% Syntax:
%   [S,SM] = KF_MCDA_SMOOTH(S,A,Q)
%
% In:
%   S  - NxNP cell array containing NP particle structures for N time steps.
%   A  - State transition matrix which can be a same numeric matrix 
%        for every target or a TxNP cell array containing separate matrices
%        for each target in each particle .      
%   Q  - Process noise covariance matrix which can be a same numeric matrix 
%        for every target or a TxNP cell array containing separate matrices
%        for each target in each particle.
%
% Out:
%   S  - NxNP cell array containing the smoothed particles for each time step.
%   SM - 1xT cell array containing smoothed means for each target as a
%        DxN matrix.
% 
% Description:
%   Perform RTS for each target and each association
%   hypothesis particle. 
%
% See also:
%   KF_MCDA_PREDICT, KF_MCDA_UPDATE, RTS_SMOOTH

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


function [SS,SM] = kf_mcda_smooth(S,A,Q,a,AW,param,same_p)
    
    
    % Default values
    if nargin < 4
        a = [];
    end
    if nargin < 5
        AW = [];
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
    
    %
    % Apply defaults
    %
    if isempty(A)
        A = eye(size(S{1}.M{1},1));
    end
    if isempty(Q)
        Q = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                Q{i,j} = zeros(size(S{j}.M{i},1));
            end
        end
    end
    
    %
    % Evaluate matrix A and
    % turn it into cell array.
    %
    if iscell(A)
        % nop
    elseif isnumeric(A)
        tmp = A;
        A = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                A{i,j} = tmp;
            end
        end
    else
        error('A is not of supported form!')
    end
    
    %
    % Turn matrix Q into cell array
    %
    if iscell(Q)
        % nop
    elseif isnumeric(Q)
        tmp = Q;
        Q = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                Q{i,j} = tmp;
            end
        end
    else
        error('Q is not of supported form!')
    end
    
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
            [SM_i, SP_i] = rts_smooth(MM,PP,A{k,i},Q{k,i});

            % Save smoothed particles
            for j = 1:n
                SS{j,i}.M{k} = SM_i(:,j);
                SS{j,i}.P{k} = SP_i(:,:,j);
            end

            % Smoothed mean 
            SM{k} = SM{k} + S{n,i}.W*SM_i;
        end
    end