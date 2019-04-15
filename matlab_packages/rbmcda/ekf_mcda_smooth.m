%EKF_MCDA_SMOOTH  RTS Smoothing of particles in EKF-RBMCDA algorithm
%
% Syntax:
%   [S,SM] = EKF_MCDA_SMOOTH(S,A,Q,a,AW,param,same_p)
%
% In:
%   S  - NxNP cell array containing NP particle structures for N time steps.
%   A  - Derivative of a() w.r.t. state 
%        as matrix if common for all targets or
%        as cell array of size TxN, for each target in
%        each particle. Can also be a inline function or name of
%        function in form A(x,i,param), where i
%        is the index of the target x.             (optional, default eye())
%   Q  - Process noise of discrete model as matrix
%        if common for all targets, or as cell array
%        of size Tx1 for all targets separately.   (optional, default zero)
%   a  - Dynamical model as cell array
%        of size TxN, for each target in each particle,
%        or inline function or name of function in
%        form a(x,i,param), where i is the
%        index of the target x.                    (optional, default A(x,i)*X)
%   AW - Derivative of a w.r.t. q
%        as matrix if common for all targets or
%        as cell array of size TxN, for each target
%        in each particle. Can also be a inline function or name
%        of function in form AW(x,i,param), where
%        i is the index of the target x.           (optional, default identity)
%  param - Parameters of A,a and AW.
%  same_p - 1 if the same parameters should be
%           used on every time step                (optional, default 1)
%
% Out:
%   S  - NxNP cell array containing the smoothed particles for each time step.
%   SM - 1xT cell array containing smoothed means for each target as a
%        DxN matrix.
% 
% Description:
%   Perform ERTS smoothing for each target
%   and each association hypothesis particle. 
% 
% See also:
%   EKF_MCDA_PREDICT, EKF_MCDA_UPDATE

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


function [SS,SM] = ekf_mcda_smooth(S,A,Q,a,AW,param,same_p)
    
    
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
    if isempty(AW)
        if iscell(Q)
            AW = cell(NT,NP);
            for j=1:NP
                for i=1:NT
                    AW{i,j} = eye(size(S{j}.M{i},1),size(Q{i,j},2));
                end
            end
        else
            AW = eye(size(S{1}.M{1},1),size(Q,2));
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
    elseif isstr(A) | strcmp(class(A),'function_handle')
        tmp = A;
        A = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                A{i,j} = feval(tmp, S{j}.M{i}, i, param);
            end
        end
    else
        tmp = A;
        A = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                A{i,j} = tmp(S{j}.M{i}, i, param);
            end
        end
    end
    
    %
    % Evaluate matrix AW and
    % turn it into cell array.
    %
    if iscell(AW)
        % nop
    elseif isnumeric(AW)
        tmp = AW;
        AW = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                AW{i,j} = tmp;
            end
        end
    elseif isstr(AW) | strcmp(class(AW),'function_handle')
        tmp = AW;
        AW = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                AW{i,j} = feval(tmp, S{j}.M{i}, i, param);
            end
        end
    else
        tmp = AW;
        AW = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                AW{i,j} = tmp(S{j}.M{i}, i, param);
            end
        end
    end
    
    %
    % Turn matrix Q into cell array
    %
    if isnumeric(Q)
        tmp = Q;
        Q = cell(NT,NP);
        for j=1:NP
            for i=1:NT
                Q{i,j} = tmp;
            end
        end
    end
    
    % Turn a into cell array
    if iscell(a)
        % nop
    elseif isempty(a)
        a = cell(NT,NP);
    elseif isnumeric(a)
        error('a cannot be a matrix');
    else
        tmp = a;
        a = cell(NT,NP)
        for j=1:NP
            for i=1:NT
                a{i,j} = tmp;
            end
        end
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
            [SM_i, SP_i] = erts_smooth1(MM,PP,A{k,i},Q{k,i},a{k,i},AW{k,i},param,same_p);

            % Save smoothed particles
            for j = 1:n
                SS{j,i}.M{k} = SM_i(:,j);
                SS{j,i}.P{k} = SP_i(:,:,j);
            end

            % Smoothed mean 
            SM{k} = SM{k} + S{n,i}.W*SM_i;
        end
    end