%KF_NMCDA_UPDATE_SA  KF/NMCDA Update step (single DA per target)
%
% Syntax:
%   [SS,EV_STRS] = kf_nmcda_update_sa(S,Y,time,H,R,pd,CD,pb,alpha,beta)
%
% Author:
%   Simo Särkkä, 2005
%
% In:
%   S     - Struct array 1xNP of particles
%   Y     - DxM measurement vector
%   time  - Time stamp of the measurement
%   H     - Measurement matrix.
%   R     - Measurement noise covariance.
%   pd    - Target detection probability             (optional, default 1)
%   CD    - Probability density of clutter measurements,
%           which could be for example 1/V, where V is
%           the volume of clutter measurement space. (optional, default 0.01)
%   pb    - Prior probability of birth               (optional, default 0.01)
%   alpha - Parameter alpha for the gamma
%           distribution model for time to death     (optional, default 1)
%   beta  - Parameter beta for the gamma
%           distribution model for time to death     (optional, default 10)
%
%
% Out:
%        SS - Predicted array of struct array of particles
%   EV_STRS - Comment strings of happened events
%
% Description:
%   Perform update step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%
%     [...]
%
%
% See also:
%   KF_NMCDA_INIT, KF_NMCDA_PREDICT, EKF_MCDA_UPDATE, KF_PREDICT

% History:
%   26.01.2008  Weights are now inside particles structures.
%   19.08.2004  Changed birth/death model to that of article
%   30.10.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: kf_nmcda_update_sa.m,v 1.1 2005/05/12 12:08:04 ssarkka Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [SS,ev_strs] = kf_nmcda_update_sa(S,Y,time,H,R,pd,CD,pb,alpha,beta)

%
% Which argument are there
%
if nargin < 7
    pd = [];
end
if nargin < 8
    CD = [];
end
if nargin < 9
    pb = [];
end
if nargin < 10
    alpha = [];
end
if nargin < 11
    beta = [];
end

%
% Apply defaults
%
if isempty(pd)
    pd = 1;
end
if isempty(CD)
    CD = 0.01;
end
if isempty(pb)
    pb = 0.01;
end
if isempty(alpha)
    alpha = 10;
end
if isempty(beta)
    beta = 1;
end

%
% Form the detection indicator vectors
%
D = cell(size(S));
for i=1:length(S)
    D{i} = zeros(1,length(S{i}.M));
end

%
% Loop over measurements on the time step
%
SS = cell(1,size(Y,2));
ev_strs = cell(length(S),size(Y,2));
for m=1:size(Y,2)

    r = size(Y,2) - m + 1; % Number of measurements left

    %
    % Loop over particles
    %
    for i=1:length(S)

        T = sum(D{i} == 0); % Number of targets not detected

        %
        % Probabilities of death
        % for the current particle.
        %
        p_death = zeros(1,length(S{i}.M));
        for j=1:length(S{i}.M)
            dt0 = time - S{i}.T(j);
            dt1 = time - S{i}.T(j) + S{i}.dt;
            if dt0 == 0
                p_death(j) = gamma_cdf(dt1,alpha,beta);
            else
                p_death(j) = 1 - (1-gamma_cdf(dt1,alpha,beta)) / ...
                    (1-gamma_cdf(dt0,alpha,beta));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Form the possible events
        % and compute probabilities
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        evt = {};  % Event descriptions
        evp = [];  % Event priors
        evl = [];  % Event likelhoods
        str = {};  % Structure after each event
        count = 0; % Event counter

        %
        % Priors with no target deaths
        %
        pt = 0; % One of targets
        pn = 0; % New target
        pc = 0; % Clutter
        for t=0:min(T,r)
            for n=0:r-t
                pp = bin_pdf(n,pb,r-t) * bin_pdf(t,pd,min(T,r));
                pt = pt + t/r * pp;
                pn = pn + n/r * pp;
                pc = pc + (r-t-n)/r * pp;
            end
        end

        pp0 = [pt pn pc];

        % Clutter
        count = count+1;
        evt{count} = 'Clutter, no deaths';
        evp(count) = prod(1-p_death)*pc;
        evl(count) = CD;
        str{count} = S{i};
        str{count}.B = 0;
        str{count}.D = 0;
        str{count}.det = 0;
        str{count}.W = S{i}.W;

        % Targets
        for j=1:length(S{i}.M)
            if D{i}(j) == 0 % Not yet detected
                %
                % Compute update result and likelihood for
                % association to signal j
                %
                [M,P,K,IM,IS,lhood] = ...
                    kf_update(S{i}.M{j},S{i}.P{j},Y(:,m),H,R);

                %
                % Assocation to target j, no target deaths
                %
                count = count+1;
                evt{count} = sprintf('Target %d, no deaths',j);
                evb(count) = 0;
                evd(count) = 0;
                evp(count) = prod(1-p_death)*pt/T;
                evl(count) = lhood;
                str{count} = S{i};
                str{count}.M{j} = M;
                str{count}.P{j} = P;
                str{count}.T(j) = time;
                str{count}.B = 0;
                str{count}.D = 0;
                str{count}.det = j;
                str{count}.W = S{i}.W;
            end
        end


        % New target
        [M,P,K,IM,IS,lhood] = ...
            kf_update(S{i}.M0,S{i}.P0,Y(:,m),H,R);

        count = count+1;
        j = length(S{i}.M)+1;
        evt{count} = sprintf('New target %d, no deaths',j);
        evp(count) = pn*prod(1-p_death);
        evl(count) = lhood;
        str{count} = S{i};
        str{count}.M{j} = M;
        str{count}.P{j} = P;
        str{count}.T(j) = time;
        str{count}.B = 1;
        str{count}.D = 0;
        str{count}.det = j;
        str{count}.W = S{i}.W;

        %
        % Priors with target death
        %
        pt = 0; % One of targets
        pn = 0; % New target
        pc = 0; % Clutter
        for t=0:min(T-1,r)
            for n=0:r-t
                pp = bin_pdf(n,pb,r-t) * bin_pdf(t,pd,min(T-1,r));
                pt = pt + t/r * pp;
                pn = pn + n/r * pp;
                pc = pc + (r-t-n)/r * pp;
            end
        end

        pp1 = [pt pn pc];

        %
        % Association to clutter, target jj dies
        %
        for jj=1:length(S{i}.M)
            %
            % Death of target jj
            %
            ind = find((1:length(S{i}.M))~=jj); % Remove index jj
            count = count+1;
            evt{count} = sprintf('Clutter, target %d dies',jj);
            if isempty(ind)
                % The only target is dying
                evp(count) = p_death(jj)*pc;
            else
                % Targets still left
                evp(count) = p_death(jj)*prod(1-p_death(ind))*pc;
            end
            evl(count) = CD;
            str{count} = S{i};
            str{count}.M = str{count}.M(ind);
            str{count}.P = str{count}.P(ind);
            str{count}.T = str{count}.T(ind);
            str{count}.B = 0;
            str{count}.D = jj;
            str{count}.det = 0;
            str{count}.W = S{i}.W;
        end

        %
        % Loop over associations to targets
        %
        for j=1:length(S{i}.M)

            if D{i}(j) == 0 % Not yet detected

                %
                % Assocation to target j, target jj dies
                %
                for jj=1:length(S{i}.M)
                    %
                    % We cannot associate to dead target
                    %
                    if jj ~= j
                        %
                        % Death of target jj
                        %
                        ind = find((1:length(S{i}.M))~=jj); % Remove index jj
                        count = count+1;
                        evt{count} = sprintf('Target %d, target %d dies',j,jj);
                        evp(count) = p_death(jj)*prod(1-p_death(ind))*pt/T;
                        evl(count) = lhood;
                        str{count} = S{i};
                        str{count}.M{j} = M;
                        str{count}.P{j} = P;
                        str{count}.T(j) = time;
                        str{count}.M = str{count}.M(ind);
                        str{count}.P = str{count}.P(ind);
                        str{count}.T = str{count}.T(ind);
                        str{count}.B = 0;
                        str{count}.D = jj;
                        str{count}.det = j;
                        str{count}.W = S{i}.W;
                    end
                end
            end
        end

        %
        % Initialization of new target
        %
        [M,P,K,IM,IS,lhood] = ...
            kf_update(S{i}.M0,S{i}.P0,Y(:,m),H,R);

        %
        % Association to new target, target jj dies
        %
        for jj=1:length(S{i}.M)
            %
            % Death of target jj
            %
            ind = find((1:length(S{i}.M))~=jj); % Remove index jj
            count = count+1;
            evt{count} = sprintf('New target %d, target %d dies',j,jj);
            if isempty(ind)
                % The only target is dying
                evp(count) = pn*p_death(jj);
            else
                % Targets still left
                evp(count) = pn*p_death(jj)*prod(1-p_death(ind));
            end
            evl(count) = lhood;
            str{count} = S{i};
            str{count}.B = 1;
            str{count}.D = jj;
            str{count}.M = str{count}.M(ind);
            str{count}.P = str{count}.P(ind);
            str{count}.T = str{count}.T(ind);
            j = length(str{count}.M)+1;
            str{count}.M{j} = M;
            str{count}.P{j} = P;
            str{count}.T(j) = time;
            str{count}.det = j;
            str{count}.W = S{i}.W;
        end

        %
        % Draw sample from importance distribution
        %
        if sum(evp) == 0
            evp = ones(size(evp)) / length(evp);
            warning('Importance distribution is zero.');
            T
            [m size(Y,2) length(evp);pp0;pp1]
        else
            evp = evp./sum(evp);
        end
        imp = evp.*evl;
        imp = imp./sum(imp);

        ev = categ_rnd(imp);  % Event index
        S{i} = str{ev};       % Copy the updated structure
        ev_comment = evt{ev}; % Event description string

        S{i}.W = S{i}.W*evl(ev)*evp(ev)/imp(ev);
        ev_strs{i,m} = ev_comment;

        if str{ev}.det ~= 0
            D{i}(str{ev}.det) = 1; % Detection
        end
    end

    S = normalize_weights(S);
    SS{m} = S;
end
