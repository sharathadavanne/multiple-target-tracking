%EKF_NMCDA_PREDICT  EKF/NMCDA Update step
%
% Syntax:
%   [S,EV_STRS] = ekf_nmcda_update(S,Y,t,H,R,h,V,CP,CD,
%                                    a_birth,l_birth,
%                                    a_death,l_death,param)
%
% In:
%   S  - Struct array 1xNP of particles
%   Y  - Dx1 measurement vector
%   t  - Time stamp of the measurement
%   H  - Derivative of h() with respect to state as matrix,
%        inline function or name of function in
%        form H(x,param).
%   R  - Measurement noise covariance.
%   h  - Mean measurement prediction as vector,
%        inline function or name of function in
%        form h(x,param). (optional, for default see EKF_UPDATE)
%   V  - Derivative of h() with respect to noise as matrix,
%        inline function or name of function in
%        form V(x,param). (optional, for default see EKF_UPDATE)
%   CP - Prior probability of a measurement being due
%        to clutter.                               (optional, default zero)
%   CD - Probability density of clutter measurements,
%        which could be for example 1/V, where V is
%        the volume of clutter measurement space.  (optional, default 0.01)
%   pb    - Prior probability of birth             (optional, default 0.01)
%   alpha - Parameter alpha for the gamma
%           distribution model for time to death   (optional, default 1)
%   beta  - Parameter beta for the gamma
%           distribution model for time to death   (optional, default 10)
%   param - Parameters of H, h and V. See, for instance, ekf_predict1 or ekf_update1
%           for more details.
%
% Out:
%         S - Predicted struct array of particles
%   EV_STRS - Comment strings of happened events
%
% Description:
%   Perform update step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%
% See also:
%   EKF_NMCDA_INIT, EKF_NMCDA_PREDICT, EKF_UPDATE

% History:
%   11.02.2008  Added support for the new particle structure
%   24.08.2004  Changed birth/death model to that of article
%   09.12.2003  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: ekf_nmcda_update.m,v 1.3 2005/05/12 12:07:59 ssarkka Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [S,ev_strs] = ekf_nmcda_update(S,Y,t,H,R,h,V,CP,CD,pb,alpha,beta,param)

%
% Which argument are there
%
if nargin < 6
    h = [];
end
if nargin < 7
    V  = [];
end
if nargin < 8
    CP = [];
end
if nargin < 9
    CD = [];
end
if nargin < 10
    pb = [];
end
if nargin < 11
    alpha = [];
end
if nargin < 12
    beta = [];
end
if nargin < 13
    param = [];
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

    %
    % Association priors to targets
    %
    if length(S{i}.M)>0
        TP0 = (1-CP)/length(S{i}.M);     % No deaths
    end
    if length(S{i}.M)>1
        TP1 = (1-CP)/(length(S{i}.M)-1); % One death
    end
    TP2 = (1-CP)/(length(S{i}.M)+1);   % One birth

    %
    % Probabilities of birth and death
    % for the current particle.
    %
    p_birth = pb;
    p_death = zeros(1,length(S{i}.M));
    for j=1:length(S{i}.M)
        dt0 = t - S{i}.T(j);
        dt1 = t - S{i}.T(j) + S{i}.dt;
        if dt0 == 0
            p_death(j) = gamma_cdf(dt1,alpha,beta);
        else
            p_death(j) = 1 - (1-gamma_cdf(dt1,alpha,beta)) / ...
                (1-gamma_cdf(dt0,alpha,beta));
        end
    end

    %
    % Possible events are:
    %
    %
    %evt = {};  % Event descriptions
    %evp = [];  % Event priors
    %evl = [];  % Event likelhoods
    %str = {};  % Structure after each event

    n_events = 2*(1 + length(S{i}.M)) + ...
        (length(S{i}.M)*(length(S{i}.M)));

    evt = cell(1,n_events); % Event descriptions
    evp = zeros(1,n_events);  % Event priors
    evl = zeros(1,n_events);  % Event likelhoods
    str = cell(1,n_events);  % Structure after each event
    count = 0; % Event counter

    %
    % Association to clutter, no target deaths
    %
    count = count+1;
    evt{count} = 'Clutter, no deaths';
    evp(count) = (1-p_birth)*prod(1-p_death)*CP;
    evl(count) = CD;
    str{count} = S{i};
    str{count}.B = 0;
    str{count}.D = 0;
    str{count}.W = S{i}.W;

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
            evp(count) = (1-p_birth)*p_death(jj)*CP;
        else
            % Targets still left
            evp(count) = (1-p_birth)*p_death(jj)*prod(1-p_death(ind))*CP;
        end
        evl(count) = CD;
        str{count} = S{i};
        str{count}.M = str{count}.M(ind);
        str{count}.P = str{count}.P(ind);
        str{count}.T = str{count}.T(ind);
        str{count}.B = 0;
        str{count}.D = jj;
        str{count}.W = S{i}.W;
    end

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
        % Assocation to target j, no target deaths
        %
        count = count+1;
        evt{count} = sprintf('Target %d, no deaths',j);
        evp(count) = (1-p_birth)*prod(1-p_death)*TP0;
        evl(count) = lhood;
        str{count} = S{i};
        str{count}.M{j} = M;
        str{count}.P{j} = P;
        str{count}.T(j) = t;
        str{count}.B = 0;
        str{count}.D = 0;
        str{count}.W = S{i}.W;

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
                evp(count) = (1-p_birth)*p_death(jj)*prod(1-p_death(ind))*TP1;
                evl(count) = lhood;
                str{count} = S{i};
                str{count}.M{j} = M;
                str{count}.P{j} = P;
                str{count}.T(j) = t;
                str{count}.M = str{count}.M(ind);
                str{count}.P = str{count}.P(ind);
                str{count}.T = str{count}.T(ind);
                str{count}.B = 0;
                str{count}.D = jj;
                str{count}.W = S{i}.W;
            end
        end
    end

    %
    % Initialization of new target
    % TODO: Invent better initialization method
    %
    m = gauss_rnd(S{i}.M0,S{i}.P0,1);
    [M,P,K,IM,IS,lhood] = ...
        ekf_update1(m,S{i}.P0,Y,H,R,h,V,param);

    %
    % Association to new target, no target deaths
    %
    count = count+1;
    j = length(S{i}.M)+1;
    evt{count} = sprintf('New target %d, no deaths',j);
    evp(count) = p_birth*prod(1-p_death)*TP2;
    evl(count) = lhood;
    str{count} = S{i};
    str{count}.M{j} = M;
    str{count}.P{j} = P;
    str{count}.T(j) = t;
    str{count}.B = 1;
    str{count}.D = 0;
    str{count}.W = S{i}.W;

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
            evp(count) = p_birth*p_death(jj)*TP0;
        else
            % Targets still left
            evp(count) = p_birth*p_death(jj)*prod(1-p_death(ind))*TP0;
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
        str{count}.T(j) = t;
        str{count}.W = S{i}.W;
    end

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
