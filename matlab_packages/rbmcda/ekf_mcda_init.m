%EKF_MCDA_INIT  EKF/MCDA Data structure initialization
%
% Syntax:
%   [S,W] = ekf_mcda_init(np,M0,P0,[W])
%
% In:
%   np - Number of particles
%   M0 - Cell array of target prior means. The size is TxN if
%        separate is used for each particle, otherwise Tx1.
%   P0 - Cell array of target prior covariances. The size is TxN if
%        separate is used for each particle, otherwise Tx1.
%   W  - Prior particle weights. (optional, default: uniform)
% 
% Out:
%   S - 1xN cell array containing $N$ particle structures.
%
% Description:
%   Initialize data structure for Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%
%   Each element of array S represents one particle, which is a
%   data structure containing the following fields:
%
%      M  : Cell array 1xT of T target means
%      P  : Cell array 1xT of T target covariances
%      W  : Importance weight of the particle
%
%
% See also:
%   EKF_MCDA_PREDICT, EKF_MCDA_UPDATE

% History:
%   17.01.2008  The first official version.
%
% Copyright (C) 2008 Jouni Hartikainen
%
% $Id: ekf_mcda_init.m,v 1.0 2008/01/17 23:33:58 jmjharti Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [S] = ekf_mcda_init(np,M0,P0,W)    
    % If weights are not given use uniform distribution. 
    if nargin < 4
        W = ones(1,np)/np;
    end
    % Same prior for all particles
    if size(M0,2) == 1  
        str = struct(...
            'M',{M0'},... % Cell array 1xT of T target means
            'P',{P0'},... % Cell array 1xT of T target covariances
            'W',0 ...
            );
        
        S = cell(1,np);
        for i=1:np
            S{i} = str;
            S{i}.W = W(i);
        end
    % Different prior for particles
    else       
        S = cell(1,np);
        for i = 1:size(M0,2)
            str = struct(...
                'M',{M0(:,i)'},... % Cell array 1xT of T target means
                'P',{P0(:,i)'},... % Cell array 1xT of T target covariances
                'W',W(i) ...
                );
            S{i} = str;
        end
    end