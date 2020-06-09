%KF_NMCDA_COLLECT  Collects and smooth trajectories from NMCDA output
%
% Syntax:
%   [FM,FP,SM,SP,Times] = kf_nmcda_collect(SS,A,Q [,silent])
%
% In:
%   SS     - NTxNP cell array of NP particle structures in NT time steps
%   A      - Transition matrix of discrete model (optional, default identity)
%   Q      - Process noise of discrete model     (optional, default zero)
%   silent - Should the function be silent
% 
% Out:
%   FM    - Filtered means of each target in each particle
%   FP    - Filtered covariances of each target in each particle
%   SM    - Smoothed means of each target in each particle
%   SP    - Smoothed covariances of each target in each particle
%   Times - Time instances when each target was alive in each particle
%
% Description:
%   Collects and smooths target trajectories from NMCDA outputs.
%
% See also:
%   NMCDA_INIT, KF_NMCDA_PREDICT, KF_NMCDA_UPDATE, 
%
% History:
%   11.02.2008 JH Added description and cleared the code a bit. 
%   30.10.2003 SS The first official version.
%
% Copyright (C) 2003 Simo S�rkk�
%               2008 Jouni Hartikainen
%
% $Id: kf_nmcda_collect.m, $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.
function [FM,FP,SM,SP,Times,Births,Deaths] = kf_nmcda_collect2(SS,A,Q,T,silent)

if nargin < 2
    A = [];
end
if nargin < 3
    Q = [];
end
if nargin < 5
    silent = [];
end
if isempty(silent)
    silent = 0;
end

FM    = {}; % Filtered means
FP    = {}; % Filtered covariances
SM    = {}; % Smoothed means
SP    = {}; % Smoothed covariances
Times = {}; % Times
Births = [];
Deaths = [];

NP = size(SS,2); % Number of particles
for ind=1:NP
    if ~silent
        fprintf('Processing particle %d/%d\n',ind,NP);
    end

    % Helper variables
    NN=[];
    MM={};
    PP={};
    TT={};
    count=0;
    n=0;
    t=0; 
    for k=1:size(SS,1)
        t = T(k); 
        if ~isempty(SS{k,ind})

            %
            % Death of a target
            %
            if (SS{k,ind}.D ~= 0)

                count = count + 1;
                %
                % Store filtered result
                %
                j = SS{k,ind}.D;
                Times{count,ind} = TT{j};
                FM{count,ind} = MM{j};
                FP{count,ind} = PP{j};

                %
                % Store smoothed result
                %
                [mm,pp] = rts_smooth(MM{j},PP{j},A,Q);
                SM{count,ind} = mm;
                SP{count,ind} = pp;

                %
                % Remove the target
                %
                if j == 1
                    TT = TT(2:end);
                    MM = MM(2:end);
                    PP = PP(2:end);
                elseif j == length(MM)
                    TT = TT(1:end-1);
                    MM = MM(1:end-1);
                    PP = PP(1:end-1);
                else
                    TT = TT([1:j-1 j+1:end]);
                    MM = MM([1:j-1 j+1:end]);
                    PP = PP([1:j-1 j+1:end]);
                end

                TT{n} = [];
                MM{n} = [];
                PP{n} = [];
                n = n - 1;
            end

            %
            % Birth
            %
            if (SS{k,ind}.B > 0)
                n = n + 1;
                MM{n} = [];
                PP{n} = [];
                TT{n} = []; 
            end

            %
            % Store the data for each alive target
            %
            if n > 0
                for j=1:n
                    MM{j} = [MM{j} SS{k,ind}.M{j}];
                    PP{j}(:,:,size(MM{j},2)) = SS{k,ind}.P{j};
                    TT{j} = [TT{j} t];
                end 
            end
            
            NN = [NN n];  
        end
    end

    %
    % Store the alive ones at the last time step
    %
    if n > 0
        for j=1:n
            count = count + 1;

            %
            % Store filtered result
            %
            Times{count,ind} = TT{j};
            FM{count,ind} = MM{j};
            FP{count,ind} = PP{j};

            %
            % Store smoothed result
            %
            [mm,pp] = rts_smooth(MM{j},PP{j},A,Q);
            SM{count,ind} = mm;
            SP{count,ind} = pp;
        end
    end
     
end
if ~silent
    fprintf('Done.\n');
end
