%CALCULATE_MEAN  Calculates the mean of the given target from particles 
%
% Syntax:
%   M = calculate_mean(S,T)
%
% In:
%   S - NTxN cell array containing particle structures
%   T - Index of target
%
% Out:
%   M - DxNT matrix containing the mean of the target
%
%
% Description:
%   Convience function for calculating the mean of the given
%   target from particle structures for a single or several time steps.

% History:
%   06.02.2008  The first official version.
%
% Copyright (C) 2008 Jouni Hartikainen
%
% $Id: calculate_mean.m $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function M = calculate_mean(S,T)
% Loop over times teps
M = zeros(size(S{1,1}.M{T},1),size(S,1));
for i = 1:size(S,1)
    for j = 1:size(S,2)
        M(:,i) = M(:,i) + S{i,j}.W*S{i,j}.M{T};
    end
%         tmp = [S{:}];
%        MM1 = [tmp.M]';
% 
%     M = sum([MM1{:}].*repmat(W,size([MM1{:}],1),1),2);
end