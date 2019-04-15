%NORMALIZE_WEIGHTS  Normalizes the weights of the given particles
%
% Syntax:
%    S = normalize_weights(S)
%
% In:
%   S   - 1xN cell array containing particle structures
%
% Out:
%   S   - Updated particle structures
%
%
% Description:
%   Convience function for normalizing the particle weights.

% History:
%   27.01.2008  The first official version.
%
% Copyright (C) 2008 Jouni Hartikainen
%
% $Id: normalize_weights.m $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S = normalize_weights(S)
   tmp = [S{:}];
   W = [tmp.W];
   W = W./sum(W);
   for i = 1:length(S)
     S{i}.W = W(i); 
   end
