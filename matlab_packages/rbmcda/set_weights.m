%SET_WEIGHTS  Sets the weights in particle structures
%
% Syntax:
%    W = set_weights(S,W)
%
% In:
%   S - NTxN cell array containing particle structures
%   W - NTxN matrix containing the particle weights
%
% Out:
%   S - NTxN cell array of updated particle structures
%
%
% Description:
%   Convience function for setting the particles weights in
%   particle structures for the values given in a regular matrix.

% History:
%   27.01.2008  The first official version.
%
% Copyright (C) 2008 Jouni Hartikainen
%
% $Id: set_weights.m $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S = set_weights(S,W)
        
    for i = 1:size(S,1)
        for j = 1:size(S,2)
            S{i,j}.W = W(i,j);
        end
    end
