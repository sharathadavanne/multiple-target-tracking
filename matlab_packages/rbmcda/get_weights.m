%GET_WEIGHTS  Extracts an array of weights from particles structures
%
% Syntax:
%    W = get_weights(S)
%
% In:
%   S - NTxN cell array containing particle structures
%
% Out:
%   W - NTxN matrix containing the particle weights
%
%
% Description:
%   Convience function for retrieving the particles weights from
%   particle structures in a regular vector.

% History:
%   27.01.2008  The first official version.
%
% Copyright (C) 2008 Jouni Hartikainen
%
% $Id: get_weights.m $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function W = get_weights(S)
    W = zeros(size(S));
    
    for i = 1:size(S,1)
        for j = 1:size(S,2)
            W(i,j) = S{i,j}.W;
        end
    end
