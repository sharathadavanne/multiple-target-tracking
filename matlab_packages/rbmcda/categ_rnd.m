%CATEG_RND  Draws samples from a given one dimensional discrete distribution
%
% Syntax:
%   C = CATEG_RND(P,N)
%
% Author:
%   Simo Särkkä, 2002
%
% In:
%   P - Discrete distribution, which can be a numeric array
%       of probabilities or a cell array of particle structures,
%       whose weights represent the distribution.
%   N - Number of samples (optional, default 1)
%
% Out:
%   C - Samples in a Nx1 vector
%
% Description:
%   Draw random category

% History:
%    3.12.2002  The first version
%
% Copyright (C) 2002 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: categ_rnd.m,v 1.1.1.1 2003/09/15 10:54:35 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function C = categ_rnd(P,N)

  %
  % Check arguments
  %
  if nargin < 1
    error('Too few arguments');
  end
  if nargin < 2
    N = [];
  end

  %
  % Default values
  %
  if isempty(N)
    N = 1;
  end
  
  %
  % If particle structures are given
  %
  if iscell(P) & isfield([P{:}],'W')
      tmp = [P{:}];
      P = [tmp.W]; 
  end
  
  %
  % Draw the categories
  %
  C = zeros(1,N);
  P = P(:) ./ sum(P(:));
  P = cumsum(P);
  for i=1:N
    C(i) = min(find(P > rand));
  end
