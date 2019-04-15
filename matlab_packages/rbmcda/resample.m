%RESAMPLE  Resampling for particle filtering
%
%  Syntax:
%    ind = RESAMPLE(W,alg)
%
% Author:
%   Simo Särkkä, 2002
%
% In:
%   W   - Positive weights descibing the relative
%         probabilities of particles. Can be a numeric
%         array or a cell array containing particle structures.
%         The weights do not not need to sum to unity.
%   alg - Algorithm index. Currently only "simple"
%         algorithm as index 1 is defined.      (optional, default 1)
%
% Out:
%   ind - Indices of resampled particles. There may be multiple
%         indices and some of them might be missing, but that is
%         the whole idea of resampling.
%
% Description:
%   Resampling for particle filtering.

% History:
%   23.01.2008  Added support for particle structures
%    3.12.2002  Pure matlab based implementation
%   20.11.2002  The first official version.
%
% Copyright (C) 2002 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: resample.m,v 1.1.1.1 2003/09/15 10:54:34 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function ind = resample(W,alg)

  %
  % Check arguments
  %
  if nargin < 1
    error('Too few arguments');
  end
  if nargin < 2
    alg = [];
  end

  %
  % Apply defaults
  %
  if isempty(alg)
    alg = 1;
  end
  
  %
  % Retrieve the weights if particle structures are used
  %
  if iscell(W) & isfield([W{:}],'W')
      tmp = [W{:}];
      W = [tmp.W]; 
  end
  
  %
  % Resampling
  %
  if alg == 1
    W = W ./ sum(W);
    W = cumsum(W);
    R = sort(rand(1,size(W,2)));
  
    B1 = [zeros(1,size(W,2)) ones(1,size(R,2))];
    [S,I] = sort([W R]);

    B1 = B1(I);                           % Where R's are in joint
                                          % ordering of W,R  
    B2 = [0 B1(1:end-1).*(B1(2:end)==0)]; % Places where R turns into W
    I2 = find(B2);                        % Corresponding indices

    WI1 = cumsum(B2);    % Indices to W's through I2 cumulatively
    WI2 = (B1.*(WI1+1)); % Indices to W's after R's placed over R's
                         % indirectly through I2
    ind = find(WI2);     % We want only nonzero indices
    WI2 = WI2(ind);
    WI3 = I2(WI2);       % Indices to actual W's in S
    ind = I(WI3);        % Indices before sorting
  else
    error(sprintf('Unknown algorithm number %d',alg));
  end

