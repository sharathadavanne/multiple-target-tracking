%AZ_DH_DX Measurement derivative function for EKF.with 2d attributes
%
%  \V{x} = [x y vx vy a1 a2]^T
%
%  dh_dx  = -(y-sy) / (x-sx)^2 * 1 / (1 + (y-sy)^2 / (x-sx)^2)
%         = -(y-sy) / ((x-sx)^2 + (y-sy)^2)
%  dh_dy  = 1 / (x-sx) * 1 / (1 + (y-sy)^2 / (x-sx)^2)
%         = (x-sx) / ((x-sx)^2 + (y-sy)^2)
%  da_da1 = 1
%  da_da2 = 1

% Copyright (C) 2003 Simo Särkkä
%
% $Id: az_dh_dx_2a.m,v 1.1 2003/12/10 18:30:12 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function dY = az_dh_dx_2a(x,s)
  dY = [-(x(2)-s(2)) / ((x(1)-s(1))^2 + (x(2)-s(2))^2) ...
	   (x(1)-s(1)) / ((x(1)-s(1))^2 + (x(2)-s(2))^2) 0 0 0 0;...
        0 0 0 0 1 0;
	    0 0 0 0 0 1];
  
