%AZ_H Azimuth measurement function for EKF.with 2d attributes
%
%  \V{x} = [x y vx vy a1 a2]^T
%
%  h = atan((y-sy) / (x-sx))
%  a = [a1 a2]^T
%  Y = [h;a]
%

% Copyright (C) 2003 Simo Särkkä
%
% $Id: az_h_2a.m,v 1.1 2003/12/10 18:30:12 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function Y = az_h_2a(x,s)
  Y = zeros(size(s,2),size(x,2));

  h = atan2(x(2)-s(2),x(1)-s(1));
  np = find(h>0.5*pi);
  nn = find(h<-0.5*pi);
  if length(nn)>length(np)
    h(np) = h(np)-2*pi;
  else
    h(nn) = h(nn)+2*pi;
  end
  a = [x(end-1);x(end)];
  Y = [h;a];
  
