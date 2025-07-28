function Dh=dwindow(h); %compute the derivative of a given window
%DWINDOW Derive a window.
%	DH=DWINDOW(H) derives a window H.
%
%	Example : 
%	 plot(dwindow(tftb_window(210,'hanning')))
%
%	See also WINDOW.
%============================================
%   window 指令不會對積分 normalize
%    h=gausswin(50,(2.5)^(0.5))
%     =exp(-x.*x*2.5/2)
%      where x=linspace(-1,1,50)
%   結論來說是在整數點上的差分
%   Dh(i)=(h(i+1)-h(i-1))/2 for i=2:n-1;
%   Dh(1) and Dh(n) 需要額外詮釋
%============================================
%	F. Auger, August 1994, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

if (nargin==0),
 error('one parameter required'); 
end;
[hrow,hcol]=size(h); 
if (hcol~=1),
 error('h must have only one column');
end;

Lh=(hrow-1)/2;
step_height=(h(1)+h(hrow))/2;
ramp=(h(hrow)-h(1))/(hrow-1); %平均增長
h2=[0;h-step_height-ramp*(-Lh:Lh).';0]; 
Dh=(h2(3:hrow+2)-h2(1:hrow))/2 + ramp; 
Dh(1)   =Dh(1)   +step_height; 
Dh(hrow)=Dh(hrow)-step_height;

