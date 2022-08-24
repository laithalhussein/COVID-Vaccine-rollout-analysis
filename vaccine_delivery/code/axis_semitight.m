function axis_semitight(x,a)
% function ax = axis_semitight(a,x)
% This is an elaboration of the built-in matlab command axis('tight') for automatically adjusting axis limits
% This function adjusts the axis limits so as to create a uniform white space boundary around it's contents
% x = boundary width (%)  [default=10] - if scalar then matched x & y boundaries, if 2-elements then treated as [x%, y%]
% a = axis handle  [default = gca]
% ----  MAS 5/17/2019 ----
if nargin <1, x=[10,10]; end
if nargin <2, a=gca; end
if length(x)==1, x = [x,x]; end
axis(a,'tight');
qx = get(a,'xlim'); qx1 = mean(qx) + (1+2*x(1)/100)*0.5*diff(qx)*[-1 1]; set(a,'xlim',qx1);
qy = get(a,'ylim'); qy1 = mean(qy) + (1+2*x(2)/100)*0.5*diff(qy)*[-1 1]; set(a,'ylim',qy1);
%qx,qx1,qy,qy1,keyboard
end  % end-function