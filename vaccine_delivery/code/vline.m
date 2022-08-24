function [] = vline(x,varargin)
%simple function to plot a vertical line, make sure to leave this last!
% other properties of the line can also be passed...
 
yb = get(gca,'ylim'); %get the bounds for y
hold on;

if isempty(varargin)
    plot([x,x],[yb(1) yb(2)])
else
    plot([x,x],[yb(1) yb(2)],varargin{:})
end
%set(gca,'ylim',yb); %redefine limits

end