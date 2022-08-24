function [] = hline(y,varargin)
%simple function to plot a horizontal line, make sure to leave this last!
% other properties of the line can also be passed...
 
xb = get(gca,'xlim'); %get the bounds for x
hold on;

if isempty(varargin)
    plot([xb(1),xb(2)],[y y])
else
    plot([xb(1),xb(2)],[y y],varargin{:})
end

%set(gca,'xlim',xb); %redefine limits
end