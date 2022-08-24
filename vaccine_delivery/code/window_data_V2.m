function [y]=window_data_V2(X,win,F)
%windows data with a sample size of win
%last window need not be a divisor
%func_handle is a function handle for the computation that should be done inside each window
%Assume X is a vector, not a matrix, though extension shouldnt be too bad

if nargin<3 || isempty(F)
    F = @nansum;
end

if isrow(X), X = X'; end
[N1,N2] = size(X);
if win~=1
    num_rem=rem(N1,win);
    if num_rem>0,
        win_flag=1;
        num_div=floor(N1/win);
        y=zeros(num_div+1,N2);
        for p=0:num_div-1,
            y(p+1)=F(X(p*win+1:(p+1)*win));
        end
        y(end)=F(X(N1-num_rem+1:N1));
    else
        win_flag=0;
        num_div=(N1/win);
        y=zeros(num_div,N2);
        for p=0:num_div-1,
            y(p+1)=F(X(p*win+1:(p+1)*win));
        end
    end
else
    y=X;
end
return

