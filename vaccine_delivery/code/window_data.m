function [new_data]=window_data(X,win,func_handle)
%windows data with a sample size of win
%last window need not be a divisor
%comb_type indicates what type of computation should be done inside each window
%for now, it only takes 'sum' (default) or 'mean'
%Assume X is a vector, not a matrix, though extension shouldnt be too bad

if nargin<3 || isempty(func_handle)
    func_handle = @nansum;
end

if isrow(X), X = X'; end
[N1,N2] = size(X);
if win~=1
    num_rem=rem(N1,win);
    if num_rem>0,
        win_flag=1;
        num_div=floor(N1/win);
        new_data=zeros(num_div+1,N2);
        for p=0:num_div-1,
            if strcmp(func_handle, 'sum')
                new_data(p+1)=nansum(X(p*win+1:(p+1)*win));
            elseif strcmp(func_handle, 'mean')
                new_data(p+1)=nanmean(X(p*win+1:(p+1)*win));
            end
        end
        new_data(:,end)=nanmean(X(:,N1-num_rem+1:N1));
    else
        win_flag=0;
        num_div=(N1/win);
        new_data=zeros(num_div,N2);
        for p=0:num_div-1,
            if strcmp(func_handle, 'sum')
                new_data(p+1)=nansum(X(p*win+1:(p+1)*win));
            elseif strcmp(func_handle, 'mean')
                new_data(p+1)=nanmean(X(p*win+1:(p+1)*win));
            end
        end
    end
else
    new_data=X;
end
return

