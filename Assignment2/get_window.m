function x_window=get_window(x,order,n)
%     x: signal
%     order: lag order
%     n: time index

    x_window=zeros(order+1,1);
    
    for j=1:order+1
        if((n-j) >= 0) 
            x_window(j,1) = x(n-j+1);
        else
            x_window(j,1) = 0; % clip here - can't consider neg. time
        end
    end
end