function p = Best_poly_fit(x,y)
n = 1;
warn = 0;
lastwarn('')
while warn == 0
    p = polyfit(x,y,n);
    [warnMsg, ~] = lastwarn;
    if ~isempty(warnMsg)
        warn = 1;
        p = polyfit(x,y,n-1);
        n = n-1;
    else
        n = n+1;
    end
end