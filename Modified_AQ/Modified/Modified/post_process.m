% This subroutine aims for post processing for Pressure, i.e.,
% to get the average value and the mid point in each element.

% It should be executed seperately after computation being complete.

xmid = 1:200;
pmid = 1:200;
for k = 1:200
    xmid(k) = (X(1,k) + X(2,k))/2;
    pmid(k) = (p(1,k) + p(2,k))/2;
end
plot(xmid, pmid, 'k-', LineWidth=1.5)