function g=piecewise_expmap(x,theta,xci)

xcihat          =dinamico_hat(xci);
if theta==0
    g           =diag([1 1 1 1])+x*xcihat;
else
    g           =diag([1 1 1 1])+x*xcihat+...
                 ((1-cos(x*theta))/(theta^2))*xcihat^2+...
                 ((x*theta-sin(x*theta))/(theta^3))*xcihat^3;
end

% eof