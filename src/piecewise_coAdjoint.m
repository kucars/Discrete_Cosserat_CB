function coAdjg=piecewise_coAdjoint(x,theta,xci)

coadjxci       =dinamico_coadj(xci);

if theta==0
    coAdjg        =diag([1 1 1 1 1 1])+x*coadjxci;
else
    coAdjg        =diag([1 1 1 1 1 1])+((3*sin(x*theta)-x*theta*cos(x*theta))/(2*theta))*coadjxci+...
                 ((4-4*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^2))*coadjxci^2+...
                 ((sin(x*theta)-x*theta*cos(x*theta))/(2*theta^3))*coadjxci^3+...
                 ((2-2*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^4))*coadjxci^4;
end

% eof