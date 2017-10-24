function ADJg=piecewise_ADJ(x,theta,xci)

adjxci       =dinamico_adj(xci);

if theta==0
    ADJg        =x*diag([1 1 1 1 1 1])+((x^2)/2)*adjxci;
else
    ADJg        =x*diag([1 1 1 1 1 1])+((4-4*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^2))*adjxci+...
                 ((4*x*theta-5*sin(x*theta)+x*theta*cos(x*theta))/(2*theta^3))*adjxci^2+...
                 ((2-2*cos(x*theta)-x*theta*sin(x*theta))/(2*theta^4))*adjxci^3+...
                 ((2*x*theta-3*sin(x*theta)+x*theta*cos(x*theta))/(2*theta^5))*adjxci^4;
end

% eof