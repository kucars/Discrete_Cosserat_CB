function se3=dinamico_hat(screw)

se3         =zeros(4,4);
se3(1:3,1:3)=dinamico_tilde(screw(1:3));
se3(1:3,4)  =screw(4:6);

% eof