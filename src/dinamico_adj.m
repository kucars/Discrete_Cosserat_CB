function adj=dinamico_adj(screw)

adj         =zeros(6,6);
adj(1:3,1:3)=dinamico_tilde(screw(1:3));
adj(4:6,1:3)=dinamico_tilde(screw(4:6));
adj(4:6,4:6)=dinamico_tilde(screw(1:3));

% eof