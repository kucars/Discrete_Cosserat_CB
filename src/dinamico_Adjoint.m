function Adj=dinamico_Adjoint(rototras)

Adj         =zeros(6,6);
Adj(1:3,1:3)=rototras(1:3,1:3);
Adj(4:6,1:3)=dinamico_tilde(rototras(1:3,4))*rototras(1:3,1:3);
Adj(4:6,4:6)=rototras(1:3,1:3);

% eof