function status = piecewise_observables(t,z,flag)
global gv

X           =gv.X;
npie        =gv.npie;
nsez        =gv.nsez;

%----------------------------------------------------------------------
% osservabili: posizione (g), velocita (eta)

switch flag
    case []
        for zz=1:length(t)
            g                =gv.g;
            eta              =gv.eta;
            nstep            =gv.nstep;
            
            Xci              =z(1:6*npie,zz);
            Xcidot           =z(6*npie+1:12*npie,zz);
%             g_r              =diag([-1 -1 1 1]);         % a testa in giu`
            g_r              =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % cantilever
            g_prec           =diag([1 1 1 1]);
            eta_prec         =zeros(6,1);

            for jj=1:npie
                xcin             =Xci(6*(jj-1)+1:6*(jj-1)+6,:);
                xcidotn          =Xcidot(6*(jj-1)+1:6*(jj-1)+6,:);
                kn               =xcin(1:3);
                thetan           =sqrt(kn'*kn);
    
                % kinematics
                for ii=1:nsez
                    invAdjgn_here    =piecewise_invAdjoint(X(ii),thetan,xcin);
                    intdAdjgn_here   =piecewise_ADJ(X(ii),thetan,xcin);
                    g(4*(nstep-1)+1:4*(nstep-1)+4,4*(jj-1)*nsez+4*(ii-1)+1:4*(jj-1)*nsez+4*(ii-1)+4)...
                        =g_r*g_prec*piecewise_expmap(X(ii),thetan,xcin);
                    eta(6*(nstep-1)+1:6*(nstep-1)+6,(jj-1)*nsez+ii)...
                        =invAdjgn_here*(eta_prec+intdAdjgn_here*xcidotn);
                end
    
                % recursive factors
                invAdjgn_last   =piecewise_invAdjoint(X(nsez),thetan,xcin);
                intdAdjgn_last  =piecewise_ADJ(X(nsez),thetan,xcin);
                g_prec          =g_prec*piecewise_expmap(X(nsez),thetan,xcin);
                ADxin           =intdAdjgn_last*xcidotn;
                eta_prec        =invAdjgn_last*(eta_prec+ADxin);
            end
            
            gv.g             =g;
            gv.eta           =eta;
            gv.nstep         =nstep;
            gv.nstep         =nstep+1;
        end
end
  
status      =0;

% eof