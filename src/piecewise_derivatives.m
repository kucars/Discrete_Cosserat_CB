function z_punto = piecewise_derivatives(t,z)
global gv

t

% global variable
L           =gv.L;
Eps         =gv.Eps;
Ipsi        =gv.Ipsi;
M           =gv.M;
xci_star    =gv.xci_star;
Gra         =gv.Gra;
dX          =gv.dX;
X           =gv.X;
nsez        =gv.nsez;
npie        =gv.npie;
tact        =gv.tact;
trel        =gv.trel;
Fax         =gv.Fax;
Fay         =gv.Fay;
Faz         =gv.Faz;
Famx        =gv.Famx;
Famy        =gv.Famy;
Famz        =gv.Famz;
Fpx         =gv.Fpx;
Fpy         =gv.Fpy;
Fpz         =gv.Fpz;
Fpmx        =gv.Fpmx;
Fpmy        =gv.Fpmy;
Fpmz        =gv.Fpmz;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% actual solution xci xcidot

Xci              =z(1:6*npie);
Xcidot           =z(6*npie+1:12*npie);             

% inizializzo coefficienti dinamica
genMasM          =zeros(6*npie,6*npie);
genCo1M          =zeros(6*npie,6*npie);
genCo2M          =zeros(6*npie,6*npie);
genForV          =zeros(6*npie,1);
genGraV          =zeros(6*npie,6);
genTipV          =zeros(6*npie,1);

% inizializzo la cnematica precedente
% g_r              =diag([-1 -1 1 1]);         % a testa in giu`
g_r              =[0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % cantilever
Jaco_prec        =diag([1 1 1 1 1 1 zeros(1, 6*(npie-1))]);
g_prec           =diag([1 1 1 1]);
eta_prec         =zeros(6,1);
adetan_prec      =zeros(6*npie,6*npie);

%-------------------------------------------------------------------------
% calcolo le componenti dei coefficienti dinamica

% mass e coriolis 1 del primo pezzo
xci1             =Xci(1:6,:);
xcidot1          =Xcidot(1:6,:);
k1               =xci1(1:3);
theta1           =sqrt(k1'*k1); 

MasX             =zeros(6,6*nsez);
LMasX            =zeros(6,6*nsez);
LRMasX           =zeros(6,6*nsez);
LRCo1X           =zeros(6,6*nsez);
Mas_prec         =zeros(6,6);
LMas_prec        =zeros(6,6);
LRMas_prec       =zeros(6,6);
LRCo1_prec       =zeros(6,6);
for ii=1:nsez
    coAdjg1_here                     =piecewise_coAdjoint(X(ii),theta1,xci1);
    invAdjg1_here                    =piecewise_invAdjoint(X(ii),theta1,xci1);
    intdAdjg1_here                   =piecewise_ADJ(X(ii),theta1,xci1);
        
    % mass
    Mas_here                         =coAdjg1_here*M*invAdjg1_here;
    trapz                            =dX*(Mas_prec+Mas_here)/2;
    MasX(:,6*(ii-1)+1:6*nsez)        =MasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    Mas_prec                         =Mas_here;
    LMas_here                        =intdAdjg1_here'*Mas_here;
    trapz                            =dX*(LMas_prec+LMas_here)/2;
    LMasX(:,6*(ii-1)+1:6*nsez)       =LMasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    LMas_prec                        =LMas_here;
    LRMas_here                       =intdAdjg1_here'*Mas_here*intdAdjg1_here;
    trapz                            =dX*(LRMas_prec+LRMas_here)/2;
    LRMasX(:,6*(ii-1)+1:6*nsez)      =LRMasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    LRMas_prec                       =LRMas_here;
        
    % coriolises 1
    LRCo1_here                       =intdAdjg1_here'*dinamico_coadj(eta_prec+intdAdjg1_here*xcidot1)*Mas_here*intdAdjg1_here;
    trapz                            =dX*(LRCo1_prec+LRCo1_here)/2;
    LRCo1X(:,6*(ii-1)+1:6*nsez)      =LRCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
    LRCo1_prec                       =LRCo1_here;
end
LMasX            =LMasX-repmat(LMasX(:,1:6),[1,nsez]);
LRMasX           =LRMasX-repmat(LRMasX(:,1:6),[1,nsez]);
LRCo1X           =LRCo1X-repmat(LRCo1X(:,1:6),[1,nsez]);
LMas             =LMasX(:,6*(nsez-1)+1:6*nsez);
LRMas            =LRMasX(:,6*(nsez-1)+1:6*nsez);
LRCo1            =LRCo1X(:,6*(nsez-1)+1:6*nsez);

% Actuation load internal load and tip load del primo pezzo
if t<=tact                                      % virata
    Fa1         =[Famx(1);Famy(1);Famz(1);Fax(1);Fay(1);Faz(1)]*t/tact;
end
if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
    Fa1         =[Famx(1);Famy(1);Famz(1);Fax(1);Fay(1);Faz(1)];
end
Fi1             =Eps*(xci1-xci_star)+Ipsi*xcidot1;
Fp1             =[Fpmx(1);Fpmy(1);Fpmz(1);Fpx(1);Fpy(1);Fpz(1)];
    
% Actuation load successiva
if npie ~= 1
    if (t<=tact)                                      % virata
        Fa1_suc     =[Famx(2);Famy(2);Famz(2);Fax(2);Fay(2);Faz(2)]*t/tact;
    end
    if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
        Fa1_suc     =[Famx(2);Famy(2);Famz(2);Fax(2);Fay(2);Faz(2)];
    end
else
    Fa1_suc     =[0;0;0;0;0;0];
end

% aggiorno coefficienti dinamica
invAdjg1_last    =piecewise_invAdjoint(X(nsez),theta1,xci1);
invAdjg1R_last   =blkdiag(invAdjg1_last(1:3,1:3),invAdjg1_last(4:6,4:6));
intdAdjg1_last   =piecewise_ADJ(X(nsez),theta1,xci1);
MasB             =blkdiag(LRMas,zeros(6*(npie-1),6*(npie-1)));
genMasM          =genMasM+Jaco_prec'*MasB*Jaco_prec;
Co1B             =blkdiag(LRCo1,zeros(6*(npie-1),6*(npie-1)));
genCo1M          =genCo1M+Jaco_prec'*Co1B*Jaco_prec;
GraB             =[LMas;zeros(6*(npie-1),6)];
genGraV          =genGraV+Jaco_prec'*GraB*dinamico_Adjoint(g_prec^-1);
% ForB             =[(invAdjg1_last*intdAdjg1_last)'*(Fa1-Fa1_suc);zeros(6*(npie-1),1)]; % pneumatico
% genForV          =genForV+Jaco_prec'*ForB-[L*Fi1;zeros(6*(npie-1),1)];                 % pneumatico
genForV          =genForV+[L*(Fa1-Fi1);zeros(6*(npie-1),1)];  % cavo tip2base
TipB             =[(invAdjg1_last*intdAdjg1_last)'*(invAdjg1R_last*Fp1);zeros(6*(npie-1),1)];
genTipV          =genTipV+Jaco_prec'*TipB;

%----------------------------------------------------------------------
% recursive factors
if npie ~= 1
    Jaco_prec       =blkdiag(invAdjg1_last*intdAdjg1_last,zeros(6*(npie-1),6*(npie-1)))*Jaco_prec+...
                     blkdiag(zeros(6,6),diag([1 1 1 1 1 1]),zeros(6*(npie-2),6*(npie-2)));
    g_prec          =g_prec*piecewise_expmap(X(nsez),theta1,xci1);
    eta_prec        =invAdjg1_last*(eta_prec+intdAdjg1_last*xcidot1);
end

%--------------------------------------------------------------------------
% masses, coriolises 1, coriolises 2 dal secondo pezzo in poi

for jj=2:npie                   
    xcin            =Xci(6*(jj-1)+1:6*(jj-1)+6,:);
    xcidotn         =Xcidot(6*(jj-1)+1:6*(jj-1)+6,:);
    kn              =xcin(1:3);
    thetan          =sqrt(kn'*kn);
    
    MasX            =zeros(6,6*nsez);
    LMasX           =zeros(6,6*nsez);
    RMasX           =zeros(6,6*nsez);
    LRMasX          =zeros(6,6*nsez);
    Co1X            =zeros(6,6*nsez);
    LCo1X           =zeros(6,6*nsez);
    RCo1X           =zeros(6,6*nsez);
    LRCo1X          =zeros(6,6*nsez);
    Co2X            =zeros(6,6*nsez);
    LCo2X           =zeros(6,6*nsez);
    Mas_prec        =zeros(6,6);
    LMas_prec       =zeros(6,6);
    RMas_prec       =zeros(6,6);
    LRMas_prec      =zeros(6,6);
    Co1_prec        =zeros(6,6);
    LCo1_prec       =zeros(6,6);
    RCo1_prec       =zeros(6,6);
    LRCo1_prec      =zeros(6,6);
    Co2_prec        =zeros(6,6);
    LCo2_prec       =zeros(6,6);
    for ii=1:nsez
        coAdjgn_here                     =piecewise_coAdjoint(X(ii),thetan,xcin);
        invAdjgn_here                    =piecewise_invAdjoint(X(ii),thetan,xcin);
        intdAdjgn_here                   =piecewise_ADJ(X(ii),thetan,xcin);
        
        % masses
        Mas_here                         =coAdjgn_here*M*invAdjgn_here;
        trapz                            =dX*(Mas_prec+Mas_here)/2;
        MasX(:,6*(ii-1)+1:6*nsez)        =MasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Mas_prec                         =Mas_here;
        LMas_here                        =intdAdjgn_here'*Mas_here;
        trapz                            =dX*(LMas_prec+LMas_here)/2;
        LMasX(:,6*(ii-1)+1:6*nsez)       =LMasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LMas_prec                        =LMas_here;
        RMas_here                        =Mas_here*intdAdjgn_here;
        trapz                            =dX*(RMas_prec+RMas_here)/2;
        RMasX(:,6*(ii-1)+1:6*nsez)       =RMasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        RMas_prec                        =RMas_here;
        LRMas_here                       =intdAdjgn_here'*Mas_here*intdAdjgn_here;
        trapz                            =dX*(LRMas_prec+LRMas_here)/2;
        LRMasX(:,6*(ii-1)+1:6*nsez)      =LRMasX(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LRMas_prec                       =LRMas_here;
        
        % coriolises 1
        Co1_here                         =dinamico_coadj(eta_prec+intdAdjgn_here*xcidotn)*Mas_here;
        trapz                            =dX*(Co1_prec+Co1_here)/2;
        Co1X(:,6*(ii-1)+1:6*nsez)        =Co1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Co1_prec                         =Co1_here;
        LCo1_here                        =intdAdjgn_here'*Co1_here;
        trapz                            =dX*(LCo1_prec+LCo1_here)/2;
        LCo1X(:,6*(ii-1)+1:6*nsez)       =LCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LCo1_prec                        =LCo1_here;
        RCo1_here                        =Co1_here*intdAdjgn_here;
        trapz                            =dX*(RCo1_prec+RCo1_here)/2;
        RCo1X(:,6*(ii-1)+1:6*nsez)       =RCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        RCo1_prec                        =RCo1_here;
        LRCo1_here                       =intdAdjgn_here'*Co1_here*intdAdjgn_here;
        trapz                            =dX*(LRCo1_prec+LRCo1_here)/2;
        LRCo1X(:,6*(ii-1)+1:6*nsez)      =LRCo1X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LRCo1_prec                       =LRCo1_here;
        
        % coriolises 2
        Co2_here                         =Mas_here*dinamico_adj(intdAdjgn_here*xcidotn);
        trapz                            =dX*(Co2_prec+Co2_here)/2;
        Co2X(:,6*(ii-1)+1:6*nsez)        =Co2X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        Co2_prec                         =Co2_here;
        LCo2_here                        =intdAdjgn_here'*Co2_here;
        trapz                            =dX*(LCo2_prec+LCo2_here)/2;
        LCo2X(:,6*(ii-1)+1:6*nsez)       =LCo2X(:,6*(ii-1)+1:6*nsez)+repmat(trapz,[1,nsez-(ii-1)]);
        LCo2_prec                        =LCo2_here;
    end
    MasX            =MasX-repmat(MasX(:,1:6),[1,nsez]);
    LMasX           =LMasX-repmat(LMasX(:,1:6),[1,nsez]);
    RMasX           =RMasX-repmat(RMasX(:,1:6),[1,nsez]);
    LRMasX          =LRMasX-repmat(LRMasX(:,1:6),[1,nsez]);
    Co1X            =Co1X-repmat(Co1X(:,1:6),[1,nsez]);
    LCo1X           =LCo1X-repmat(LCo1X(:,1:6),[1,nsez]);
    RCo1X           =RCo1X-repmat(RCo1X(:,1:6),[1,nsez]);
    LRCo1X          =LRCo1X-repmat(LRCo1X(:,1:6),[1,nsez]);
    Co2X            =Co2X-repmat(Co2X(:,1:6),[1,nsez]);
    LCo2X           =LCo2X-repmat(Co2X(:,1:6),[1,nsez]);
    Mas             =MasX(:,6*(nsez-1)+1:6*nsez);
    LMas            =LMasX(:,6*(nsez-1)+1:6*nsez);
    RMas            =RMasX(:,6*(nsez-1)+1:6*nsez);
    LRMas           =LRMasX(:,6*(nsez-1)+1:6*nsez);
    Co1             =Co1X(:,6*(nsez-1)+1:6*nsez);
    LCo1            =LCo1X(:,6*(nsez-1)+1:6*nsez);
    RCo1            =RCo1X(:,6*(nsez-1)+1:6*nsez);
    LRCo1           =LRCo1X(:,6*(nsez-1)+1:6*nsez);
    Co2             =Co2X(:,6*(nsez-1)+1:6*nsez);
    LCo2            =LCo2X(:,6*(nsez-1)+1:6*nsez);
    
    % Actuation and internal load
    if t<=tact                                      % virata
        Fan         =[Famx(jj);Famy(jj);Famz(jj);Fax(jj);Fay(jj);Faz(jj)]*t/tact;
    end
    if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
        Fan         =[Famx(jj);Famy(jj);Famz(jj);Fax(jj);Fay(jj);Faz(jj)];
    end
    Fin             =Eps*(xcin-xci_star)+Ipsi*xcidotn;
    Fpn             =[Fpmx(jj);Fpmy(jj);Fpmz(jj);Fpx(jj);Fpy(jj);Fpz(jj)];
    
    % Actuation and internal load successiva
    if jj~=npie
        if t<=tact                                      % virata
            Fan_suc     =[Famx(jj+1);Famy(jj+1);Famz(jj+1);Fax(jj+1);Fay(jj+1);Faz(jj+1)]*t/tact;
        end
        if (t>tact && t<=trel+tact)                     % rilascio/mantenimento  
            Fan_suc     =[Famx(jj+1);Famy(jj+1);Famz(jj+1);Fax(jj+1);Fay(jj+1);Faz(jj+1)];
        end
    else
        Fan_suc     =[0;0;0;0;0;0];
    end
    
    % aggiorno coefficienti dinamica
    invAdjgn_last   =piecewise_invAdjoint(X(nsez),thetan,xcin);
    invAdjgnR_last  =blkdiag(invAdjgn_last(1:3,1:3),invAdjgn_last(4:6,4:6));
    invAdjgprec     =dinamico_Adjoint(g_prec^-1);
    invAdjgprecR    =blkdiag(invAdjgprec(1:3,1:3),invAdjgprec(4:6,4:6));
    intdAdjgn_last  =piecewise_ADJ(X(nsez),thetan,xcin);
    MasB            =blkdiag([repmat(Mas,[jj-1 jj-1]) repmat(RMas,[jj-1 1]); repmat(LMas,[1 jj-1]) LRMas],...
                     zeros(6*(npie-jj),6*(npie-jj)));
    genMasM         =genMasM+Jaco_prec'*MasB*Jaco_prec;
    Co1B            =blkdiag([repmat(Co1,[jj-1 jj-1]) repmat(RCo1,[jj-1 1]); repmat(LCo1,[1 jj-1]) LRCo1],...
                     zeros(6*(npie-jj),6*(npie-jj)));
    genCo1M         =genCo1M+Jaco_prec'*Co1B*Jaco_prec;
    Co2B            =blkdiag([repmat(Co2,[jj-1 jj-1]) zeros(6*(jj-1),6); repmat(LCo2,[1 jj-1]) zeros(6,6)],...
                     zeros(6*(npie-jj),6*(npie-jj)))+MasB*adetan_prec;
    genCo2M         =genCo2M+Jaco_prec'*Co2B*Jaco_prec;
    GraB            =[repmat(Mas,[jj-1 1]);LMas;zeros(6*(npie-jj),6)];
    genGraV         =genGraV+Jaco_prec'*GraB*dinamico_Adjoint(g_prec^-1);
    ForB            =[repmat(invAdjgn_last'*(Fan-Fan_suc),[jj-1 1]);(invAdjgn_last*intdAdjgn_last)'*(Fan-Fan_suc);...
                      zeros(6*(npie-jj),1)];                                                 % pneumatico
    genForV         =genForV+Jaco_prec'*ForB-[zeros(6*(jj-1),1);L*Fin;zeros(6*(npie-jj),1)]; % pneumatico
%     genForV         =genForV+[repmat(L*Fan,[jj-1 1]);L*(Fan-Fin);zeros(6*(npie-jj),1)];       % cavo tip2base
    TipB            =[repmat(invAdjgn_last'*(invAdjgnR_last*invAdjgprecR*Fpn),[jj-1 1]);...
                     (invAdjgn_last*intdAdjgn_last)'*(invAdjgnR_last*invAdjgprecR*Fpn);zeros(6*(npie-jj),1)];
    genTipV         =genTipV+Jaco_prec'*TipB;
    
    %----------------------------------------------------------------------
    % recursive factors
    prec2prec       =[];
    for ii=1:jj-1
        prec2prec   =blkdiag(prec2prec,invAdjgn_last);
    end
    Jaco_prec       =blkdiag(prec2prec,invAdjgn_last*intdAdjgn_last,zeros(6*(npie-jj),6*(npie-jj)))*Jaco_prec+...
                     blkdiag(zeros(6*jj,6*jj),repmat(eye(6),[jj~=npie jj~=npie]),zeros(6*(npie-jj-1),6*(npie-jj-1)));
    g_prec          =g_prec*piecewise_expmap(X(nsez),thetan,xcin);
    ADxin           =intdAdjgn_last*xcidotn;
    eta_prec        =invAdjgn_last*(eta_prec+ADxin);
    prec2prec       =[];
    prec2prec_inv   =[];
    for ii=1:npie
        prec2prec     =blkdiag(prec2prec,invAdjgn_last);
        prec2prec_inv =blkdiag(prec2prec_inv,invAdjgn_last^-1);
        for zz=1:npie
            if (ii+zz == jj)
                adetan_prec(6*(ii-1)+1:6*(ii-1)+6,6*(zz-1)+1:6*(zz-1)+6)  =dinamico_adj(ADxin);
            end
        end
    end
    adetan_prec     =prec2prec*adetan_prec*prec2prec_inv;
end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% calcolo le derivate

Xci_punto       =Xcidot;
Xcidot_punto    =genMasM^-1*(genForV+genGraV*dinamico_Adjoint(g_r^-1)*Gra+genTipV-(genCo1M-genCo2M)*Xcidot);

%----------------------------------------------------------------------
z_punto         =[Xci_punto;Xcidot_punto];
              
% eof