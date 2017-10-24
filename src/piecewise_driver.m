close all
format long
clc

global gv

tic

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
disp('Pre-processing')

%-------------------------------------------------------------------------
% beginning of input section

% Geometrical input del braccio siliconico (section)
E        =110e3;                         % [Pa] modulo di young 110e3
eta      =3e3;                           % [Pa*s] viscosità 5e3
Poi      =0;                           % [-] modulo di Poisson 0.5
G        =E/(2*(1+Poi));                 % [Pa] modulo di taglio
R        =10e-3;                         % [m] Raggio braccio 10e-3
L        =62.5e-3;                        % [m] Lunghezza del braccio 250e-3
nsez     =floor(L*2e2+1);                % una sezione per mezzo centimetro floor(L*2e2+1)
X        =linspace(0,L,nsez);            % [m] curvilinear abscissa
A        =pi*R^2;                        % [m^2]
J        =pi*R^4/4;                      % [m^4]
I        =pi*R^4/2;                      % [m^4]

%-------------------------------------------------------------------------
% initial configuration t=0

xci_star    =[0;0;0;1;0;0];

%-------------------------------------------------------------------------
% dinamic parameters

ro_arm      =2000;                             % [Kg/m^3] densità nominale 1080
Gra         =[0;0;0;0;0;0];                    % [m/s^2] vettore gravità [0;0;0;-9.81;0;0]
Eps         =diag([G*I E*J E*J E*A G*A G*A]);  % stifness matrix
Ipsi        =eta*diag([I 3*J 3*J 3*A A A]);    % viscosity matrix eta*diag([I 3*J 3*J 3*A A A]);
M           =ro_arm*diag([I J J A A A]);       % inertia matrix

%-------------------------------------------------------------------------
% numerical setting

time        =20;                     % [s]
nsol        =time*10^2+1;            % una soluzione ogni centisecondo
tspan       =linspace(0,time,nsol);  % [s] time
npie        =4;                      % numero di pezzi
dX          =L/(nsez-1);             % delta X

%-------------------------------------------------------------------------
% actuation load (body coordinate)

tact        =1;                      % [s] tempo torque in dir z o y
trel        =19.0;                    % [s] tempo rilassamento
Fax         =0*[0 0 0 0];              % [N] contraction load
Fay         =0*[0 0 0 0];              % [N] lateral y load 0.1
Faz         =0*[0 0 0 0];              % [N] lateral z load 0.01
Famx        =0*[0 0 0 0];              % [Nm] torsion torque 0.001
Famy        =0*[0 0 0 0];              % [Nm] bending torque
Famz        =0*[0 0 0 0];              % [Nm] bending torque 0.005
%-------------------------------------------------------------------------
% external tip load (base (X=0) coordinate)

Fpx         =0*[0 0 0 0];              % [N] contraction load
Fpy         =0.01*[0 0 0 1];              % [N] lateral y load
Fpz         =0*[0 0 0 0];              % [N] lateral z load
Fpmx        =0*[0 0 0 0];              % [Nm] torsion torque
Fpmy        =0*[0 0 0 0];              % [Nm] bending torque
Fpmz        =0*[0 0 0 0];              % [Nm] bending torque

%-------------------------------------------------------------------------
% osservabili

g           =zeros(4*nsol,4*nsez*npie);
eta         =zeros(6*nsol,nsez*npie);
nstep       =1;

% global variable
gv.ro_arm      =ro_arm;
gv.Gra         =Gra;
gv.L           =L;
gv.X           =X;
gv.R           =R;
gv.xci_star    =xci_star;
gv.A           =A;
gv.Eps         =Eps;
gv.Ipsi        =Ipsi;
gv.M           =M;
gv.nsol        =nsol;
gv.nsez        =nsez;
gv.npie        =npie;
gv.dX          =dX;
gv.time        =time;
gv.tspan       =tspan;
gv.tact        =tact;
gv.trel        =trel;
gv.Fax         =Fax;
gv.Fay         =Fay;
gv.Faz         =Faz;
gv.Famx        =Famx;
gv.Famy        =Famy;
gv.Famz        =Famz;
gv.Fpx         =Fpx;
gv.Fpy         =Fpy;
gv.Fpz         =Fpz;
gv.Fpmx        =Fpmx;
gv.Fpmy        =Fpmy;
gv.Fpmz        =Fpmz;

% osservabili
gv.g           =g;
gv.eta         =eta;
gv.nstep       =nstep;

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% solution initialization 
% sol=[xci*npie xci_dot*npie]

disp('Time-advancing')
myopt          =odeset('RelTol',1e-4,'OutputFcn',@piecewise_observables);

%-------------------------------------------------------------------------
% condizioni iniziali temporali

xci_0          =[0;0;0;1;0;0];
xcidot_0       =[0;0;0;0;0;0];
            
ini_cond        =[repmat(xci_0',[1,npie]) repmat(xcidot_0',[1,npie])];

% integrate
[t,z]          =ode45(@piecewise_derivatives,tspan,ini_cond,myopt);

toc
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% postproc
disp('Post-processing')

nsol=size(z,1);

% piecewise_postproc(t,z)
piecewise_postprocsolo(t,z)

toc

% end
