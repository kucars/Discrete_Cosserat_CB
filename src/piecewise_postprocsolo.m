function piecewise_postprocsolo(t,z)
global gv

% global variable
L       =gv.L;
R       =gv.R;
g       =gv.g;
eta     =gv.eta;
X       =gv.X;
nsez    =gv.nsez;
npie    =gv.npie;
nsol    =gv.nsol;

%-------------------------------------------------------------------------
% pre-processing

mkdir('.\LAST RUN\');

% get the solution
tau       =zeros(nsol,npie);
xci       =zeros(nsol,npie);
k         =zeros(nsol,npie);
q         =zeros(nsol,npie);
p         =zeros(nsol,npie);
r         =zeros(nsol,npie);
for ii=1:npie
    tau(:,ii) =z(:,6*(ii-1)+1);
    xci(:,ii) =z(:,6*(ii-1)+2);
    k(:,ii)   =z(:,6*(ii-1)+3);
    q(:,ii)   =z(:,6*(ii-1)+4);
    p(:,ii)   =z(:,6*(ii-1)+5);
    r(:,ii)   =z(:,6*(ii-1)+6);
end
xcidot    =z(:,6*npie+1:12*npie);

%-------------------------------------------------------------------------
% save risults

save('.\LAST RUN\postproc','t','z')
save('.\LAST RUN\long strain','t','q');
save('.\LAST RUN\tras n strain','t','p');
save('.\LAST RUN\tras b strain','t','r');
save('.\LAST RUN\torsione','t','tau');
save('.\LAST RUN\curavtura su n','t','xci');
save('.\LAST RUN\curvatura su b','t','k');
save('.\LAST RUN\strain velocity','t','xcidot');

save('.\LAST RUN\Rototraslation','t','g');
save('.\LAST RUN\velocity','t','eta');

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% plots

for ii=1:npie
    % deformations

    figure
    plot(t,tau(:,ii))
    grid on
    title(strcat('torsion of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('tau [1/m]')
    auxstr  =strcat('.\LAST RUN\torsione',num2str(ii),'.png');
    print('-dpng',auxstr)

    figure
    plot(t,xci(:,ii))
    grid on
    title(strcat('curvature on y of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('xci [1/m]')
    auxstr  =strcat('.\LAST RUN\curvature',num2str(ii),'_on_y.png');
    print('-dpng',auxstr)

    figure
    plot(t,k(:,ii))
    grid on
    title(strcat('curvature on z of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('k [1/m]')
    auxstr  =strcat('.\LAST RUN\curvature',num2str(ii),'_on_z.png');
    print('-dpng',auxstr)

    figure
    plot(t,q(:,ii))
    grid on
    title(strcat('longitudinal strain of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('q [-]')
    auxstr  =strcat('.\LAST RUN\longitudinal_strain',num2str(ii),'.png');
    print('-dpng',auxstr)

    figure
    plot(t,p(:,ii))
    grid on
    title(strcat('tras y strain of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('p [-]')
    auxstr  =strcat('.\LAST RUN\tras_y_strain',num2str(ii),'.png');
    print('-dpng',auxstr)

    figure
    plot(t,r(:,ii))
    grid on
    title(strcat('tras z strain of piece',num2str(ii)))
    xlabel('t [s]')
    ylabel('r [-]')
    auxstr  =strcat('.\LAST RUN\tras_z_strain',num2str(ii),'.png');
    print('-dpng',auxstr)

end

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% Video Poseidrone

mov               =VideoWriter(strcat('.\LAST RUN\Section Dynamics'),'MPEG-4');
mov.FrameRate     =10^2;                                    % fps
open(mov)
scrsz             =get(0,'ScreenSize');
figure('color','w','Position',[scrsz(3)/24 2*scrsz(4)/48 11*scrsz(3)/12 9*scrsz(4)/10])

% figure invariants
ang               =linspace(0,2*pi,180);

for ii=1:nsol                               % per ogni istante
    
    clf
    
    % get rototraslation
    g_now              =g(4*(ii-1)+1:4*(ii-1)+4,:);
    
    % set the graph options
%     set(gca,'CameraPosition',[-0.5*npie*L npie*L -npie*L],...
%         'CameraTarget',[0 0 0],...
%         'CameraUpVector',[1 0 0])          % testa in giu`
    set(gca,'CameraPosition',[0 0 -L],...
        'CameraTarget',[0 0 0],...
        'CameraUpVector',[1 0 0])           % cantilever
    axis equal
    grid on
    hold on
    xlabel('E1 [m]')
    ylabel('E2 [m]')
    zlabel('E3 [m]')
    title(strcat('t= ',num2str(t(ii))))
%     axis ([-1.5*npie*L 2*L -1.5*npie*L 1.5*npie*L -1.5*npie*L 1.5*npie*L]) % testa in giu`
    axis ([-npie*L npie*L 0 1.5*npie*L -L L])  % cantilever

    % disegno la sezione
    sez     =[zeros(1,180) 0;R*sin(ang) 0;R*cos(ang) 0;ones(1,180) 1];
    for zz=1:npie
        for jj=1:nsez
            sez_qui  =g_now(:,4*nsez*(zz-1)+4*(jj-1)+1:4*nsez*(zz-1)+4*(jj-1)+4)*sez;
            plot3(sez_qui(1,:),sez_qui(2,:),sez_qui(3,:),'Color',[1-mod(zz,2),0,mod(zz,2)])
        end
        drawnow
    end
    
    % force drawing
    drawnow
    
    % for movie
    F   = getframe(gcf);
    writeVideo(mov,F);
end

close(mov);

% eof