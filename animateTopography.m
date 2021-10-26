%Script to make a video of a flume experiment

%parameters to set
fileName = 'FlumeExperiment';
skpFrame = 10; %number of scans to skip (1 = no skipping)


%%


% find max and min elevation
zMin = nanmin(Z(:));
zMax = nanmax(Z(:));

skpFrame = 4;

idx = 1;

figure('position',[100 100 1800 400])

v= VideoWriter([fileName '.mp4']);
v.Quality = 100;
open(v);

y = Y(:,1)';
x = X(1,:)';
gxMin = min(x);
gxMax = max(x);
gyMin = min(y);
gyMax = max(y);

for jdx = 1:skpFrame:size(Z,3)

       imageplot(x,y,squeeze(Z(:,:,jdx))')
       view(2)
       ylabel('Flow ⊥ (mm)')
       xlabel('Flow || (mm)')

       shading interp
       axis equal
       axis tight
       xlim([gxMin gxMax])
       ylim([gyMin gyMax])
       box on
       colormap viridis
       caxis([zMin zMax])

       set(gcf,'color','w')
       set(gca,'fontsize',14)
       title(['Duration: ' char(T(jdx).t-T(1).t)])

       drawnow

       frame = getframe(gcf);
       writeVideo(v,frame);



end

close(v);
