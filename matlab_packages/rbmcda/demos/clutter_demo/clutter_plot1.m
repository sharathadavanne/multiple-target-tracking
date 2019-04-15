  % Draws the filtered estimate of the clutter demo

  
  % Set the figure properties
  set(gcf,'PaperType','a4');
  set(gcf,'PaperPosition',[0.25 2.5 3.8 3]);
  set(gca,'Position',[0.25 2.5 3.8 3]); 

  % Find the clutter and non-clutter measurements
  I1 = find(TC);
  I2 = find(TC == 0);
  
  % Plot the measurements
  clf;
  h=plot(Y(1,I1),Y(2,I1),'ko',...
         Y(1,I2),Y(2,I2),'kx');
  set(h(1),'color',[0.3 0.3 0.3]);
  set(h,'markersize',4);
  hold on;


  % Plot the true trajectory and mean estimate
  h=plot(X_r(1,:),X_r(2,:),'-',...
         EST_mc1(1,:),EST_mc1(2,:),'--');
 
  set(h(1),'color',[0.5 0.5 0.5]);
  set(h(2),'color',[0.0 0.0 0.0]);

  hold on
  % Plot particles
  for i = 1:size(X_r,2)
      nvis = 50;
      MV = zeros(m,N1);
      PV = zeros(m,m,N1);
      for j=1:N1
          MV(:,j)   = SS{i,j}.M{1};
          PV(:,:,j) = SS{i,j}.P{1};
      end
      S1 = SS(i,:);
      tmp = [S1{:}];
      W1 = [tmp.W];
      samp = gmm_rnd(MV,PV,W1,nvis);
      h = plot(samp(1,:),samp(2,:),'r.');
      set(h(1),'color',[0.0 0.0 0.0]);
      set(h,'markersize',2);
      set(h,'linewidth',2);
  end

  
  set(h,'markersize',2);
  set(h,'linewidth',2);
  set(gca,'FontSize',4);
  
  legend('Measurement',...
         'Clutter Measurement',...
         'True Signal',...
         'RBMCDA Estimate')
  
  xlim([cx_min cx_max]);
  ylim([cy_min cy_max]);
  
  hold off