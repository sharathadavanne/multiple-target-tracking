  % Draws the smoothed estimate

  set(gcf,'PaperType','a4');
  set(gcf,'PaperPosition',[0.25 2.5 3.8 3]);
  set(gca,'Position',[0.25 2.5 3.8 3]);

  I1 = find(TC);
  I2 = find(TC == 0);
  
  clf;
  h=plot(Y(1,I1),Y(2,I1),'ko',...
         Y(1,I2),Y(2,I2),'kx');
  set(h(1),'color',[0.3 0.3 0.3]);
  set(h,'markersize',4);
  hold on;
  
  h=plot(X_r(1,:),X_r(2,:),'-',...
         S_EST1{1}(1,:),S_EST1{1}(2,:),'--');

  set(h(1),'color',[0.5 0.5 0.5]);
  set(h(2),'color',[0.0 0.0 0.0]);

  hold on

  % Plot particles
  for k = 1:NT
      for i = 1:size(X_r,2)
          nvis = 50;
          MV = zeros(m,N1);
          PV = zeros(m,m,N1);
          for j=1:N1
              MV(:,j)   = SM{i,j}.M{k};
              PV(:,:,j) = SM{i,j}.P{k};
          end
          %W1 = W_mc1(i,:);
          S1 = SS(i,:);
          tmp = [S1{:}];
          W1 = [tmp.W];
          samp = gmm_rnd(MV,PV,W1,nvis);
          h = plot(samp(1,:),samp(2,:),'r.');
          set(h(1),'color',[0.0 0.0 0.0]);
          set(h,'markersize',2);
          set(h,'linewidth',2);
      end
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