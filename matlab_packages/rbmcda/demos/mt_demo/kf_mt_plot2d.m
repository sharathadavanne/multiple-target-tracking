  % Script for plotting 2D Gaussian signals estimated with RBMCDA

  %
  %       Copyright (C) 2005 Simo Särkkä
  %                     2008 Jouni Hartikainen
  %
  % This software is distributed under the GNU General Public
  % Licence (version 2 or later); please refer to the file
  % Licence.txt, included with the software, for details.

  %
  % Load and collect the trajectories
  %
  [FM,FP,SM,SP,Times] = kf_nmcda_collect(SS,A,Q);

  
  save_figures = 0;

  %
  % Plot the filtered data
  %
  hold off;

  set(gcf,'PaperType','a4');
  set(gcf,'PaperPosition',[0.25 2.5 3.8 3]);
  set(gca,'Position',[0.25 2.5 3.8 3]);
  clf;
  axis off;

  for j=1:4
    XX = [];
    for kk=1:nsteps
      if ~isempty(states{j,kk})
	XX = [XX states{j,kk}];
      end
    end
    if ~isempty(XX)
      h=plot(XX(1,:),XX(2,:),'--');
      set(h,'linewidth',2);
      set(h(1),'color',[0.5 0.5 0.5]);
      hold on;
    end
  end
  h=plot(Y(1,:),Y(2,:),'o');
  set(h,'markersize',2);
  set(h(1),'color',[0.5 0.5 0.5]);
  hold on;

  cols = repmat((0:0.1:0.8)',1,3);

  c = 0;
  tmp = [SS{end,:}];
  W = [tmp.W];
  [mx,ind] = max(W);

  for j=1:size(FM,1)
    if ~isempty(FM{j,ind}) & j~= 7
      t = Times{j,ind};
      m = FM{j,ind};
      c = c + 1;
      if c == size(cols,1)
	c = 1;
      end

      h=plot(m(1,:),m(2,:),'-');
      set(h(1),'color',[0 0 0]);
      set(h(1),'linewidth',2);

      h=plot(m(1,1),m(2,1),'ko');
      set(h(1),'markersize',10);
    end
  end
  if save_figures
      print('-dpsc',sprintf('filtered2-%.2f-%.2f.ps',alpha,beta));
  end
  pause;

  %
  % Plot the number of signals in
  % the most likely particle
  %
  [mx,ind] = max(W);
  clf;
  count = [];
  for k=1:size(SS,1)
    count = [count length(SS{k,ind}.M)];
  end
  size(T)
  size(NT)
  size(count)
  h=plot(T,NT,'--',T,count,'-');
  set(h,'linewidth',2);
  set(h(1),'color',[0.5 0.5 0.5]);
  set(h(2),'color',[0.0 0.0 0.0]);
  legend('True Number of Targets',...
	 'Estimated Number of Targets');
  xlabel('Time');

  axis([min(T) max(T) 0 max(count)+1]);
  fprintf('Number of targets \\alpha=%.2f \\beta=%.2f\n',alpha,beta);
  if save_figures
      print('-dpsc',sprintf('number2-%.2f-%.2f.ps',alpha,beta));
  end
  pause;
  
  %
  % Plot the smoothed data
  %
  hold off;
  clf;
  axis off;

  for j=1:4
    XX = [];
    for kk=1:nsteps
      if ~isempty(states{j,kk})
	XX = [XX states{j,kk}];
      end
    end
    if ~isempty(XX)
      h=plot(XX(1,:),XX(2,:),'--');
      set(h,'linewidth',2);
      set(h(1),'color',[0.5 0.5 0.5]);
      hold on;
    end
  end
  h=plot(Y(1,:),Y(2,:),'ro');
  set(h,'markersize',2);
  set(h(1),'color',[0.5 0.5 0.5]);
  hold on;

  cols = repmat((0:0.1:0.8)',1,3);

  c = 0;
  [mx,ind] = max(W);

  for j=1:size(FM,1)
    if ~isempty(FM{j,ind})       
      t = Times{j,ind};
      m = SM{j,ind};
      c = c + 1;
      if c == size(cols,1)
	c = 1;
      end
      h=plot(m(1,:),m(2,:),'-');
      set(h(1),'color',[0 0 0]);
      set(h(1),'linewidth',2);

      h=plot(m(1,1),m(2,1),'ko');
      set(h(1),'markersize',10);
    end
  end
  if save_figures
      print('-dpsc',sprintf('smoothed2-%.2f-%.2f.ps',alpha,beta));
  end
  pause;

