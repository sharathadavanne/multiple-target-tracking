
%
% Load and collect the trajectories
%
[FM,FP,SM,SP,Times] = kf_nmcda_collect(SS,A,Q);

save_figures = 0;
%
% Plot the data
%
set(gcf,'PaperType','a4');
set(gcf,'PaperPosition',[0.25 2.5 3.8 3]);
set(gca,'Position',[0.25 2.5 3.8 3]);

clf;
h=plot(T,Y,'o');
set(h(1),'color',[0.5 0.5 0.5]);
set(h,'markersize',2);
set(gca,'FontSize',8);
hold on;
for i=1:length(X)
    x = X{i};
    ind = find(diff(C(:,i) ~= 0) ~= 0);
    ind = [ind;size(C,1)];
    k = 1;
    for j=1:length(ind)
        if C(ind(j),i) ~= 0
            h=plot(T(k:ind(j)),x(k:ind(j)));
            set(h(1),'color',[0 0 0]);
            set(h,'linewidth',2);
        end
        k = ind(j)+1;
    end
end
fprintf('This is the measured/true data.\n<Press enter>\n');
%  print('-depsc',sprintf('data1-%.2f-%.2f.eps',alpha,beta));
if save_figures
    print('-dpsc',sprintf('multisignal_data.ps',alpha,beta));
end
pause;

%
% Plot the filtered results in
% best particle
%
clf;
h=plot(T,Y,'o','markersize',1.0);
set(h(1),'color',[0.5 0.5 0.5]);
set(h,'markersize',2);
set(gca,'FontSize',8);
xlabel('Time');
fprintf('Filtered \\alpha=%.2f \\beta=%.2f\n',alpha,beta);
hold on;

cols = repmat((0:0.1:0.5)',1,3);
cols = cols(randperm(size(cols,1)),:);

W = get_weights(S);
[mx,ind] = max(W);
c = 0;
for j=1:size(FM,1)
    if ~isempty(FM{j,ind})
        t = Times{j,ind};
        m = FM{j,ind}(1,:);
        s = sqrt(reshape(FP{j,ind}(1,1,:),1,size(FP{j,ind},3)));
        c = c + 1;
        if c == size(cols,1)
            c = 1;
        end
        h=plot(t,m);
        %      set(h(1),'color',cols(c,:));
        set(h(1),'color',[0 0 0]);
        set(h,'linewidth',2);

        h=plot(t(1),m(1),'ko');
        set(h(1),'markersize',10);
        set(h(1),'color', [0.0 0.0 0.0]);
    end
end
if save_figures
    print('-dpsc',sprintf('multisignal_filtered1-%.2f-%.2f.ps',alpha,beta));
end
pause;


%
% Plot the filtered results in
% each particle
%
for i=1:length(W)
    clf;
    h=plot(T,Y,'o');
    set(h(1),'color',[0.5 0.5 0.5]);
    set(h,'markersize',2);
    set(gca,'FontSize',8);
    xlabel('Time');
    fprintf('(%d/%d) Filtered \\alpha=%.2f \\beta=%.2f\n',i,length(W),alpha,beta);
    hold on;

    cols = repmat((0:0.1:0.5)',1,3);
    cols = cols(randperm(size(cols,1)),:);

    ind = i;
    c = 0;
    for j=1:size(FM,1)
        if ~isempty(FM{j,ind})
            t = Times{j,ind};
            m = FM{j,ind}(1,:);
            s = sqrt(reshape(FP{j,ind}(1,1,:),1,size(FP{j,ind},3)));
            c = c + 1;
            if c == size(cols,1)
                c = 1;
            end
            h=plot(t,m);
            %      set(h(1),'color',cols(c,:));
            set(h(1),'color',[0 0 0]);
            set(h,'linewidth',2);

            h=plot(t(1),m(1),'ko');
            set(h(1),'markersize',10);
            set(h(1),'color', [0.0 0.0 0.0]);
        end
    end

    % Uncomment this if you want to plot each particle
    %if save_figures
    %    print('-dpsc',sprintf('filtered1-%.2f-%.2f-%02d.ps',alpha,beta,i));
    %end
end

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

h=plot(T,NT,'--',T,count,'-');
set(h,'linewidth',2);
set(h(1),'color',[0.5 0.5 0.5]);
set(h(2),'color',[0.0 0.0 0.0]);
legend('True Number of Signals',...
    'Estimated Number of Signals');
xlabel('Time');

axis([min(T) max(T) 0 max(count)+1]);
fprintf('Number of signals \\alpha=%.2f \\beta=%.2f\n',alpha,beta);
if save_figures
    print('-depsc',sprintf('number1-%.2f-%.2f.ps',alpha,beta));
end
pause;

%
% Plot the smoothed result in
% the most likely particle
%
clf;
h=plot(T,Y,'ko');
set(h,'markersize',2);
set(h(1),'color',[0.5 0.5 0.5]);
xlabel('Time');
fprintf('Smoothed \\alpha=%.2f \\beta=%.2f\n',alpha,beta);
hold on;

[mx,ind] = max(W);
c = 0;
for j=1:size(SM,1)
    if ~isempty(SM{j,ind})
        t = Times{j,ind};
        m = SM{j,ind}(1,:);
        s = sqrt(reshape(SP{j,ind}(1,1,:),1,size(SP{j,ind},3)));
        c = c + 1;
        if c == size(cols,1)
            c = 1;
        end
        %      plot(t,m,'-',t,m-2*s,'--',t,m+2*s,'--');

        h=plot(t,m);
        %      set(h(1),'color',cols(c,:));
        set(h(1),'color',[0 0 0]);
        set(h,'linewidth',2);
        h=plot(t(1),m(1),'ko');
        set(h(1),'markersize',10);
        set(h(1),'color', [0.0 0.0 0.0]);
    end
end
if save_figures
    print('-dpsc',sprintf('multisignal_smoothed1-%.2f-%.2f.ps',alpha,beta));
end
pause;
