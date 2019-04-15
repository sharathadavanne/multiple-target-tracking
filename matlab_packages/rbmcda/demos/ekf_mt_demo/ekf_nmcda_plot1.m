
%
% Collect the trajectories
%
%  [FM,FP,SM,SP,Times] = kf_nmcda_collect(SS,W,A,Q);

%
% Find the number of targets
%
[mx,ind] = max(W);
N = zeros(1,size(SS,1));
for k=1:size(SS,1)
    N(k) = length(SS{k,ind}.M);
end
clf;
plot(N);
title('Number of targets');
pause;

%
% Plot the filtered data
%
hold off;
for j=1:4
    XX = [];
    for kk=1:nsteps
        if ~isempty(states{j,kk})
            XX = [XX states{j,kk}];
        end
    end
    if ~isempty(XX)
        plot(XX(1,:),XX(2,:),'k--');
        hold on;
    end
end
%  h=plot(Y(1,:),Y(2,:),'rx');
%  set(h,'markersize',4);
%  hold on;

cols = repmat((0:0.1:0.8)',1,3);
cols = cols(randperm(size(cols,1)),:);
c = 0;
[mx,ind] = max(W);

for j=1:size(FM,1)
    if ~isempty(FM{j,ind})
        t = Times{j,ind};
        m = FM{j,ind};
        c = c + 1;
        if c == size(cols,1)
            c = 1;
        end
        h=plot(m(1,:),m(2,:),'-');
        set(h(1),'color',cols(c,:));
        set(h(1),'linewidth',2);
        h=plot(m(1,1),m(2,1),'ko');
        set(h(1),'markersize',10);
    end
end
pause;

%
% Plot the smoothed data
%
hold off;
for j=1:4
    XX = [];
    for kk=1:nsteps
        if ~isempty(states{j,kk})
            XX = [XX states{j,kk}];
        end
    end
    if ~isempty(XX)
        plot(XX(1,:),XX(2,:),'k--');
        hold on;
    end
end
%  h=plot(Y(1,:),Y(2,:),'rx');
%  set(h,'markersize',4);
%  hold on;

cols = repmat((0:0.1:0.8)',1,3);
cols = cols(randperm(size(cols,1)),:);
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
        set(h(1),'color',cols(c,:));
        set(h(1),'linewidth',2);
        h=plot(m(1,1),m(2,1),'ko');
        set(h(1),'markersize',10);
    end
end

