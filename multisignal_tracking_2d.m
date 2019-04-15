%% Tracking unknown number of 2D targets using Rao-Blackwellized particle filtering 
% 
% Read the documentation in - http://becs.aalto.fi/en/research/bayes/rbmcda/
% Specifically this script is modified version of - http://becs.aalto.fi/en/research/bayes/rbmcda/mt_demo.html
%   
% This script reads the frame-wise 2D target location from a CSV file. Each 
% row of the CSV file consists of the time in seconds, and estimated
% location. In this script, as the estimated location, we use the direction 
% of arrival (DOA) of a sound source represented using the azimuth (in 
% 0-360 degree range) and elevation angles (in 0-180 degree range). If more 
% than one source occurs in the same frame, the  consecutive rows will 
% contain the spatial location of each sources with an identical time stamp
%
% The output of the script is the list of DOA tracks in mat file.
% Additionally, this script also visualizes the estimated tracks, and its 
% respective reference. 
%
% For custom setup, you will need to tune PARAMETERS. I have made notes of
% what I have understood of these parameters when tuning myself, and from 
% the original documentation. 

addpath(genpath('matlab_packages/rbmcda/'))
addpath(genpath('matlab_packages/ekfukf/'))

clear all
close all
clc

%% I/O

example = 2; % options - [1, 2]
if example ==1
    % % % Example 1 - Stationary source
    % NOTE: ideally the stationary source location in the input data will be in
    % the range 0:360 and 0:180 for azimuth and elevation. In the test data
    % being used the range is in 0:36, and 0:18, which truly represents the
    % angles 0:10:360 and 0:10:180. Hence we use the spatial_resolution
    % variable to handle this. For your test case, this may not be required.
    % 
    Y_gt = csvread('test_files/stationarly_reference.txt')';
    Y = csvread('test_files/stationary_estimated.txt')';
    spatial_resolution = 10;  
    output_file = 'outputs/test_0_desc_30_100_GT_pks_tracked.mat';
    % Tuned variables
    V_azi = 10;
    V_ele = 1;
    in_sd = 10;
    in_sdn = 100;
    in_cp = 0.25;
elseif example == 2
    % % % Example 2 - Moving source
    Y_gt = csvread('test_files/moving_reference.txt')';
    Y = csvread('test_files/moving_estimated.txt')';
    spatial_resolution = 1;
    output_file = 'outputs/test_0_desc_30_100_GT_pks_tracked.mat';
    % % Tuned variables
    V_azi = 20;
    V_ele = 10;
    in_sd = 5;
    in_sdn = 50;
    in_cp = 0.25;
end


T = Y(1, 1:end);
Y = Y(2:3, 1:end)*spatial_resolution;
%% PARAMETERS

dt =1/33; % Hop length of frames

% Sound signal prior - statistics of the sound sources. 
% We tell the tracker where all a sound source can exist in the spatial 
% grid, and how fast can it move around.    

M0 = [0;0;0;0]; % M0(1, 2) Mean position of the sound sources in the 
                % respective 2D axes. 
                % We are assuming the sources have a mean azimuth and 
                % elevation of 0 and 0 respectively.
                % M0(3, 4) Mean velocity of the sound sources in the 
                % respective 2D axes. 
                % We are assuming the overall velocity average along 
                % azimuth and elevation is 0 and 0 respectively.
                
P0 = diag([360^2 180^2 V_azi^2 V_ele^2]); 
    % P0(0, 1) is the variance of the signal along azimuth and elevation.
    % p0(0)= 360 means that the sound sources can exist anywhere in 0:360
    % degrees in azimuth and similarly 0:180 for elevation P(1)
    %
    % P0(2, 3) is the derivative/velocity of the signal in azimuth and 
    % elevation respectively. For example if a sound source moves 50 degrees 
    % in roughly 40 frames( 8 seconds) then P0(2) has to be set to 50/8 = 6
    % 
    % IDEALLY the velocity along the two axes will have to be tuned on your
    % dataset


% Noise prior
sd = in_sd; % standard deviation of measurement noise - [1 50] range is good
           % Approximately the width of Y axis data spread. 
           % This says that any DOA in the range of gt +/- sd belongs to the
           % same track. Lets say that the gt of the current track is 30 deg
           % and sd is 5, then if the DOA in the next frame is in the range of
           % 30-5 and 30+5 it will be part of the current track.
           % When two sources are very close by, if sd is large it can merge
           % this sources into one, and the final DOA track will be the average
           % of these two tracks.
           
R = diag([sd^2 (sd/2)^2]);

% noise spectral density along x and y axis
% q along with sd decides how smooth the tracked signal is. 
qx = in_sdn;
qy = qx/2;

% Probability of birth and death - Tuning not mandatory
init_birth = 0.1; % value between [0 1] - Prior probability of birth 
alpha_death = 1; % always >= 1; 1 is good  
beta_death = 1; % always >= 1; 1 is good 
                

% Prior probabilities of noise
cd = 1/(360*180);      % Density of noise [max_azi - min_azi] *[max_ele - min_ele]
                       % Assuming the noise is uniformly distributed in the entire spatial grid
                       
cp = in_cp;      % Noise prior - estimate of percentage of noise in the 
                 % measurement data

                 
% Initialize filter
N = 30; %Number of Monte Carlo samples - [10 100] range good   
S = kf_nmcda_init(N,M0,P0,dt);


% Dont need to tune this.    
F = [0 0 1 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];
[A,Q] = lti_disc(F,[],diag([0 0 qx qy]),dt);
H = [1 0 0 0; 0 1 0 0];

%% TRACKING SCRIPT STARTS
% Tracking unknown number of sources
SS = cell(size(Y,2),size(S,2));
for k=1:size(Y,2)
    S = kf_nmcda_predict_dp(S,A,Q,[],[],T(k),alpha_death,beta_death);
    [S,E] = kf_nmcda_update_dp(S,Y(:,k),T(k),H,R,cp,cd,init_birth);

    fprintf('%d/%d: %s\n',k,size(Y,2),E{1});

    SS(k,:) = S;

    %
    % Resample if needed
    %
    S = normalize_weights(S);
    W = get_weights(S);
    n_eff = eff_particles(W);
    if n_eff < N/4
        ind = resampstr(W);
        S = S(ind);
        SS = SS(:,ind);
        W = ones(1,N)/N;
        S = set_weights(S,W);
        fprintf('Resampling done on time step %d\n',k);
    end
end
[FM,FP,SM,SP,Times] = kf_nmcda_collect(SS,A,Q);

%% VISUALIZATON

% Visualize Azimuth and elevation separately
for isazi=1:2
    figure,
    hold on;
    h = plot(Y_gt(1, :), Y_gt(1+isazi,:)*spatial_resolution, 'ro');        
    h=plot(T,Y(isazi,:),'bx');
    set(h,'markersize',2);    

    tmp = [SS{end,:}];
    W = [tmp.W];
    [~,ind] = max(W);
    for j=1:size(SM,1)
        if ~isempty(SM{j,ind})
            t = Times{j,ind};
            m = SM{j,ind};

            h=plot(t,m(isazi,:),'g-');
            set(h(1),'linewidth',2);
        end
    end
    grid on;    
    h = xlabel('Time (seconds) -->');    
    if isazi == 1
        h = ylabel('Azimuth (degree) -->');
        h = title('AZIMUTH - RED : Groundtruth, BLUE: Predicted, GREEN: Tracked');
        ylim([0 360])
    else        
        h = ylabel('Elevation (degree) -->');
        h = title('ELEVATION - RED : Groundtruth, BLUE: Predicted, GREEN: Tracked');
        ylim([0 180])
    end
end

% Let us visualize both azimuth and elevation together on a 2D plot
figure,
h = plot(T,(Y(1,:)*(180/spatial_resolution)+Y(2,:))/spatial_resolution,'b.');
h = title(sprintf(' AZIMUTH*%d+Elevation - RED : Groundtruth, BLUE: Predicted, GREEN: Tracked ', (180/spatial_resolution)));
hold on;
h = plot(Y_gt(1, :), (Y_gt(2,:)*(180/spatial_resolution)+Y_gt(3,:)), 'ro');

tmp = [SS{end,:}];
W = [tmp.W];
[mx,ind] = max(W);
c = 1;
for j=1:size(SM,1)
    if ~isempty(SM{j,ind})
        t = Times{j,ind};
        m = SM{j,ind};        
        
        prev_t = -10;
        new_t = [];        
        new_m = [];
        for k=1:length(t)
            if t(k) ~= prev_t
                new_t = [new_t t(k)];                
                new_m = [new_m m(1:2, k)];
                prev_t = t(k);
            end
        end
        
        h=plot(new_t,(new_m(1,:)*(180/spatial_resolution)+new_m(2,:))/spatial_resolution,'g-','LineWidth',4);
        set(h(1),'linewidth',2);

        tracks(c).t = new_t;        
        tracks(c).m = new_m;
        c = c +1;
    end
end
h = xlabel('time (seconds) -->');
h = ylabel(sprintf('Azimuth*%d+Elevation -->', 180/spatial_resolution));
grid on;

% Dump the tracks to a file. For future processing.
save(output_file, 'tracks');



