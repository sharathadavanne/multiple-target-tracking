% Script for Monte Carlo comparison of EKF/MCDA and UKF/MCDA


mcn = 100;

EKF1  = zeros(1,mcn);
EKF2  = zeros(1,mcn);
ERTS1 = zeros(1,mcn);
ERTS2 = zeros(1,mcn);

UKF1  = zeros(1,mcn);
UKF2  = zeros(1,mcn);
URTS1 = zeros(1,mcn);
URTS2 = zeros(1,mcn);

for run = 1:mcn
    run
    ekfukf_tt_demo;
    
    EKF1(run)  = rmse_ekf1;
    EKF2(run)  = rmse_ekf2;
    ERTS1(run) = rmse_erts1;
    ERTS2(run) = rmse_erts2;
    
    UKF1(run)  = rmse_ukf1;
    UKF2(run)  = rmse_ukf2;
    URTS1(run) = rmse_urts1;
    URTS2(run) = rmse_urts2;
end
    
E = zeros(8,mcn);
E(1,:) = EKF1;
E(2,:) = EKF2;
E(3,:) = ERTS1;
E(4,:) = ERTS2;
E(5,:) = UKF1;
E(6,:) = UKF2;
E(7,:) = URTS1;
E(8,:) = URTS2;

    
