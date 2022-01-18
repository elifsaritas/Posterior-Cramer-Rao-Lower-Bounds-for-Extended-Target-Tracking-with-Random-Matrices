% Implementation of PCRLB recursion for target state (both kinematic and extent states) in 
% E. Sar?ta? and U. Orguner, “Posterior Cramer-Rao Lower Bounds for Extended Target Tracking with Random Matrices,” 
% in Proceedings of 19th International Conference on Information Fusion (FUSION'16), Jul. 2016.

clear;
close all

%% Set Model Parameters

MCRun = 1000;%number of MC runs to use to calculate expectations
T = 1;%sampling time  
runTime = 100;%number of samples in the trajectories


model.T = T;
model.F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
model.B = [T^2/2 0 T 0; 0 T^2/2 0 T]';
model.Q = 0.01*[100 0; 0 100];
model.H = [1 0 0 0; 0 1 0 0];
model.R = 1*[1000 0; 0 1000];
model.s = 1;
model.n0 = 10;
model.n = 20000;
model.angle = 45;    
model.semiMajor = 300;  
model.semiMinor = 100;  
measProcess.mbar = 20;

  
%% Calculate Posterior CRLB
tic
[posCRLB] = ETT_PCRLB(runTime,MCRun,model,measProcess);
toc 


%% Draw the results
figure(1)
suptitle('PCRLB for position')
subplot(2,1,1);
plot(sqrt(squeeze(posCRLB.kinematic(1,1,1:end-1))),'r--','LineWidth',2)
xlabel('time(s)')
ylabel('x(m)')
subplot(2,1,2);
plot(sqrt(squeeze(posCRLB.kinematic(2,2,1:end-1))),'r--','LineWidth',2)
xlabel('time(s)')
ylabel('y(m)')

figure(2)
suptitle('PCRLB for velocity')
subplot(2,1,1);
plot(sqrt(squeeze(posCRLB.kinematic(3,3,1:end-1))),'r--','LineWidth',2)
xlabel('time(s)')
ylabel('v_x(m)')
subplot(2,1,2);
plot(sqrt(squeeze(posCRLB.kinematic(4,4,1:end-1))),'r--','LineWidth',2)
xlabel('time(s)')
ylabel('v_y(m)')

figure(3)
suptitle('PCRLB for extent')
subplot(3,1,1);
plot(sqrt(squeeze(posCRLB.extension(1,1,1:end-1))),'r--','LineWidth',2);xlabel('time(s)')
xlabel('time(s)')
ylabel('[X]_{11}')
subplot(3,1,2);
plot(sqrt(squeeze(posCRLB.extension(2,2,1:end-1))),'r--','LineWidth',2);
xlabel('time(s)')
ylabel('[X]_{21}')
subplot(3,1,3);
plot(sqrt(squeeze(posCRLB.extension(3,3,1:end-1))),'r--','LineWidth',2);
xlabel('time(s)')
ylabel('[X]_{22}')




