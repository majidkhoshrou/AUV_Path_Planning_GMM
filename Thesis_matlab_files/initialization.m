
clc ; clear all; close all ;

format long ;

clear  getTime1 getTime2 getTime3 
clear DefineHypotheses1 DefineHypotheses2 DefineHypotheses3 
  
BestHypothesis_V1 = [] ;
BestHypothesis_V2 = [] ;
BestHypothesis_V3 = [] ;

save BestHypothesis_V1 BestHypothesis_V1
save BestHypothesis_V2 BestHypothesis_V2
save BestHypothesis_V3 BestHypothesis_V3

Sample_Time = 0.05 ;   % change everywhere!!! ... S-Functions 
% load the depth values 
load('z_path/z_one_dive.mat')
% to create environment with convex combination of the left most and the
% right most values!
y1_d = 250 ;
y_left_most = 1 ;
y_right_most = 500 ;

y2_0 = 220 ;
y1_0 = y1_d ;
y3_0 = 290 ;
theta0 = 0*pi/180 ;

lin_v_1 = 1 ;
eps = .5 ;

run_time = 4000 ;
x_max =  run_time ;
disp('Starting ...')

return

figure; plot(position_V1(:,2),position_V1(:,1),'c',position_V2(:,2),position_V2(:,1),'r',position_V3(:,2),position_V3(:,1),'k')
legend('V1','V2','V3','Location','NorthWest');
xlabel('y');  ylabel('x');


tout(1:abs(length(tout)-length(position_V2)))= [] ;
figure; plot(tout,position_V1(:,1),'c',tout,position_V2(:,1),'r',tout,position_V3(:,1),'k')
legend('V1','V2','V3','Location','SouthEast');title('x');

figure; plot(tout,position_V1(:,2),'c',tout,position_V2(:,2),'r',tout,position_V3(:,2),'k')
legend('V1','V2','V3','Location','SouthEast');title('y');

save position_V1 position_V1
save position_V2 position_V2
save position_V3 position_V3

save CTD1 CTD1
save CTD2 CTD2
save CTD3 CTD3

save K_V1 K_V1
save K_V2 K_V2
save K_V3 K_V3
figure;
K_max = max(max([K_V1,K_V2,K_V3]));
subplot(3,1,1);plot(K_V1);title('K_{V1}'); ylim([0 K_max+1]) ;
subplot(3,1,2);plot(K_V2);title('K_{V2}'); ylim([0 K_max+1]) ;
subplot(3,1,3);plot(K_V3);title('K_{V3}'); ylim([0 K_max+1]);
set(findall(gcf,'type','text'),'FontSize',10,'fontWeight','bold')

save lin_V2 lin_V2
save lin_V3 lin_V3
save ang_V1 ang_V1
save ang_V2 ang_V2
save ang_V3 ang_V3

save tout tout
save err_gmm_12 err_gmm_12
save err_gmm_31 err_gmm_31

save FLAG2 FLAG2
save FLAG3 FLAG3

tout(1:abs(length(tout)-length(err_gmm_31)))= [] ;
figure; plot( tout, err_gmm_12, tout, err_gmm_31 ,'LineWidth',2.2);
% legend('err\_gmm\_12','err\_gmm\_31','Location','SouthEast');
% legend('Left Follower','Right Follower'); 
xlabel('Time (s)')
s= sprintf('Dissimilarity Measure');title(s)
% set(gca,'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',10,'fontWeight','bold')
% grid minor; 
hold on
err_gmmm_goal = 3 ;
gmm_goal_left = err_gmmm_goal*ones(size(tout));
plot(tout,gmm_goal_left,'c','LineWidth',2)
err_gmmm_goal = 3.5 ;
gmm_goal_right = err_gmmm_goal*ones(size(tout));
plot(tout,gmm_goal_right,'m','LineWidth',2)
legend('Left Follower','Right Follower','gmm\_goal\_left','gmm\_goal\_right'); 


bh_dist = [bh_dist_12,bh_dist_13,bh_dist_21,bh_dist_31];
figure:plot(bh_dist);title('Bhattacharyya distance');
legend('bh\_dist\_12','bh\_dist\_13','bh\_dist\_21','bh\_dist\_31');
set(findall(gcf,'type','text'),'FontSize',10,'fontWeight','bold')
grid minor; 

var_dist = [var_dist_12,var_dist_13] ;
figure;plot(var_dist);title('Variational distance');
legend('var\_dist\_12','var\_dist\_13');
set(findall(gcf,'type','text'),'FontSize',10,'fontWeight','bold')
grid minor; 

lin_V1 = ones(size(lin_V2));
figure;
subplot(3,1,1);plot(lin_V1);title('V1'); ylim([0 2]) ;
subplot(3,1,2);plot(lin_V2);title('V2'); ylim([0 2]) ;
subplot(3,1,3);plot(lin_V3);title('V3'); ylim([0 2]);

angular_velocity = [ang_V1,ang_V2,ang_V3] ;
save angular_velocity angular_velocity

lin_V1 = ones(size(lin_V2)) ;
linear_velocity = [lin_V1,lin_V2,lin_V3] ;
save linear_velocity linear_velocity

temperature = [CTD1(:,2),CTD2(:,2),CTD3(:,2)] ;
figure; plot(temperature)
legend('V1','V2','V3')
