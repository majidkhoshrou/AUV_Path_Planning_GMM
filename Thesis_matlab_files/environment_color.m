
clc ;
close all ; %clear ;
% format short
% x [0 1]
% y [y_left_d y_right_d]
% z [0 1]
x_min = -5 ;
x_max = run_time+5 ;
% 
y_min = 160 ;
y_max = 350 ;

[x,y] = meshgrid( x_min :2: x_max , y_min : y_max );
y_left = ones(size(y)) * y_min ;
y_right = ones(size(y)) * y_max ;
z_val = [0.0001 .5 1]
% z_val = [ 0.0001 ]
for j = 1 : length(z_val)
    
    z = ones(size(x))* z_val(j);
    in = [z(:),x(:),y(:),y_left(:),y_right(:)] ;
    C = zeros(size(x(:))) ;
    T = zeros(size(x(:))) ;
    D = zeros(size(x(:))) ;
    for i = 1 : size(in,1)
        [C(i),T(i),D(i)] =  T_C_estimate(in(i,:)) ;
    end
    fig_h = surf(y,x,-reshape(D,size(z)),reshape(T,size(z))) ;
% % %      set(fig_h, 'edgecolor','none','box','on'),
    s = sprintf('min T : %2f',min(T));
    text(-12,60,min(-D)+1,s,'color','m') ;
    s = sprintf('max T : %2f',max(T));
    text(35,60,min(-D)+1,s) ;
    set(gca,'FontSize',8,'fontWeight','bold')
    set(findall(gcf,'type','text'),'FontSize',10,'fontWeight','bold')
    hold on    
end


axis([y_min y_max , x_min x_max , -12 0])
colormap;
colorbar;
shading interp ;
h = findobj('tag', 'Colorbar');
set(get(h,'Title'),'String','Temperature \circ C')
set(gca,'FontSize',8,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',8,'fontWeight','bold')
xlabel('y (m)'); ylabel('x (m)') ; zlabel('z (m)') ;
% title('Temperature \circC')
% view(2)
% return

% points = [3 614 1242 1870 2498 3127 3755 4383 ...
%     5012 5640 6268 6897 7525 8153 8782 9410] ;

% points= [1:4000:27000] ;
lin_V1 = ones(length(position_V1),1) ;
lin_V2 = lin_V1 ;
lin_V3 = lin_V2 ;
lin_V4 = lin_V3 ;
lin_V5 = lin_V4 ;
......
x1 = position_V1(:,1) ;
y1 = position_V1(:,2) ;
theta1 = position_V1(:,3) ;
u1 = lin_V1.*cos(theta1) ;
v1 = lin_V1.*sin(theta1) ;
......
x2 = position_V2(:,1) ;
y2 = position_V2(:,2) ;
theta2 = position_V2(:,3) ;
u2 = lin_V2.*cos(theta2) ;
v2 = lin_V2.*sin(theta2)  ;
......
x3 = position_V3(:,1) ;
y3 = position_V3(:,2) ;
theta3 = position_V3(:,3) ;
u3 = lin_V3.*cos(theta3) ;
v3 = lin_V3.*sin(theta3) ;
......
% x4 = position_V4(:,1) ;
% y4 = position_V4(:,2) ;
% theta4 = position_V4(:,3) ;
% u4 = lin_V4.*cos(theta4) ;
% v4 = lin_V4.*sin(theta4) ;
% ......
% x5 = position_V5(:,1) ;
% y5 = position_V5(:,2) ;
% theta5 = position_V5(:,3) ;
% u5 = lin_V5.*cos(theta5) ;
% v5 = lin_V5.*sin(theta5) ;
% close all;
hold on
 plot(y1,x1,'.k',y2,x2,'.r',y3,x3,'.g',y4,x4,'.k',y5,x5,'.y')
legend('','Leader','Follower1','Follower2','Follower3','Follower4','location','southwest')
scale = 0.12 ;

plot(y1(points),x1(points),'S','LineWidth',2, 'MarkerSize',10,'color',[.5 .5 .5],'MarkerFaceColor','w')
quiver(y1(points),x1(points),v1(points),u1(points),scale,'linewidth',2,'color','w')

plot(y2(points),x2(points),'S','LineWidth',2, 'MarkerSize',10, 'color',[.5 .5 .5],'MarkerFaceColor','w')
quiver(y2(points),x2(points),v2(points),u2(points),scale,'linewidth',2,'color','w')

plot(y3(points),x3(points),'S','LineWidth',2, 'MarkerSize',10, 'color',[.5 .5 .5], 'MarkerFaceColor','w')
quiver(y3(points),x3(points),v3(points),u3(points),scale,'linewidth',2,'color','w')

% plot(y4(points),x4(points),'S','LineWidth',2, 'MarkerSize',10, 'color',[.5 .5 .5], 'MarkerFaceColor','w')
% quiver(y4(points),x4(points),v4(points),u4(points),scale,'linewidth',2,'color','w')
% 
% plot(y5(points),x5(points),'S','LineWidth',2, 'MarkerSize',10, 'color',[.5 .5 .5], 'MarkerFaceColor','w')
% quiver(y5(points),x5(points),v5(points),u5(points),scale,'linewidth',2,'color','w')
% 

