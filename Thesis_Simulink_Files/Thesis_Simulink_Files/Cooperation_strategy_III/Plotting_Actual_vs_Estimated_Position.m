% Plotting of Acutal and Estimated Position

figure(1)
h = plot(Position1(:,1),Position1(:,2),Sensors1(:,3),Sensors1(:,4),Position2(:,1),Position2(:,2),Sensors2(:,3),Sensors2(:,4), Position3(:,1),Position3(:,2),Sensors3(:,3),Sensors3(:,4) );
set(h,'linewidth',1)
grid minor

s = 1; e= length(Sensors1(:,1));
figure(2)
subplot(3,1,1)
stairs(simtime(s:e),[Position1(s:e,1),Position1(s:e,2),Sensors1(s:e,3),Sensors1(s:e,4)])
grid minor
%axis equal
subplot(3,1,2)
stairs(simtime(s:e),[Position2(s:e,1),Position2(s:e,2),Sensors2(s:e,3),Sensors2(s:e,4)])
grid minor
%axis equal
subplot(3,1,3)
stairs(simtime(s:e),[Position3(s:e,1),Position3(s:e,2),Sensors3(s:e,3),Sensors3(s:e,4)])
grid minor
%axis equal


sigma1=mean(((Sensors1(s:e,3)-Position1(s:e,1)).^2+(Sensors1(s:e,4)-Position1(s:e,2)).^2))
sigma2=mean(((Sensors2(s:e,3)-Position2(s:e,1)).^2+(Sensors2(s:e,4)-Position2(s:e,2)).^2))
sigma3=mean(((Sensors3(s:e,3)-Position3(s:e,1)).^2+(Sensors3(s:e,4)-Position3(s:e,2)).^2))
Graphics
sigma1+sigma2+sigma3