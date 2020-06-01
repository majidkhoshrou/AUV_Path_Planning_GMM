
clc; clear all; close all;
Sample_Time = 1 ;

save a a
save b b
save c c 

depth = [a;b;c] ;
plot(depth)

z_one_dive.time = [] ;
z_one_dive.signals.values = depth ;
z_one_dive.signals.dimensions = 1 ;

save z_one_dive z_one_dive


a = 0.001;
b = 0.1;
r = (b-a).*rand(1,1) + a;
