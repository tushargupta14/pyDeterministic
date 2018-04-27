
clear all
clc

global start_temp

y0 = [0 0 0 0 0 0.073];


start_temp = 45;
[t,y] = ode45(@aspgrowthmom,[0 240],y0);


m = y(:,6);

plot(t,m)
hold on
alcp = [0 36 60 120 144 180 204 240];  % Time (sec)
blcp = [0.0730 0.0724 0.0700 0.0524 0.0418 0.0384 0.0332 0.0273]; % Concentration (g/gsolvent)
%for ccp
accp = [0 111.3981 151.1905 176.8335 196.6371 207.8974 218.0545 227.3444 240.0000];
bccp = [0.0730 0.0695 0.0692 0.0604 0.0553 0.0426 0.0425 0.0408 0.0398];
%for occp
acncp = [0 10.2900 30.0000 51.8400 82.3200 101.2500 147.3900 174.9600 205.7700 240.0000];
bcncp = [0.0730 0.0725 0.0682 0.0660 0.0447 0.0400 0.0371 0.0334 0.0323 0.0275];
%for clcp
aclcp = [0    30    60    90   120   140   170   190   210   240];
bclcp = [ 0.0730    0.0729    0.0727    0.0724    0.0706    0.0693    0.0555    0.0488    0.0410    0.0358];

%plot(alcp,blcp,'r*') %linear cooling profile
%plot(accp,bccp,'r*') %controlled cooling profile
plot(acncp,bcncp,'r*') %natural cooling profile
%plot(aclcp,bclcp,'r*') %bilinear cooling profile



(y(end,5)/y(end,4))
(y(end,4)/y(end,3))






