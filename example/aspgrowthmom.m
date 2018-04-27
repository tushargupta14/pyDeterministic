function dy = aspgrowthmom(t,y)


dy = zeros(6,1);


%temperature profile


T = 45 - 20*(t/(240))^(1/3); %natural cooling profile



M = 150.14;

% Get supersaturation
% First get solubility, then concentration at t
cstar = (5*10^(-5)*T^2-0.001*T+0.0236);%*10^6/M;     % solubility (g/g sol)
%S = c/cstar;
if y(6)>= cstar
S = (y(6)-cstar)/cstar;
else 
    S = 0;
end

% Calculate moment
%average
kg = exp(-14.5);
g = 1.5475;
kJ1 =exp(20.205);
kJ2 = 0.02165;



Roc = 1568*10^3; %gm/m3

cc = Roc/M; %mol/m3

J = kJ1*(real(y(6)/cstar))*exp(-kJ2*((log(cc/(cstar*10^6/M))^3)/(log(real(y(6)/cstar)))^2))*60;
G = kg*S^g*60*10^6;


kv = pi/6;
           

dy(1)=J;                               % 0 moment 
dy(2)=y(1)*G;                          % 1st moment
dy(3)=2*y(2)*G;                        % 2nd moment
dy(4)=3*y(3)*G;                        % 3rd moment
dy(5)=4*y(4)*G;                        % 4th moment
dy(6)= -3*kv*Roc*G*y(3)*1e-18/10^6;                % concentration
end

