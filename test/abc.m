%%%%%%%%%%%%%%%%%%%%%%%
%               Calc for ARSMP
%   (Ag Remote Sensing Mobile Platform)
%
%%%%%%%%%%%%%%%%%%%%%%%
function abc(v,acc)
m = 1000;   % mass--(kg)
g = 9.8;             % Gravity--(N/kg)
vmax = 0.5;        % MAX velocity--(m/s)
accmax = 0.5;     % MAX acceleration--(m/s^2)
%v = vmax;           % velocity--(m/s)
%acc = 0;             % ratio of acceleration--(m/s^2)
f = 0.03;            % coefficient of groud resistence
alpha = 0;           % the angle of slop--(degree)
c = 0;                % coefficient of areodynamic resistence
r= 0.3;               % radius of the wheel--(m)

a = 1;    %
b = 1;    %
L = a+b; %
H = 2;    %

Ff = f*m*g*cos(alpha);  % force of ground resistence--(N)--F=u*cos(a)*mg
Fi = sin(alpha)*m*g;      % force of climb the slope--(N)--F=sin(a)mg
Fw = c*v^2;                % force of areodynamic resistence--(N)F=cv^2
Fj = acc*m;                  % inertia force of acceleretion--(N)--F=am
Fp = Ff;                % force of the plant--(N)--

Ft = Ff + Fi + Fw + Fj + Fp;    %the traction force--(N)
Tt = Ft * r;                   %the traction touque--(Nm)
power = Ft * v;               %the request power rate--(w)

Fz1 = cos(alpha) * m * g * b / L - H / L * (Fi + Fj + Fp);   %
Fz2 = cos(alpha) * m * g * a / L + H / L * (Fi + Fj + Fp);   %

Fx1 = Fz1 / (Fz1 +Fz2) * Ft;  %
Fx2 = Fz2 / (Fz1 +Fz2) * Ft;  %

T1 = Fx1*r;  %
T2 = Fx2*r;  %

X=[Ft,Tt,power,Fz1,Fz2,Fx1,Fx2,T1,T2];
Y=1:9;
bar(Y,X);







