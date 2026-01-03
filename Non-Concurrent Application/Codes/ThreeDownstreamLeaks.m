%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LPV Kalman Filter Design for Quasi-LPV Systems 
% with Unmeasurable Gain Scheduling Functions 
% Application to leak diagnosis 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean
clear all
clear classes
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       YALMIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('Insert YALMIP direction here.'));

%Example:
addpath(genpath('C:\Users\adria\Desktop\DOCUMENTOS TESIS\solvers\YALMIP-master'));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SOLVER/S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('Insert solver/s direction here.'));

%Example:
addpath(genpath('C:\Program Files\Mosek'));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       IMPORT DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Datos = importdata('Insert the data base direction here\datos_reales_2.csv');

%Example:
Datos = importdata('E:\Tesis_Codigos_DatosReales\datos_reales_2.csv');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = Datos.data(:,1)*1e-03;
QinFiltrado = Datos.data(:,2)*1e-03;
HinFiltrado = Datos.data(:,5);
HoutFiltrado = Datos.data(:,6);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DATA 
%%                  SINGLE LEAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = QoutFiltrado((1:70000),1);
QinFiltrado = QinFiltrado((1:70000),1);
HoutFiltrado = HoutFiltrado((1:70000),1);
HinFiltrado = HinFiltrado((1:70000),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       TIME VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;
fs = 1/dt;
muestras = length(QinFiltrado);
tiempo = [dt:dt:(muestras/fs)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       PIPELINE PHYSICAL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lr=68.147;                            %Length (m)
g=9.81;                              %Gravity (m/s^2)
D=0.06271;%                          %Inner Diameter (m)
A=pi*(D/2)^2;                        %Inner Area (m^2)
e=0.013095;                           %Wall thickness (m)
eps=7e-6;                        %Roughness (m)
E=6.8646e+08; %% How can the temperature affect the performance of a classical pipeline model when plastic pipes are used?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Average Temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TemperaturaFiltrada = Datos.data(:, 9);
Temp_av = sum(TemperaturaFiltrada((1:66000),1))/length(TemperaturaFiltrada((1:66000),1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Fluid properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATOS DEL LIBRO PAG 288 1 BAR (COMPRENSIBILIDAD)
Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
Comprensibilidad = [0.50881e-9 0.49167e-9 0.4779e-9 0.46691e-9 0.45825e-9 0.45157e-9 0.4466e-9 0.44314e-9 0.441e-9 0.44005e-9 0.4402e-9 0.44341e-9 0.45014e-9 0.4601e-9 0.47311e-9];
%% m^2/N
%% DATOS DEL LIBRO PAG 296 (\emph{v} Viscocidad cinemática)
Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
v = [1.7920 1.5182 1.3063 1.1386 1.0034 0.89266 0.80070 0.72344 0.65785 0.60166 0.55313 0.474 0.41273 0.36433 0.32547];
v = v*1e-06; %% 
%% m^2/s

%% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (vapor–liquid phase)
Temperatura1 = [273.160 278 284 288 294 298 304 308 314 318 324 336 346 356 366];
Temperatura1 = Temperatura1-(ones*273.15);
Densidad1 = [999.793 999.919 999.575 999.079 997.983 997.042 995.346 994.042  991.848  990.235 987.610  981.671  976.086 969.972 963.363];
%% kg/m^3
%% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (one-single phase)
Temperatura2 = [273.15 275 280 285 295 300 305 310 315 325 335 340 345 355 365];
Temperatura2 = Temperatura2-(ones*273.15);
Densidad2 = [999.843 999.937 999.910  999.517 997.807  996.556 995.075 993.383  991.496  987.187  982.233  979.535  976.699  970.628 964.057];
%% kg/m^3

for i=1:1:length(Comprensibilidad)
MEA(i)=1/Comprensibilidad(i);
Visabs(i)=v(i)*Densidad1(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Values are Interpolated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Viscosidad_av = interp1(Temperatura,  v,Temp_av, 'spline','extrap'); % Interpolación si es necesario
Comprensibilidad_av = interp1(Temperatura, Comprensibilidad, Temp_av, 'spline','extrap'); % Interpolación si es necesario
Densidad_av = interp1(Temperatura1, Densidad1, Temp_av, 'spline','extrap'); % Interpolación si es necesario
MEA_av = interp1(Temperatura,  MEA,Temp_av, 'spline','extrap'); % Interpolación si es necesario

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Initial Conditions Qin0, Qout0, Hin0, Hout0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q0in_av= mean(QinFiltrado(1:6100,1));
Q0out_av= mean(QoutFiltrado(1:6100,1));
H0in_av= mean(HinFiltrado(1:6100,1));
H0out_av= mean(HoutFiltrado(1:6100,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Rein,Reout for f1 y f2, miu1, miu2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rein=((Q0in_av/A)*D)/(Viscosidad_av);
Reout=((Q0in_av/A)*D)/(Viscosidad_av);


f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;

miu1=f1/(2*D*A);
miu2=f2/(2*D*A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Leq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V1=2*(H0in_av-H0out_av)*D*g*A^2/(f1*Q0in_av^2);
V2=2*(H0in_av-H0out_av)*D*g*A^2/(f2*Q0out_av^2);
Leq=(V1+V2)/2;
z1 = Leq*0.3;
H20=15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Declare again E,ro,D,e,K to obtain b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=E;
ro = Densidad_av;
D = D;
e=e;
K=MEA_av;

b=sqrt((K/ro)/(1+K*D/(e*E)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Observer initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X1=Q0in_av;
X2=17;
X3=Q0out_av;
X4=z1;
X5=0;

X1i=X1;
X2i=X2;
X3i=X3;
X4i=X4;%z1
X5i=X5;%0
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Bounds for each varying parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1max=max(QinFiltrado);
x1min=min(QinFiltrado);

x2max=19.5;
x2min=12;

x3max=max(QoutFiltrado);
x3min=min(QoutFiltrado)-0.0001;

x4max=85;
x4min=15;
tic
Cd=[1 0 0 0 0;0 0 1 0 0];
Dd=0;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LMIs Declaration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P     = sdpvar(5,5);
Y     = sdpvar(5,5,'symmetric');
W1    = sdpvar(2,5);
W2    = sdpvar(2,5);
W3    = sdpvar(2,5);
W4    = sdpvar(2,5);
W5    = sdpvar(2,5);
W6    = sdpvar(2,5);
W7    = sdpvar(2,5);
W8    = sdpvar(2,5);
W9    = sdpvar(2,5);
W10    = sdpvar(2,5);
W11    = sdpvar(2,5);
W12    = sdpvar(2,5);
W13    = sdpvar(2,5);
W14    = sdpvar(2,5);
W15    = sdpvar(2,5);
W16    = sdpvar(2,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Local Linear Models (Vertex Systems)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ar=A;
L=Leq;
A1 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min)
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];

B1=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];

SYS1=ss(A1,B1,Cd,Dd);
SYSd1=c2d(SYS1,dt);
[adis1,bdis1,cdis1,ddis1]=ssdata(SYSd1);

%%           
A2 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min)
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0]; 
B2=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS2=ss(A2,B2,Cd,Dd);
SYSd2=c2d(SYS2,dt);
[adis2,bdis2,cdis2,ddis2]=ssdata(SYSd2);

%%
A3 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B3=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS3=ss(A3,B3,Cd,Dd);
SYSd3=c2d(SYS3,dt);
[adis3,bdis3,cdis3,ddis3]=ssdata(SYSd3);

%%              
A4 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0]; 
              
B4=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS4=ss(A4,B4,Cd,Dd);
SYSd4=c2d(SYS4,dt);
[adis4,bdis4,cdis4,ddis4]=ssdata(SYSd4);
%%
A5 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B5=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS5=ss(A5,B5,Cd,Dd);
SYSd5=c2d(SYS5,dt);
[adis5,bdis5,cdis5,ddis5]=ssdata(SYSd5);
%%
A6 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B6=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS6=ss(A6,B6,Cd,Dd);
SYSd6=c2d(SYS6,dt);
[adis6,bdis6,cdis6,ddis6]=ssdata(SYSd6);
%%              
A7 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B7=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS7=ss(A7,B7,Cd,Dd);
SYSd7=c2d(SYS7,dt);
[adis7,bdis7,cdis7,ddis7]=ssdata(SYSd7);
%%              
A8 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B8=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS8=ss(A8,B8,Cd,Dd);
SYSd8=c2d(SYS8,dt);
[adis8,bdis8,cdis8,ddis8]=ssdata(SYSd8);
%%              
A9 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B9=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS9=ss(A9,B9,Cd,Dd);
SYSd9=c2d(SYS9,dt);
[adis9,bdis9,cdis9,ddis9]=ssdata(SYSd9);
%%
A10 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B10=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS10=ss(A10,B10,Cd,Dd);
SYSd10=c2d(SYS10,dt);
[adis10,bdis10,cdis10,ddis10]=ssdata(SYSd10);
%%
A11 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B11=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS11=ss(A11,B11,Cd,Dd);
SYSd11=c2d(SYS11,dt);
[adis11,bdis11,cdis11,ddis11]=ssdata(SYSd11);
%%
A12 = [-miu1*x1max, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2min/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2min)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B12=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS12=ss(A12,B12,Cd,Dd);
SYSd12=c2d(SYS12,dt);
[adis12,bdis12,cdis12,ddis12]=ssdata(SYSd12);
%%
A13 = [-miu1*x1max, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
B13=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS13=ss(A13,B13,Cd,Dd);
SYSd13=c2d(SYS13,dt);
[adis13,bdis13,cdis13,ddis13]=ssdata(SYSd13);
%%
A14 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2min/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2min)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B14=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS14=ss(A14,B14,Cd,Dd);
SYSd14=c2d(SYS14,dt);
[adis14,bdis14,cdis14,ddis14]=ssdata(SYSd14);
%%
A15 = [-miu1*x1min, (-g*Ar/x4max)*0.5, 0, (-g*Ar*x2max/(x4max)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4max), 0, -b^2/(g*Ar*x4max), 0, -b^2*sqrt(x2max)/(g*Ar*x4max) %x2 x4
                  0, g*Ar/(L-x4max), -miu2*x3min, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B15=[g*Ar/x4max, 0
    0, 0
    0, -g*Ar/(L-x4max)
    0, 0
    0, 0];
SYS15=ss(A15,B15,Cd,Dd);
SYSd15=c2d(SYS15,dt);
[adis15,bdis15,cdis15,ddis15]=ssdata(SYSd15);
%%
A16 = [-miu1*x1min, (-g*Ar/x4min)*0.5, 0, (-g*Ar*x2max/(x4min)^2)*0.5, 0 % x1 x4
      b^2/(g*Ar*x4min), 0, -b^2/(g*Ar*x4min), 0, -b^2*sqrt(x2max)/(g*Ar*x4min) %x2 x4
                  0, g*Ar/(L-x4min), -miu2*x3max, 0, 0 %x3 x4
                  0, 0, 0, 0, 0
                  0, 0, 0, 0, 0];
              
B16=[g*Ar/x4min, 0
    0, 0
    0, -g*Ar/(L-x4min)
    0, 0
    0, 0];
SYS16=ss(A16,B16,Cd,Dd);
SYSd16=c2d(SYS16,dt);
[adis16,bdis16,cdis16,ddis16]=ssdata(SYSd16);
              
Ad1 = adis1;
Ad2=adis2;
Ad3=adis3;
Ad4=adis4;
Ad5=adis5;
Ad6=adis6;
Ad7=adis7;
Ad8=adis8;
Ad9=adis9;
Ad10=adis10;
Ad11=adis11;
Ad12=adis12;
Ad13=adis13;
Ad14=adis14;
Ad15=adis15;
Ad16=adis16;
C= [1 0 0 0 0; 0 0 1 0 0]


% %% Alternative
% Qk=[1*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-5 0 0;0 0 0 2000 0;0 0 0 0 1*10^-7];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Covariances Matrices Q and R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative
%Qk=[7.5*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-5 0 0;0 0 0 3000 0;0 0 0 0 1*10^-7];

Qk=[7.5*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-5 0 0;0 0 0 3000 0;0 0 0 0 7.5*10^-8];

Rk=[1*10^-5 0;0 1*10^-5]; %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       LQR dual LMIs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_obs=Rk;
H = Qk^0.5;

sdpvar gamma_LQR_obs;

Trace_cond  = [gamma_LQR_obs*eye(5) eye(5); eye(5) Y];


H1     = [      -Y       ,  Y*Ad1-W1'*C  ,   Y*H'     ,      W1';
                  Ad1'*Y-C'*W1   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W1        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H2     = [      -Y       ,  Y*Ad2-W2'*C  ,   Y*H'     ,      W2';
                  Ad2'*Y-C'*W2   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W2        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H3     = [      -Y       ,  Y*Ad3-W3'*C  ,   Y*H'     ,      W3';
                  Ad3'*Y-C'*W3   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W3        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H4     = [      -Y       ,  Y*Ad4-W4'*C  ,   Y*H'     ,      W4';
                  Ad4'*Y-C'*W4   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W4        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H5     = [      -Y       ,  Y*Ad5-W5'*C  ,   Y*H'     ,      W5';
                  Ad5'*Y-C'*W5   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W5        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H6     = [      -Y       ,  Y*Ad6-W6'*C  ,   Y*H'     ,      W6';
                  Ad6'*Y-C'*W6   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W6        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H7     = [      -Y       ,  Y*Ad7-W7'*C  ,   Y*H'     ,      W7';
                  Ad7'*Y-C'*W7   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W7        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H8     = [      -Y       ,  Y*Ad8-W8'*C  ,   Y*H'     ,      W8';
                  Ad8'*Y-C'*W8   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W8        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H9     = [      -Y       ,  Y*Ad9-W9'*C  ,   Y*H'     ,      W9';
                  Ad9'*Y-C'*W9   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W9        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H10     = [      -Y       ,  Y*Ad10-W10'*C  ,   Y*H'     ,      W10';
                  Ad10'*Y-C'*W10   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W10        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H11     = [      -Y       ,  Y*Ad11-W11'*C  ,   Y*H'     ,      W11';
                  Ad11'*Y-C'*W11   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W11        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H12     = [      -Y       ,  Y*Ad12-W12'*C  ,   Y*H'     ,      W12';
                  Ad12'*Y-C'*W12   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W12        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H13     = [      -Y       ,  Y*Ad13-W13'*C  ,   Y*H'     ,      W13';
                  Ad13'*Y-C'*W13   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W13        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H14     = [      -Y       ,  Y*Ad14-W14'*C  ,   Y*H'     ,      W14';
                  Ad14'*Y-C'*W14   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W14        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H15     = [      -Y       ,  Y*Ad15-W15'*C  ,   Y*H'     ,      W15';
                  Ad15'*Y-C'*W15   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W15        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];
                   
H16     = [      -Y       ,  Y*Ad16-W16'*C  ,   Y*H'     ,      W16';
                  Ad16'*Y-C'*W16   ,     -Y         , zeros(5,5) ,   zeros(5,2);
                       H*Y       ,  zeros(5,5)    ,  -eye(5)   ,   zeros(5,2) ;
                       W16        ,   zeros(2,5)   , zeros(2,5) , -inv(R_obs)];

Constraints = [Y>=0,H1<=0,H2<=0,H3<=0,H4<=0,H5<=0,H6<=0,H7<=0,H8<=0,H9<=0,H10<=0,H11<=0,H12<=0,H13<=0,H14<=0,H15<=0,H16<=0,Trace_cond >= 0];
ops= sdpsettings('solver','mosek','verbose',0);
[opt_info]  = optimize(Constraints,[gamma_LQR_obs],ops);
tic
Y_sol  = value(Y);%  Y solution -> WARNING: LQR solves for -Y!!
        
        W1_sol = value(W1); % W1 solution
        W2_sol = value(W2); % W2 solution
        W3_sol = value(W3); % W3 solution
        W4_sol = value(W4); % W4 solution
        
        W5_sol = value(W5); % W5 solution
        W6_sol = value(W6); % W6 solution
        W7_sol = value(W7); % W7 solution
        W8_sol = value(W8); % W8 solution
        
        W9_sol = value(W9); % W9 solution
        W10_sol = value(W10); % W10 solution
        W11_sol = value(W11); % W11 solution
        W12_sol = value(W12); % W12 solution
        
        W13_sol = value(W13); % W13 solution
        W14_sol = value(W14); % W14 solution
        W15_sol = value(W15); % W15 solution
        W16_sol = value(W16); % W16 solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Poles for each vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        L1 = inv(Y_sol)*W1_sol'
        L2 = inv(Y_sol)*W2_sol'
        L3 = inv(Y_sol)*W3_sol'
        L4 = inv(Y_sol)*W4_sol'
        
        L5 = inv(Y_sol)*W5_sol'
        L6 = inv(Y_sol)*W6_sol'
        L7 = inv(Y_sol)*W7_sol'
        L8 = inv(Y_sol)*W8_sol'
        
        L9 = inv(Y_sol)*W9_sol'
        L10 = inv(Y_sol)*W10_sol'
        L11 = inv(Y_sol)*W11_sol'
        L12 = inv(Y_sol)*W12_sol'
        
        L13 = inv(Y_sol)*W13_sol'
        L14 = inv(Y_sol)*W14_sol'
        L15 = inv(Y_sol)*W15_sol'
        L16 = inv(Y_sol)*W16_sol'

        polos1_obs  = eig(adis1-L1*C)
        polos2_obs  = eig(adis2-L2*C)
        polos3_obs  = eig(adis3-L3*C)
        polos4_obs  = eig(adis4-L4*C)
        polos5_obs  = eig(adis5-L5*C)
        polos6_obs  = eig(adis6-L6*C)
        polos7_obs  = eig(adis7-L7*C)
        polos8_obs  = eig(adis8-L8*C)
        polos9_obs  = eig(adis9-L9*C)
        polos10_obs  = eig(adis10-L10*C)
        polos11_obs  = eig(adis11-L11*C)
        polos12_obs  = eig(adis12-L12*C)
        polos13_obs  = eig(adis13-L13*C)
        polos14_obs  = eig(adis14-L14*C)
        polos15_obs  = eig(adis15-L15*C)
        polos16_obs  = eig(adis16-L16*C)

        n=1
enclave=0;
errorcuadraticoz1 = 0;
MAE = 0;
tamano = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   OBSERVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (n <= muestras)
    %Inputs and Outputs
    Qin = QinFiltrado(n, 1);
    Qout = QoutFiltrado(n, 1);
    Hin = HinFiltrado(n, 1);
    Hout = HoutFiltrado(n, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Leak Alarm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Qin - Qout) > 1e-04 && enclave ==0
        enclave = 1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Observer triggered if leak alarm = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if enclave ==1

             %Scheduling Variable Bounds
            Q1max = max(QinFiltrado);
            Q1min = min(QinFiltrado);
            H2max = 19.5;
            H2min = 12;
            Q2max = max(QoutFiltrado);
            Q2min = min(QoutFiltrado)-0.0001;
            z1max = 85;
            z1min = 15;

            % Normalization
            Q1Sch = (Q1max - x1) / (Q1max - Q1min);
            H2Sch = (H2max - x2) / (H2max - H2min);
            Q2Sch = (Q2max - x3) / (Q2max - Q2min);
            z1Sch = (z1max - x4) / (z1max - z1min);
            tic

            if n>=69677
                tic
            end

            % System States
            Q1 = x1;
            H2 = x2;
            Q2 = x3;
            z1 = x4;
            lam = x5;
            Hk = [1 0 0 0 0; 0 0 1 0 0];

            % COMPUTE Fluid and pipeline parameters
            A = Ar;
            promediolongitud = Leq;
            Temp_av = TemperaturaFiltrada(n);
            Viscosidad_av = interp1(Temperatura,  v,Temp_av, 'spline','extrap'); % Interpolación si es necesario
            Comprensibilidad_av = interp1(Temperatura, Comprensibilidad, Temp_av, 'spline','extrap'); % Interpolación si es necesario
            Densidad_av = interp1(Temperatura1, Densidad1, Temp_av, 'spline','extrap'); % Interpolación si es necesario
            MEA_av = interp1(Temperatura,  MEA,Temp_av, 'spline','extrap'); % Interpolación si es necesario
            E=E;
            ro = Densidad_av;
            D = D;
            e=e;
            K=MEA_av;
            Rein=((Qin/A)*D)/(Viscosidad_av);
            Reout=((Qout/A)*D)/(Viscosidad_av);

            f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
            f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;
            miu1=f1/(2*D*A);
            miu2=f2/(2*D*A);

            b=sqrt((K/ro)/(1+K*D/(e*E)));

            %%  PREDICTED STATE VARIABLES (Heuns)
            Q11=(-g*A*(H2-Hin)/(z1)-f1*Q1^2/(2*D*A))*dt+Q1;
            H21=(-b^2*(Q2-Q1+lam*sqrt(H2))/(g*A*z1))*dt+H2; 
            Q21=(-g*A*(Hout-H2)/(promediolongitud-z1)-f2*Q2^2/(2*D*A))*dt+Q2;
            tic
            %% 1- Predicted LPV (Heuns)
            Aact1 = [-miu1 * Q11, (-g * A / z1) * 0.5, 0, (-g * A * H21 / z1^2) * 0.5, 0;
                     b^2 / (g * A * z1), 0, -b^2 / (g * A * z1), 0, (-b^2 * sqrt(abs(H21))) / (g * A * z1);
                     0, g * A / (promediolongitud - z1), -miu2 * Q21, 0, 0;
                     0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0];
            xact1 = [Q11, H21, Q21, z1, lam]' * (dt / 2);
            Bact1 = [(g * A) / z1, 0;
                     0, 0;
                     0, (-g * A) / (promediolongitud - z1);
                     0, 0;
                     0, 0];
            u = [Hin * (dt / 2); Hout * (dt / 2)];
            Ax1 = Aact1 * xact1;
            Bu1 = Bact1 * u;

            %% 2 - Current LPV (Heuns)
            Aact2 = [-miu1 * Q1, (-g * A / z1) * 0.5, 0, (-g * A * H2 / z1^2) * 0.5, 0;
                     b^2 / (g * A * z1), 0, -b^2 / (g * A * z1), 0, (-b^2 * sqrt(abs(H2))) / (g * A * z1);
                     0, g * A / (promediolongitud - z1), -miu2 * Q2, 0, 0;
                     0, 0, 0, 0, 0;
                     0, 0, 0, 0, 0];
            xact2 = [Q1, H2, Q2, z1, lam]' * (dt / 2);
            Bact2 = [(g * A) / z1, 0;
                     0, 0;
                     0, (-g * A) / (promediolongitud - z1);
                     0, 0;
                     0, 0];
            Ax2 = Aact2 * xact2;
            Bu2 = Bact2 * u;

            %% States
            xx = [Q1; H2; Q2; z1; lam];

            %% System without gain
            dXmat = (Ax1 + Ax2) + (Bu1 + Bu2) + xx;
            dQ1 = dXmat(1);
            dH2 = dXmat(2);
            dQ2 = dXmat(3);
            dXm = [dQ1; dH2; dQ2; z1; lam];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Estimated Scheduling Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_1(n)= Q1Sch * H2Sch * Q2Sch * z1Sch; % 0 0 0 0
mu_2(n)= (1-Q1Sch) * H2Sch * Q2Sch * z1Sch; % 1 0 0 0
mu_3(n)= Q1Sch * (1-H2Sch) * Q2Sch * z1Sch; % 0 1 0 0
mu_4(n)= Q1Sch * H2Sch * (1-Q2Sch) * z1Sch; % 0 0 1 0

mu_5(n)= Q1Sch * H2Sch * Q2Sch * (1-z1Sch); % 0 0 0 1
mu_6(n)= (1-Q1Sch) * (1-H2Sch) * (1-Q2Sch) * (1-z1Sch); % 1 1 1 1
mu_7(n)= Q1Sch * (1-H2Sch) * (1-Q2Sch) * (1-z1Sch); % 0 1 1 1 
mu_8(n)= (1-Q1Sch) * H2Sch * (1-Q2Sch) * (1-z1Sch); % 1 0 1 1

mu_9(n)= (1-Q1Sch) * (1-H2Sch) * Q2Sch * (1-z1Sch); % 1 1 0 1
mu_10(n)= (1-Q1Sch) * (1-H2Sch) * (1-Q2Sch) * z1Sch; % 1 1 1 0
mu_11(n)= (1-Q1Sch) * (1-H2Sch) * Q2Sch * z1Sch; % 1 1 0 0
mu_12(n)= (1-Q1Sch) * H2Sch * (1-Q2Sch) * z1Sch; % 1 0 1 0
 
mu_13(n)= (1-Q1Sch) * H2Sch * Q2Sch * (1-z1Sch); % 1 0 0 1
mu_14(n)= Q1Sch * H2Sch * (1-Q2Sch) * (1-z1Sch); % 0 0 1 1
mu_15(n)= Q1Sch * (1-H2Sch) * Q2Sch * (1-z1Sch); % 0 1 0 1
mu_16(n)= Q1Sch * (1-H2Sch) * (1-Q2Sch) * z1Sch; % 0 1 1 0           


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Convex Property
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
musL1(n)=mu_1(n)+mu_2(n)+mu_3(n)+mu_4(n)+mu_5(n)+mu_6(n)+mu_7(n)+mu_8(n)+mu_9(n)+mu_10(n)+mu_11(n)+mu_12(n)+mu_13(n)+mu_14(n)+mu_15(n)+mu_16(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Interpolated Gain and Observer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 L_interp            = mu_1(n)*L1 + mu_2(n)*L2 + mu_3(n)*L3 + mu_4(n)*L4...
                       + mu_5(n)*L5 + mu_6(n)*L6 + mu_7(n)*L7 + mu_8(n)*L8...
                       + mu_9(n)*L9 + mu_10(n)*L10 + mu_11(n)*L11 + mu_12(n)*L12...
                       + mu_13(n)*L13 + mu_14(n)*L14 + mu_15(n)*L15 + mu_16(n)*L16;
                   
          Kk=L_interp;
            dX = dXm + Kk * ([Qin; Qout] - [dQ1; dQ2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   System updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x1 = dX(1);
            x2 = dX(2);
            x3 = dX(3);
            x4 = dX(4);
            x5 = dX(5);

    else
        %% NO  LEAK CASE
        %% INITIAL CONDITIONS REMAINS
        x1 = X1i;
        x2 = X2i;
        x3 = X3i;
        x4 = X4i;
        x5 = X5i;

        % Inicialización de vectores
        dX = [x1; x2; x3; x4; x5];
        mus = 0;
        mu_vals = zeros(16, 1);  % Inicializamos todos los `mu_i` a 0
        tic
    end

    Q1Fuga1Observador(n)=dX(1);
    H2Fuga1Observador(n)=dX(2);
    Q2Fuga1Observador(n)=dX(3);
    z1Fuga1Observador(n)=dX(4);
    lamda1Fuga1Observador(n)=dX(5);
            errorz1(n) = abs(17.04-(z1Fuga1Observador(n)*(Lr/Leq)));
            if enclave == 1 
            errorcuadraticoz1 = errorcuadraticoz1 + (17.04-(z1Fuga1Observador(n)*(Lr/Leq)))^2;
            MAE = MAE + abs(17.04-(z1Fuga1Observador(n)*(Lr/Leq)));
            end
        n=n+1;
end
%% Estimation start time
RMSEz1 = sqrt(errorcuadraticoz1 / (70001-10121));
MAEz1 = 1/(70001-10121) * MAE;
tic
Ql = lamda1Fuga1Observador.* sqrt(H2Fuga1Observador)


Ts = 0.01; % tiempo de muestreo
t = 0:Ts:(length(QinFiltrado)-1)*Ts;



%% GRAFICO LA EVOLUCION DE MU'S
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','MUS EVOLUTION    ','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 700]); 
ylim(axes1,[-0.01 0.4]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Genera colores RGB
colors = lines(16);  % Genera 64 colores distintos
% Plotea cada mu_Leak2 en un color diferente usando un ciclo for
colors = lines(16); % 16 colores distintos

for l = 1:16
    mu_name = sprintf('mu_%d', l);       % genera el nombre: 'mu_1', 'mu_2', ...
    y = eval(mu_name);                   % obtiene los datos de mu_l
    plot(t, y, 'Color', colors(l,:), 'LineWidth', 0.5); hold on;
end
grid on;
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
%print(gcf, fullfile(output_folder, 'musEvolutionL1.eps'), '-depsc');




%% GRÁFICA DE MUS, LEAK 1, LEAK2, LEAK3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','SUMAMUS','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Configuración de los límites y ticks del eje X
xlim(axes1,[0 700]); 
ylim(axes1,[-0.1 1.1]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');



% Colores en tonos de gris
plot(t, musL1, 'color', [0, 0.45, 0.74], 'linewidth', 1.5), hold on  % Gris oscuro
% grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$\sum \psi_{Leak_1}(t)$';
leyenda2 = '$\sum \psi_{Leak_2}(t)$';



%% GRÁFICA DE HIN Y HOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Hin&Hout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 700]); 
ylim(axes1,[8 22]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,HinFiltrado, 'color', [0, 0.45, 0.74],'linewidth',0.5), hold on
plot(t,HoutFiltrado,'color', [0.85, 0.33, 0.10] ,'linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$H_{in}(t)$';
leyenda2 = '$H_{out}(t)$';



%% GRÁFICA DE QIN Y QOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Qin&Qout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 700]); 
ylim(axes1,[0.0088 0.0095]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,QinFiltrado, 'color', [0, 0.45, 0.74],'linewidth',0.5), hold on
plot(t,QoutFiltrado,'color', [0.85, 0.33, 0.10],'linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^3/s)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$Q_{in}(t)$';
leyenda2 = '$Q_{out}(t)$';


%% GRÁFICA DE QinFiltrado-QoutFiltrado
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','FLOWSDIFFERENCE','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 700]); 
ylim(axes1,[-0.0002 0.00065]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
plot(t, QinFiltrado-QoutFiltrado,'color', [0, 0.45, 0.74], 'linewidth', 0.5), hold on
plot(t,Q1Fuga1Observador-Q2Fuga1Observador,'--','color',[0.85, 0.33, 0.10],'linewidth',1), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^{3}/s)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$|Q_{in}-Q_{out}|$';



%% GRÁFICA DE LAM1,LAM2,LAM3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','magnitudleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);

Lam1Ext = [lamda1Fuga1Observador];

%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 700]); 
ylim(axes1,[0 1.5e-04]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
% Línea azul punteada
plot(t, Lam1Ext, 'Color', [0, 0.45, 0.74], 'LineWidth', 0.5), hold on

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^{5/2}/s)','FontName','Times New Roman','FontSize', 20);






%% GRÁFICA DE z1,z2,z3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 700]); 
ylim(axes1,[10 25]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
z1Real = ones(1,70000)*17;
z1Obs = [z1Fuga1Observador]*(Lr/Leq);
z1Obs_filtrado = movmean(z1Obs, 2000, 'Endpoints', 'shrink');  % este es el más común


plot(t, z1Real, '--','color',[0.2 0.2 0.2], 'linewidth', 2), hold on
plot(t, z1Obs,'color', [0,0.5,1], 'linewidth', 1), hold on
plot(t, z1Obs_filtrado,'color',  [0.85, 0.65, 0.13], 'linewidth', 1.5), hold on

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$z_{1}(t)$';
leyenda2 =  '$\hat{z}_{1}(t)$';
leyenda3 =  '$z_{2}(t)$';
leyenda4 =  '$\hat{z}_{2}(t)$';

z1prom = sum(z1Fuga1Observador(1,(10100:end)))/length(z1Fuga1Observador(1,(10100:end)));

tic
vectorexpandz1 = ones(1,60000)*z1prom;
z1expand = [z1Fuga1Observador,vectorexpandz1];
z1transformado = z1expand * (Lr/Leq);

tic
tamano = length(lamda1Fuga1Observador);
for i =1 :1 : tamano
    Ql(i)=lamda1Fuga1Observador(i)*sqrt(H2Fuga1Observador(i));
end


lambda1prom = sum(lamda1Fuga1Observador(1,(12000:end)))/length(lamda1Fuga1Observador(1,(12000:end)));



%%-------------------------------------------------------------------------%%
%%                 FIRST LEAK DIAGNOSING ENDS HERE                                        
%%-------------------------------------------------------------------------%%
%%                 REPEAT THE METHODOLOGY FOR THE
%%                        SECOND LEAK
%%-------------------------------------------------------------------------%%


%% Z1 IS IN LEQ


z1=z1prom;
lam1 = lambda1prom;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       IMPORT ALL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = Datos.data(:,1)*1e-03;
QinFiltrado = Datos.data(:,2)*1e-03;
HinFiltrado = Datos.data(:,5);
HoutFiltrado = Datos.data(:,6);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          DATA FOR
%%                       SECOND LEAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = QoutFiltrado((1:130100),1);
QinFiltrado = QinFiltrado((1:130100),1);
HoutFiltrado = HoutFiltrado((1:130100),1);
HinFiltrado = HinFiltrado((1:130100),1);

iteraciones = length(QinFiltrado);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Saco Temperatura Promedio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TemperaturaFiltrada = Datos.data(:, 9);
Temp_av = sum(TemperaturaFiltrada((1:130100),1))/length(TemperaturaFiltrada(1:130100));
% Parametro = ((Temp_av-55)/30.277);
% 
% YoungsModulus = (34.247*Parametro.^5 + 94.875*Parametro.^4 - 323.75*Parametro.^3 + 799.98*Parametro.^2 - 2549.7*Parametro + 3368.6)*g*1e4;
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Hago Tablas Para Interpolar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATOS DEL LIBRO PAG 288 1 BAR (COMPRENSIBILIDAD)
Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
Comprensibilidad = [0.50881e-9 0.49167e-9 0.4779e-9 0.46691e-9 0.45825e-9 0.45157e-9 0.4466e-9 0.44314e-9 0.441e-9 0.44005e-9 0.4402e-9 0.44341e-9 0.45014e-9 0.4601e-9 0.47311e-9];
%% m^2/N
%% DATOS DEL LIBRO PAG 296 (\emph{v} Viscocidad cinemática)
Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
v = [1.7920 1.5182 1.3063 1.1386 1.0034 0.89266 0.80070 0.72344 0.65785 0.60166 0.55313 0.474 0.41273 0.36433 0.32547];
v = v*1e-06; %% 
%% m^2/s

%% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (vapor–liquid phase)
Temperatura1 = [273.160 278 284 288 294 298 304 308 314 318 324 336 346 356 366];
Temperatura1 = Temperatura1-(ones*273.15);
Densidad1 = [999.793 999.919 999.575 999.079 997.983 997.042 995.346 994.042  991.848  990.235 987.610  981.671  976.086 969.972 963.363];
%% kg/m^3
%% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (one-single phase)
Temperatura2 = [273.15 275 280 285 295 300 305 310 315 325 335 340 345 355 365];
Temperatura2 = Temperatura2-(ones*273.15);
Densidad2 = [999.843 999.937 999.910  999.517 997.807  996.556 995.075 993.383  991.496  987.187  982.233  979.535  976.699  970.628 964.057];
%% kg/m^3

for i=1:1:length(Comprensibilidad)
MEA(i)=1/Comprensibilidad(i);
Visabs(i)=v(i)*Densidad1(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Interpolo los valores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Viscosidad_av = interp1(Temperatura,  v,Temp_av, 'spline','extrap'); % Interpolación si es necesario
Comprensibilidad_av = interp1(Temperatura, Comprensibilidad, Temp_av, 'spline','extrap'); % Interpolación si es necesario
Densidad_av = interp1(Temperatura, Densidad1, Temp_av, 'spline','extrap'); % Interpolación si es necesario
MEA_av = interp1(Temperatura,  MEA,Temp_av, 'spline','extrap'); % Interpolación si es necesario



z0=Leq*0.40;
x1 = Q1Fuga1Observador(end);
x2 = 17;
x3 = Q2Fuga1Observador(end);
x4 = 15;
x5 = Q2Fuga1Observador(end);
x6 = z0;
x7 = 0;

x1i = x1;
x2i = x2;
x3i = x3;
x4i = x4;
x5i = x5;
x6i = x6;
x7i = x7;

tic
%% METO ECUACIONES PARA LMI
P       = sdpvar(7,7);
Y       = sdpvar(7,7,'symmetric');
W1  = sdpvar(2,7);
W2  = sdpvar(2,7);
W3  = sdpvar(2,7);
W4  = sdpvar(2,7);
W5  = sdpvar(2,7);
W6  = sdpvar(2,7);
W7  = sdpvar(2,7);
W8  = sdpvar(2,7);
W9  = sdpvar(2,7);
W10 = sdpvar(2,7);
W11 = sdpvar(2,7);
W12 = sdpvar(2,7);
W13 = sdpvar(2,7);
W14 = sdpvar(2,7);
W15 = sdpvar(2,7);
W16 = sdpvar(2,7);
W17 = sdpvar(2,7);
W18 = sdpvar(2,7);
W19 = sdpvar(2,7);
W20 = sdpvar(2,7);
W21 = sdpvar(2,7);
W22 = sdpvar(2,7);
W23 = sdpvar(2,7);
W24 = sdpvar(2,7);
W25 = sdpvar(2,7);
W26 = sdpvar(2,7);
W27 = sdpvar(2,7);
W28 = sdpvar(2,7);
W29 = sdpvar(2,7);
W30 = sdpvar(2,7);
W31 = sdpvar(2,7);
W32 = sdpvar(2,7);
W33 = sdpvar(2,7);
W34 = sdpvar(2,7);
W35 = sdpvar(2,7);
W36 = sdpvar(2,7);
W37 = sdpvar(2,7);
W38 = sdpvar(2,7);
W39 = sdpvar(2,7);
W40 = sdpvar(2,7);
W41 = sdpvar(2,7);
W42 = sdpvar(2,7);
W43 = sdpvar(2,7);
W44 = sdpvar(2,7);
W45 = sdpvar(2,7);
W46 = sdpvar(2,7);
W47 = sdpvar(2,7);
W48 = sdpvar(2,7);
W49 = sdpvar(2,7);
W50 = sdpvar(2,7);
W51 = sdpvar(2,7);
W52 = sdpvar(2,7);
W53 = sdpvar(2,7);
W54 = sdpvar(2,7);
W55 = sdpvar(2,7);
W56 = sdpvar(2,7);
W57 = sdpvar(2,7);
W58 = sdpvar(2,7);
W59 = sdpvar(2,7);
W60 = sdpvar(2,7);
W61 = sdpvar(2,7);
W62 = sdpvar(2,7);
W63 = sdpvar(2,7);
W64 = sdpvar(2,7);

x1max=max(QinFiltrado)+0.0001;
x1min=min(QinFiltrado)-0.0001;

x2max=18.5;
x2min=16;


x3max=x1max-(lam1*sqrt(x2max));
x3min= x1min-(lam1*sqrt(x2min));

x4max=17.5;
x4min=13;

x5max = max(QoutFiltrado)+0.0001;
x5min = min(QoutFiltrado)-0.0001;

x6max = 90;
x6min = 35;

%%90,39 es el mínimo. Pero las estimaciones salen con mucho ruido
%%80,50 es lo mejor

% x3max = (x1max+x5max)/2;
% x3min = (x1min+x5min)/2;
tic
Cd=[1 0 0 0 0 0 0;0 0 0 0 1 0 0];
Dd=0;
Combinatoria = matrizVerdad(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Saco Rein,Reout para f1 y f2, miu1, miu2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qintermediate = Q0in_av-(lam1*sqrt(H2Fuga1Observador(end)));
Qintermediate = (Q0in_av+Q0out_av)/2;

Rein=((Q0in_av/A)*D)/(Viscosidad_av);
Reinter = ((Qintermediate/A)*D)/(Viscosidad_av)
Reout=((Q0in_av/A)*D)/(Viscosidad_av);


f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;
f3=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reinter) ) ).^2;
miu1=f1/(2*D*A);
miu2=f2/(2*D*A);
miu3=f3/(2*D*A);

b=sqrt((K/ro)/(1+K*D/(e*E)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Leq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V1=2*(H0in_av-H0out_av)*D*g*A^2/(f1*Q0in_av^2);
V2=2*(H0in_av-H0out_av)*D*g*A^2/(f2*Q0out_av^2);
V3=2*(H0in_av-H0out_av)*D*g*A^2/(f3*Qintermediate^2);
Leq=(V1+V2+V3)/3;

[adis,bdis,cdis,ddis,matricesA,matricesB] = evaluarF2(D,Leq,b,f1,f2,f3,z1,lam1,x1max,x1min,x2max,x2min,x3max,x3min,x4max,x4min,x5max,x5min,x6max,x6min,Combinatoria,6,Cd,Dd,dt);


Ad1 = adis{1,1}; 
Ad2 = adis{2,1}; 
Ad3 = adis{3,1}; 
Ad4 = adis{4,1}; 
Ad5 = adis{5,1}; 
Ad6 = adis{6,1}; 
Ad7 = adis{7,1}; 
Ad8 = adis{8,1}; 
Ad9 = adis{9,1}; 
Ad10 = adis{10,1}; 
Ad11 = adis{11,1}; 
Ad12 = adis{12,1}; 
Ad13 = adis{13,1}; 
Ad14 = adis{14,1}; 
Ad15 = adis{15,1}; 
Ad16 = adis{16,1}; 
Ad17 = adis{17,1}; 
Ad18 = adis{18,1}; 
Ad19 = adis{19,1}; 
Ad20 = adis{20,1}; 
Ad21 = adis{21,1}; 
Ad22 = adis{22,1}; 
Ad23 = adis{23,1}; 
Ad24 = adis{24,1}; 
Ad25 = adis{25,1}; 
Ad26 = adis{26,1}; 
Ad27 = adis{27,1}; 
Ad28 = adis{28,1}; 
Ad29 = adis{29,1}; 
Ad30 = adis{30,1}; 
Ad31 = adis{31,1}; 
Ad32 = adis{32,1}; 
Ad33 = adis{33,1}; 
Ad34 = adis{34,1}; 
Ad35 = adis{35,1}; 
Ad36 = adis{36,1}; 
Ad37 = adis{37,1}; 
Ad38 = adis{38,1}; 
Ad39 = adis{39,1}; 
Ad40 = adis{40,1}; 
Ad41 = adis{41,1}; 
Ad42 = adis{42,1}; 
Ad43 = adis{43,1}; 
Ad44 = adis{44,1}; 
Ad45 = adis{45,1}; 
Ad46 = adis{46,1}; 
Ad47 = adis{47,1}; 
Ad48 = adis{48,1}; 
Ad49 = adis{49,1}; 
Ad50 = adis{50,1}; 
Ad51 = adis{51,1}; 
Ad52 = adis{52,1}; 
Ad53 = adis{53,1}; 
Ad54 = adis{54,1}; 
Ad55 = adis{55,1}; 
Ad56 = adis{56,1}; 
Ad57 = adis{57,1}; 
Ad58 = adis{58,1}; 
Ad59 = adis{59,1}; 
Ad60 = adis{60,1}; 
Ad61 = adis{61,1}; 
Ad62 = adis{62,1}; 
Ad63 = adis{63,1}; 
Ad64 = adis{64,1};
C=Cd;

%Qk=[1*10^-8 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1*10^-8 0 0 0 0;0 0 0 2*10^-6 0 0 0;0 0 0 0 1*10^-10 0 0; 0 0 0 0 0 8*10^3 0; 0 0 0 0 0 0 5*10^-12];
%% QK ORIGINAL
Qk=[1*10^-5 0 0 0 0 0 0;
    0 1*10^-2 0 0 0 0 0;
    0 0 1*10^-5 0 0 0 0;
    0 0 0 1*10^-2 0 0 0;
    0 0 0 0 1*10^-5 0 0;
    0 0 0 0 0 2000 0;
    0 0 0 0 0 0 1*10^-7];

%% MATRIX PAPER
Qk=[1*10^-5 0 0 0 0 0 0;
    0 0.01 0 0 0 0 0;
    0 0 1*10^-6 0 0 0 0;
    0 0 0 1*10^2 0 0 0;
    0 0 0 0 1*10^-5 0 0;
    0 0 0 0 0 4000 0;
    0 0 0 0 0 0 1*10^-7];


% Qk=[1*10^-5 0 0 0 0 0 0;
%     0 1*10^-2 0 0 0 0 0;
%     0 0 1*10^-5 0 0 0 0;
%     0 0 0 1*10^-2 0 0 0;
%     0 0 0 0 1*10^-5 0 0;
%     0 0 0 0 0 8000 0;
%     0 0 0 0 0 0 1*10^-8];

%% PROPUESTA
Qk=[1*10^-6 0 0 0 0 0 0;
    0 1*10^-2 0 0 0 0 0;
    0 0 1*10^-6 0 0 0 0;
    0 0 0 1*10^-2 0 0 0;
    0 0 0 0 1*10^-6 0 0;
    0 0 0 0 0 5000 0;
    0 0 0 0 0 0 1*10^-7];

Rk=[1*10^-5 0;0 1*10^-5]; %funciona bien para muestro de 100 Hz

%% EXPERIMENTAL
%% FUNCIONA BIEN, ESTÁ BUENA
% Qk=[1*10^-8 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1*10^-8 0 0 0 0;0 0 0 2*10^-6 0 0 0;0 0 0 0 1*10^-10 0 0; 0 0 0 0 0 8*10^5 0; 0 0 0 0 0 0 5*10^-12];


%Qk=[1*10^-8 0 0 0 0 0 0;0 0.001 0 0 0 0 0;0 0 1*10^-16 0 0 0 0;0 0 0 2*10^-12 0 0 0;0 0 0 0 1*10^-24 0 0; 0 0 0 0 0 8*10^5 0; 0 0 0 0 0 0 5*10^-36];


R_obs=Rk;
H = Qk.^0.5;
tic
sdpvar gamma_LQR_obs;

Trace_cond  = [gamma_LQR_obs*eye(7) eye(7); eye(7) Y];

H1  = [ -Y , Y*Ad1-W1'*C , Y*H' , W1';
        Ad1'*Y-C'*W1 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W1 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H2  = [ -Y , Y*Ad2-W2'*C , Y*H' , W2';
        Ad2'*Y-C'*W2 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W2 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H3  = [ -Y , Y*Ad3-W3'*C , Y*H' , W3';
        Ad3'*Y-C'*W3 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W3 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H4  = [ -Y , Y*Ad4-W4'*C , Y*H' , W4';
        Ad4'*Y-C'*W4 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W4 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H5  = [ -Y , Y*Ad5-W5'*C , Y*H' , W5';
        Ad5'*Y-C'*W5 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W5 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H6  = [ -Y , Y*Ad6-W6'*C , Y*H' , W6';
        Ad6'*Y-C'*W6 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W6 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H7  = [ -Y , Y*Ad7-W7'*C , Y*H' , W7';
        Ad7'*Y-C'*W7 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W7 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H8  = [ -Y , Y*Ad8-W8'*C , Y*H' , W8';
        Ad8'*Y-C'*W8 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W8 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H9  = [ -Y , Y*Ad9-W9'*C , Y*H' , W9';
        Ad9'*Y-C'*W9 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W9 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H10 = [ -Y , Y*Ad10-W10'*C , Y*H' , W10';
        Ad10'*Y-C'*W10 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W10 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H11 = [ -Y , Y*Ad11-W11'*C , Y*H' , W11';
        Ad11'*Y-C'*W11 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W11 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H12 = [ -Y , Y*Ad12-W12'*C , Y*H' , W12';
        Ad12'*Y-C'*W12 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W12 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H13 = [ -Y , Y*Ad13-W13'*C , Y*H' , W13';
        Ad13'*Y-C'*W13 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W13 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H14 = [ -Y , Y*Ad14-W14'*C , Y*H' , W14';
        Ad14'*Y-C'*W14 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W14 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H15 = [ -Y , Y*Ad15-W15'*C , Y*H' , W15';
        Ad15'*Y-C'*W15 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W15 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H16 = [ -Y , Y*Ad16-W16'*C , Y*H' , W16';
        Ad16'*Y-C'*W16 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W16 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H17 = [ -Y , Y*Ad17-W17'*C , Y*H' , W17';
        Ad17'*Y-C'*W17 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W17 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H18 = [ -Y , Y*Ad18-W18'*C , Y*H' , W18';
        Ad18'*Y-C'*W18 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W18 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H19 = [ -Y , Y*Ad19-W19'*C , Y*H' , W19';
        Ad19'*Y-C'*W19 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W19 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H20 = [ -Y , Y*Ad20-W20'*C , Y*H' , W20';
        Ad20'*Y-C'*W20 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W20 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H21 = [ -Y , Y*Ad21-W21'*C , Y*H' , W21';
        Ad21'*Y-C'*W21 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W21 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H22 = [ -Y , Y*Ad22-W22'*C , Y*H' , W22';
        Ad22'*Y-C'*W22 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W22 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H23 = [ -Y , Y*Ad23-W23'*C , Y*H' , W23';
        Ad23'*Y-C'*W23 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W23 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H24 = [ -Y , Y*Ad24-W24'*C , Y*H' , W24';
        Ad24'*Y-C'*W24 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W24 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H25 = [ -Y , Y*Ad25-W25'*C , Y*H' , W25';
        Ad25'*Y-C'*W25 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W25 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H26 = [ -Y , Y*Ad26-W26'*C , Y*H' , W26';
        Ad26'*Y-C'*W26 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W26 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H27 = [ -Y , Y*Ad27-W27'*C , Y*H' , W27';
        Ad27'*Y-C'*W27 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W27 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H28 = [ -Y , Y*Ad28-W28'*C , Y*H' , W28';
        Ad28'*Y-C'*W28 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W28 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H29 = [ -Y , Y*Ad29-W29'*C , Y*H' , W29';
        Ad29'*Y-C'*W29 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W29 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H30 = [ -Y , Y*Ad30-W30'*C , Y*H' , W30';
        Ad30'*Y-C'*W30 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W30 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H31 = [ -Y , Y*Ad31-W31'*C , Y*H' , W31';
        Ad31'*Y-C'*W31 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W31 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H32 = [ -Y , Y*Ad32-W32'*C , Y*H' , W32';
        Ad32'*Y-C'*W32 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W32 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H33 = [ -Y , Y*Ad33-W33'*C , Y*H' , W33';
        Ad33'*Y-C'*W33 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W33 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H34 = [ -Y , Y*Ad34-W34'*C , Y*H' , W34';
        Ad34'*Y-C'*W34 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W34 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H35 = [ -Y , Y*Ad35-W35'*C , Y*H' , W35';
        Ad35'*Y-C'*W35 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W35 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H36 = [ -Y , Y*Ad36-W36'*C , Y*H' , W36';
        Ad36'*Y-C'*W36 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W36 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H37 = [ -Y , Y*Ad37-W37'*C , Y*H' , W37';
        Ad37'*Y-C'*W37 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W37 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H38 = [ -Y , Y*Ad38-W38'*C , Y*H' , W38';
        Ad38'*Y-C'*W38 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W38 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H39 = [ -Y , Y*Ad39-W39'*C , Y*H' , W39';
        Ad39'*Y-C'*W39 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W39 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H40 = [ -Y , Y*Ad40-W40'*C , Y*H' , W40';
        Ad40'*Y-C'*W40 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W40 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];
H41 = [ -Y , Y*Ad41-W41'*C , Y*H' , W41';
        Ad41'*Y-C'*W41 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W41 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H42 = [ -Y , Y*Ad42-W42'*C , Y*H' , W42';
        Ad42'*Y-C'*W42 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W42 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H43 = [ -Y , Y*Ad43-W43'*C , Y*H' , W43';
        Ad43'*Y-C'*W43 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W43 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H44 = [ -Y , Y*Ad44-W44'*C , Y*H' , W44';
        Ad44'*Y-C'*W44 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W44 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H45 = [ -Y , Y*Ad45-W45'*C , Y*H' , W45';
        Ad45'*Y-C'*W45 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W45 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H46 = [ -Y , Y*Ad46-W46'*C , Y*H' , W46';
        Ad46'*Y-C'*W46 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W46 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H47 = [ -Y , Y*Ad47-W47'*C , Y*H' , W47';
        Ad47'*Y-C'*W47 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W47 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H48 = [ -Y , Y*Ad48-W48'*C , Y*H' , W48';
        Ad48'*Y-C'*W48 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W48 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H49 = [ -Y , Y*Ad49-W49'*C , Y*H' , W49';
        Ad49'*Y-C'*W49 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W49 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H50 = [ -Y , Y*Ad50-W50'*C , Y*H' , W50';
        Ad50'*Y-C'*W50 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W50 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H51 = [ -Y , Y*Ad51-W51'*C , Y*H' , W51';
        Ad51'*Y-C'*W51 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W51 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H52 = [ -Y , Y*Ad52-W52'*C , Y*H' , W52';
        Ad52'*Y-C'*W52 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W52 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H53 = [ -Y , Y*Ad53-W53'*C , Y*H' , W53';
        Ad53'*Y-C'*W53 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W53 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H54 = [ -Y , Y*Ad54-W54'*C , Y*H' , W54';
        Ad54'*Y-C'*W54 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W54 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H55 = [ -Y , Y*Ad55-W55'*C , Y*H' , W55';
        Ad55'*Y-C'*W55 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W55 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H56 = [ -Y , Y*Ad56-W56'*C , Y*H' , W56';
        Ad56'*Y-C'*W56 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W56 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H57 = [ -Y , Y*Ad57-W57'*C , Y*H' , W57';
        Ad57'*Y-C'*W57 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W57 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H58 = [ -Y , Y*Ad58-W58'*C , Y*H' , W58';
        Ad58'*Y-C'*W58 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W58 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H59 = [ -Y , Y*Ad59-W59'*C , Y*H' , W59';
        Ad59'*Y-C'*W59 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W59 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H60 = [ -Y , Y*Ad60-W60'*C , Y*H' , W60';
        Ad60'*Y-C'*W60 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W60 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H61 = [ -Y , Y*Ad61-W61'*C , Y*H' , W61';
        Ad61'*Y-C'*W61 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W61 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H62 = [ -Y , Y*Ad62-W62'*C , Y*H' , W62';
        Ad62'*Y-C'*W62 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W62 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H63 = [ -Y , Y*Ad63-W63'*C , Y*H' , W63';
        Ad63'*Y-C'*W63 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W63 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];

H64 = [ -Y , Y*Ad64-W64'*C , Y*H' , W64';
        Ad64'*Y-C'*W64 , -Y , zeros(7,7) , zeros(7,2);
        H*Y , zeros(7,7) , -eye(7) , zeros(7,2);
        W64 , zeros(2,7) , zeros(2,7) , -inv(R_obs)];



tic
ops= sdpsettings('solver','mosek','verbose',0);

[opt_info] = optimize([[H1<=0] + [H2<=0] + [H3<=0] + [H4<=0] + [H5<=0] + [H6<=0] + [H7<=0] + [H8<=0] + [H9<=0] + [H10<=0] + [H11<=0] + [H12<=0] + [H13<=0] + [H14<=0] + [H15<=0] + [H16<=0] + ...
                       [H17<=0] + [H18<=0] + [H19<=0] + [H20<=0] + [H21<=0] + [H22<=0] + [H23<=0] + [H24<=0] + [H25<=0] + [H26<=0] + [H27<=0] + [H28<=0] + [H29<=0] + [H30<=0] + [H31<=0] + [H32<=0] + ...
                       [H33<=0] + [H34<=0] + [H35<=0] + [H36<=0] + [H37<=0] + [H38<=0] + [H39<=0] + [H40<=0] + [H41<=0] + [H42<=0] + [H43<=0] + [H44<=0] + [H45<=0] + [H46<=0] + [H47<=0] + [H48<=0] + ...
                       [H49<=0] + [H50<=0] + [H51<=0] + [H52<=0] + [H53<=0] + [H54<=0] + [H55<=0] + [H56<=0] + [H57<=0] + [H58<=0] + [H59<=0] + [H60<=0] + [H61<=0] + [H62<=0] + [H63<=0] + [H64<=0] + ...
                       [Y>=0] + [Trace_cond>=0] ], [gamma_LQR_obs],ops);
%gamma_LQR_obs<=2e+05
tic
tic
Y_sol  = value(Y);%  Y solution -> WARNING: LQR solves for -Y!!
W1_sol  = value(W1); % W1 solution
W2_sol  = value(W2); % W2 solution
W3_sol  = value(W3); % W3 solution
W4_sol  = value(W4); % W4 solution
W5_sol  = value(W5); % W5 solution
W6_sol  = value(W6); % W6 solution
W7_sol  = value(W7); % W7 solution
W8_sol  = value(W8); % W8 solution
W9_sol  = value(W9); % W9 solution
W10_sol = value(W10); % W10 solution
W11_sol = value(W11); % W11 solution
W12_sol = value(W12); % W12 solution
W13_sol = value(W13); % W13 solution
W14_sol = value(W14); % W14 solution
W15_sol = value(W15); % W15 solution
W16_sol = value(W16); % W16 solution
W17_sol = value(W17); % W17 solution
W18_sol = value(W18); % W18 solution
W19_sol = value(W19); % W19 solution
W20_sol = value(W20); % W20 solution
W21_sol = value(W21); % W21 solution
W22_sol = value(W22); % W22 solution
W23_sol = value(W23); % W23 solution
W24_sol = value(W24); % W24 solution
W25_sol = value(W25); % W25 solution
W26_sol = value(W26); % W26 solution
W27_sol = value(W27); % W27 solution
W28_sol = value(W28); % W28 solution
W29_sol = value(W29); % W29 solution
W30_sol = value(W30); % W30 solution
W31_sol = value(W31); % W31 solution
W32_sol = value(W32); % W32 solution
W33_sol = value(W33); % W33 solution
W34_sol = value(W34); % W34 solution
W35_sol = value(W35); % W35 solution
W36_sol = value(W36); % W36 solution
W37_sol = value(W37); % W37 solution
W38_sol = value(W38); % W38 solution
W39_sol = value(W39); % W39 solution
W40_sol = value(W40); % W40 solution
W41_sol = value(W41); % W41 solution
W42_sol = value(W42); % W42 solution
W43_sol = value(W43); % W43 solution
W44_sol = value(W44); % W44 solution
W45_sol = value(W45); % W45 solution
W46_sol = value(W46); % W46 solution
W47_sol = value(W47); % W47 solution
W48_sol = value(W48); % W48 solution
W49_sol = value(W49); % W49 solution
W50_sol = value(W50); % W50 solution
W51_sol = value(W51); % W51 solution
W52_sol = value(W52); % W52 solution
W53_sol = value(W53); % W53 solution
W54_sol = value(W54); % W54 solution
W55_sol = value(W55); % W55 solution
W56_sol = value(W56); % W56 solution
W57_sol = value(W57); % W57 solution
W58_sol = value(W58); % W58 solution
W59_sol = value(W59); % W59 solution
W60_sol = value(W60); % W60 solution
W61_sol = value(W61); % W61 solution
W62_sol = value(W62); % W62 solution
W63_sol = value(W63); % W63 solution
W64_sol = value(W64); % W64 solution



L1  = inv(Y_sol)*W1_sol.';
L2  = inv(Y_sol)*W2_sol.';
L3  = inv(Y_sol)*W3_sol.';
L4  = inv(Y_sol)*W4_sol.';
L5  = inv(Y_sol)*W5_sol.';
L6  = inv(Y_sol)*W6_sol.';
L7  = inv(Y_sol)*W7_sol.';
L8  = inv(Y_sol)*W8_sol.';
L9  = inv(Y_sol)*W9_sol.';
L10 = inv(Y_sol)*W10_sol.';
L11 = inv(Y_sol)*W11_sol.';
L12 = inv(Y_sol)*W12_sol.';
L13 = inv(Y_sol)*W13_sol.';
L14 = inv(Y_sol)*W14_sol.';
L15 = inv(Y_sol)*W15_sol.';
L16 = inv(Y_sol)*W16_sol.';
L17 = inv(Y_sol)*W17_sol.';
L18 = inv(Y_sol)*W18_sol.';
L19 = inv(Y_sol)*W19_sol.';
L20 = inv(Y_sol)*W20_sol.';
L21 = inv(Y_sol)*W21_sol.';
L22 = inv(Y_sol)*W22_sol.';
L23 = inv(Y_sol)*W23_sol.';
L24 = inv(Y_sol)*W24_sol.';
L25 = inv(Y_sol)*W25_sol.';
L26 = inv(Y_sol)*W26_sol.';
L27 = inv(Y_sol)*W27_sol.';
L28 = inv(Y_sol)*W28_sol.';
L29 = inv(Y_sol)*W29_sol.';
L30 = inv(Y_sol)*W30_sol.';
L31 = inv(Y_sol)*W31_sol.';
L32 = inv(Y_sol)*W32_sol.';
L33 = inv(Y_sol)*W33_sol.';
L34 = inv(Y_sol)*W34_sol.';
L35 = inv(Y_sol)*W35_sol.';
L36 = inv(Y_sol)*W36_sol.';
L37 = inv(Y_sol)*W37_sol.';
L38 = inv(Y_sol)*W38_sol.';
L39 = inv(Y_sol)*W39_sol.';
L40 = inv(Y_sol)*W40_sol.';
L41 = inv(Y_sol)*W41_sol.';
L42 = inv(Y_sol)*W42_sol.';
L43 = inv(Y_sol)*W43_sol.';
L44 = inv(Y_sol)*W44_sol.';
L45 = inv(Y_sol)*W45_sol.';
L46 = inv(Y_sol)*W46_sol.';
L47 = inv(Y_sol)*W47_sol.';
L48 = inv(Y_sol)*W48_sol.';
L49 = inv(Y_sol)*W49_sol.';
L50 = inv(Y_sol)*W50_sol.';
L51 = inv(Y_sol)*W51_sol.';
L52 = inv(Y_sol)*W52_sol.';
L53 = inv(Y_sol)*W53_sol.';
L54 = inv(Y_sol)*W54_sol.';
L55 = inv(Y_sol)*W55_sol.';
L56 = inv(Y_sol)*W56_sol.';
L57 = inv(Y_sol)*W57_sol.';
L58 = inv(Y_sol)*W58_sol.';
L59 = inv(Y_sol)*W59_sol.';
L60 = inv(Y_sol)*W60_sol.';
L61 = inv(Y_sol)*W61_sol.';
L62 = inv(Y_sol)*W62_sol.';
L63 = inv(Y_sol)*W63_sol.';
L64 = inv(Y_sol)*W64_sol.';



polos1_obs  = eig(Ad1 - L1*C);
polos2_obs  = eig(Ad2 - L2*C);
polos3_obs  = eig(Ad3 - L3*C);
polos4_obs  = eig(Ad4 - L4*C);
polos5_obs  = eig(Ad5 - L5*C);
polos6_obs  = eig(Ad6 - L6*C);
polos7_obs  = eig(Ad7 - L7*C);
polos8_obs  = eig(Ad8 - L8*C);
polos9_obs  = eig(Ad9 - L9*C);
polos10_obs = eig(Ad10 - L10*C);
polos11_obs = eig(Ad11 - L11*C);
polos12_obs = eig(Ad12 - L12*C);
polos13_obs = eig(Ad13 - L13*C);
polos14_obs = eig(Ad14 - L14*C);
polos15_obs = eig(Ad15 - L15*C);
polos16_obs = eig(Ad16 - L16*C);
polos17_obs  = eig(Ad17 - L17*C);
polos18_obs  = eig(Ad18 - L18*C);
polos19_obs  = eig(Ad19 - L19*C);
polos20_obs  = eig(Ad20 - L20*C);
polos21_obs  = eig(Ad21 - L21*C);
polos22_obs  = eig(Ad22 - L22*C);
polos23_obs  = eig(Ad23 - L23*C);
polos24_obs  = eig(Ad24 - L24*C);
polos25_obs  = eig(Ad25 - L25*C);
polos26_obs  = eig(Ad26 - L26*C);
polos27_obs  = eig(Ad27 - L27*C);
polos28_obs  = eig(Ad28 - L28*C);
polos29_obs  = eig(Ad29 - L29*C);
polos30_obs  = eig(Ad30 - L30*C);
polos31_obs  = eig(Ad31 - L31*C);
polos32_obs  = eig(Ad32 - L32*C);
polos33_obs  = eig(Ad33 - L33*C);
polos34_obs  = eig(Ad34 - L34*C);
polos35_obs  = eig(Ad35 - L35*C);
polos36_obs  = eig(Ad36 - L36*C);
polos37_obs  = eig(Ad37 - L37*C);
polos38_obs  = eig(Ad38 - L38*C);
polos39_obs  = eig(Ad39 - L39*C);
polos40_obs  = eig(Ad40 - L40*C);
polos41_obs  = eig(Ad41 - L41*C);
polos42_obs  = eig(Ad42 - L42*C);
polos43_obs  = eig(Ad43 - L43*C);
polos44_obs  = eig(Ad44 - L44*C);
polos45_obs  = eig(Ad45 - L45*C);
polos46_obs  = eig(Ad46 - L46*C);
polos47_obs  = eig(Ad47 - L47*C);
polos48_obs  = eig(Ad48 - L48*C);
polos49_obs  = eig(Ad49 - L49*C);
polos50_obs  = eig(Ad50 - L50*C);
polos51_obs  = eig(Ad51 - L51*C);
polos52_obs  = eig(Ad52 - L52*C);
polos53_obs  = eig(Ad53 - L53*C);
polos54_obs  = eig(Ad54 - L54*C);
polos55_obs  = eig(Ad55 - L55*C);
polos56_obs  = eig(Ad56 - L56*C);
polos57_obs  = eig(Ad57 - L57*C);
polos58_obs  = eig(Ad58 - L58*C);
polos59_obs  = eig(Ad59 - L59*C);
polos60_obs  = eig(Ad60 - L60*C);
polos61_obs  = eig(Ad61 - L61*C);
polos62_obs  = eig(Ad62 - L62*C);
polos63_obs  = eig(Ad63 - L63*C);
polos64_obs  = eig(Ad64 - L64*C);

tic

% Vector de valores propios (polos)
% Graficar el círculo unitario
theta = linspace(0, 2*pi, 100);
x_circulo = cos(theta); % Parte real
y_circulo = sin(theta); % Parte imaginaria

figure;
hold on;
plot(x_circulo, y_circulo, 'r--', 'LineWidth', 2); % Círculo unitario
for i = 1:64
    poloeval = eval(['polos' num2str(i) '_obs'])
    plot(real(poloeval), ...
         imag(poloeval), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Valores propios
end
hold off

% Configuración del gráfico
axis equal;
xlabel('Parte real');
ylabel('Parte imaginaria');
title('Valores propios y el círculo unitario');
legend('Círculo unitario', 'Valores propios');
grid on;

tic
n=1
enclave=0;

            errorcuadraticoz2 = 0;
MAE = 0;
while(n<=iteraciones)
tic
    Qin = QinFiltrado(n,1);
    Qout = QoutFiltrado(n,1);
    Hin = HinFiltrado(n,1);
    Hout = HoutFiltrado(n,1);

    if (Qin-Qout>=6e-04)  && enclave ==0 && n>=66000
        enclave =1;
        tic
    end

    if enclave ==1
        Q1max=max(QinFiltrado)+0.0001;
        Q1min=min(QinFiltrado)-0.0001;

        H2max=18.5;
        H2min=16;


        Q2max=Q1max-(lam1*sqrt(H2max));
        Q2min= Q1min-(lam1*sqrt(H2min));


        H3max=17.5;
        H3min=13;

        Q3max = max(QoutFiltrado)+0.0001;
        Q3min = min(QoutFiltrado)-0.0001;


        % Q2max=(Q1max+Q3max)/2;
        % Q2min= ((Q1min+Q3min)/2);

        z2max = 90;
        z2min = 35;

        Q1Sch = (Q1max - x1) / (Q1max - Q1min);
        H2Sch = (H2max - x2) / (H2max - H2min);
        Q2Sch = (Q2max - x3) / (Q2max - Q2min);
        H3Sch = (H3max - x4) / (H3max - H3min);
        Q3Sch = (Q3max - x5) / (Q3max - Q3min);
        z2Sch = (z2max - x6) / (z2max - z2min);

        if Q1Sch < 0 || H2Sch <0 || Q2Sch <0 || H3Sch<0 || Q3Sch <0 || z2Sch <0
        tic
        end
        if n >=67870
            Q1Sch;
            H2Sch;
            Q2Sch;
            H3Sch;
            Q3Sch;
            z2Sch;
            tic
        end
        Q1 = x1;
        H2 = x2;
        Q2 = x3;
        H3 = x4;
        Q3 = x5;
        z2 = x6;
        lam = x7;
        Hk = [1 0 0 0 0 0 0; 0 0 0 0 1 0 0];

        A = Ar;
        promediolongitud = Leq;

            Temp_av = TemperaturaFiltrada(n);
            Viscosidad_av = interp1(Temperatura,  v,Temp_av, 'spline','extrap'); % Interpolación si es necesario
            Comprensibilidad_av = interp1(Temperatura, Comprensibilidad, Temp_av, 'spline','extrap'); % Interpolación si es necesario
            Densidad_av = interp1(Temperatura1, Densidad1, Temp_av, 'spline','extrap'); % Interpolación si es necesario
            MEA_av = interp1(Temperatura,  MEA,Temp_av, 'spline','extrap'); % Interpolación si es necesario
            E=E;
            ro = Densidad_av;
            D = D;
            e=e;
            K=MEA_av;
            Rein=((Qin/A)*D)/(Viscosidad_av);
            Reout=((Qout/A)*D)/(Viscosidad_av);
% % 
% % 
% % 
% % % f1=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Rein)))^-1)^2;
 f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
% % % f2=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Reout)))^-1)^2;
 f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;
  f3=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;

miu1=f1/(2*D*A);
miu2=f2/(2*D*A);
miu3=f3/(2*D*A);
% Varvis1(n)=miu1;
% Varvis2(n)=miu2;
 %E=(34.247*((Temp_av-55)/30.277)^5 + 94.874*((Temp_av-55)/30.277)^4 - 323.75*((Temp_av-55)/30.277)^3 + 799.98*((Temp_av-55)/30.277)^2 - 2549.7*((Temp_av-55)/30.277) + 3368.6)*g*1e4;

            b=sqrt((K/ro)/(1+K*D/(e*E)));
        dz1 = z1;
        dz2 = (z2-z1);
        dz3 = (Leq-z2);
        %% HAGO MIS Q ESTIMADAS (Qi gorro j+1)
        Q11=Q1+dt*((((-g*A)*(H2-Hin))/(dz1))-miu1*Q1^2);
        Q21=Q2+dt*((((-g*A)*(H3-H2))/(dz2))-miu2*Q2^2);
        Q31=Q3+dt*((((-g*A)*(Hout-H3))/(dz3))-miu3*Q3^2);

        %% HAGO MIS H ESTIMADAS (Hi gorro j+1)
        H21=H2+dt*((-b^2)*(Q2-Q1+lam1*sqrt((abs(H2))))/(g*A*dz1));
        H31=H3+dt*((-b^2)*(Q3-Q2+lam*sqrt((abs(H3))))/(g*A*dz2));

        phi1 = -miu1*Q11;
        phi2 = (g*Ar*(H21-1))/z1;
        phi3 = (-g*Ar*H21^2)/(z1*z2);
        phi4 = (b^2/(g*Ar*z1))*(1-((lam1*sqrt(abs(H21)))/(2*Q11)));
        phi5 = (-b^2/(g*Ar*z1))*(1+((lam1*sqrt(abs(H21)))/(2*Q21)));
        phi6 = (g*Ar)/(z2-z1);
        phi7 = -miu2*Q21;
        phi8 = (-g*Ar)/(z2-z1);
        phi9 = (b^2)/(g*Ar*(z2-z1));
        phi10 = (-b^2)/(g*Ar*(z2-z1));
        phi11 = (-b^2*sqrt(abs(H31)))/(g*Ar*(z2-z1));
        phi12 = (g*Ar)/(L-z2);
        phi13 = -miu3*Q31;

        Aact1 = [phi1 ,phi2, 0, 0 ,0 ,phi3, 0;
        phi4 ,0, phi5, 0, 0, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        0, 0, phi9 ,0 ,phi10, 0, phi11;
        0, 0, 0, phi12, phi13, 0, 0;
        zeros(2,7)];

        xact1 = [Q11, H21, Q21, H31, Q31, z2, lam]' * (dt / 2);
        Bact1=[(g*Ar)/z1, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z2);
        0, 0;
        0, 0];
        u = [Hin * (dt / 2); Hout * (dt / 2)];
        Ax1 = Aact1 * xact1;
        Bu1 = Bact1 * u;

        phi1 = -miu1*Q1;
        phi2 = (g*Ar*(H2-1))/z1;
        phi3 = (-g*Ar*H2^2)/(z1*z2);
        phi4 = (b^2/(g*Ar*z1))*(1-((lam1*sqrt(abs(H2)))/(2*Q1)));
        phi5 = (-b^2/(g*Ar*z1))*(1+((lam1*sqrt(abs(H2)))/(2*Q2)));
        phi6 = (g*Ar)/(z2-z1);
        phi7 = -miu2*Q2;
        phi8 = (-g*Ar)/(z2-z1);
        phi9 = (b^2)/(g*Ar*(z2-z1));
        phi10 = (-b^2)/(g*Ar*(z2-z1));
        phi11 = (-b^2*sqrt(abs(H3)))/(g*Ar*(z2-z1));
        phi12 = (g*Ar)/(L-z2);
        phi13 = -miu3*Q3;

        Aact2 = [phi1 ,phi2, 0, 0 ,0 ,phi3, 0;
        phi4 ,0, phi5, 0, 0, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        0, 0, phi9 ,0 ,phi10, 0, phi11;
        0, 0, 0, phi12, phi13, 0, 0;
        zeros(2,7)];

        xact2 = [Q1, H2, Q2, H3, Q3, z2, lam]' * (dt / 2);
        Bact2=[(g*Ar)/z1, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z2);
        0, 0;
        0, 0];
        u = [Hin * (dt / 2); Hout * (dt / 2)];
        Ax2 = Aact2 * xact2;
        Bu2 = Bact2 * u;

        xx = [Q1; H2; Q2; H3; Q3; z2; lam];
        dXmat = (Ax1 + Ax2) + (Bu1 + Bu2) + xx;
        dQ1 = dXmat(1);
        dH2 = dXmat(2);
        dQ2 = dXmat(3);
        dH3 = dXmat(4);
        dQ3 = dXmat(5);
        dXm = [dQ1; dH2; dQ2; dH3; dQ3; z2; lam];

        %% Cálculo de los coeficientes de pertenencia `mu`
        Q1SchMax= (1 - Q1Sch);
        Q1SchMin= Q1Sch;

        H2SchMax= (1 - H2Sch);
        H2SchMin= H2Sch;

        Q2SchMax= (1 - Q2Sch);
        Q2SchMin= Q2Sch;

        H3SchMax= (1 - H3Sch);
        H3SchMin= H3Sch;

        Q3SchMax= (1 - Q3Sch);
        Q3SchMin= Q3Sch;

        z2SchMax= (1-z2Sch);
        z2SchMin = z2Sch;

        mu_vals = obtenermusF2(Q1SchMax, Q1SchMin, H2SchMax, H2SchMin, Q2SchMax, Q2SchMin,H3SchMax,H3SchMin,Q3SchMax,Q3SchMin, z2SchMax, z2SchMin, Combinatoria, 6);
        tic

        mu(n) = sum(mu_vals);
        tic
        L_interp=0;
        L_interp = L1 * mu_vals(1) + ...
           L2 * mu_vals(2) + ...
           L3 * mu_vals(3) + ...
           L4 * mu_vals(4) + ...
           L5 * mu_vals(5) + ...
           L6 * mu_vals(6) + ...
           L7 * mu_vals(7) + ...
           L8 * mu_vals(8) + ...
           L9 * mu_vals(9) + ...
           L10 * mu_vals(10) + ...
           L11 * mu_vals(11) + ...
           L12 * mu_vals(12) + ...
           L13 * mu_vals(13) + ...
           L14 * mu_vals(14) + ...
           L15 * mu_vals(15) + ...
           L16 * mu_vals(16) + ...
           L17 * mu_vals(17) + ...
           L18 * mu_vals(18) + ...
           L19 * mu_vals(19) + ...
           L20 * mu_vals(20) + ...
           L21 * mu_vals(21) + ...
           L22 * mu_vals(22) + ...
           L23 * mu_vals(23) + ...
           L24 * mu_vals(24) + ...
           L25 * mu_vals(25) + ...
           L26 * mu_vals(26) + ...
           L27 * mu_vals(27) + ...
           L28 * mu_vals(28) + ...
           L29 * mu_vals(29) + ...
           L30 * mu_vals(30) + ...
           L31 * mu_vals(31) + ...
           L32 * mu_vals(32) + ...
           L33 * mu_vals(33) + ...
           L34 * mu_vals(34) + ...
           L35 * mu_vals(35) + ...
           L36 * mu_vals(36) + ...
           L37 * mu_vals(37) + ...
           L38 * mu_vals(38) + ...
           L39 * mu_vals(39) + ...
           L40 * mu_vals(40) + ...
           L41 * mu_vals(41) + ...
           L42 * mu_vals(42) + ...
           L43 * mu_vals(43) + ...
           L44 * mu_vals(44) + ...
           L45 * mu_vals(45) + ...
           L46 * mu_vals(46) + ...
           L47 * mu_vals(47) + ...
           L48 * mu_vals(48) + ...
           L49 * mu_vals(49) + ...
           L50 * mu_vals(50) + ...
           L51 * mu_vals(51) + ...
           L52 * mu_vals(52) + ...
           L53 * mu_vals(53) + ...
           L54 * mu_vals(54) + ...
           L55 * mu_vals(55) + ...
           L56 * mu_vals(56) + ...
           L57 * mu_vals(57) + ...
           L58 * mu_vals(58) + ...
           L59 * mu_vals(59) + ...
           L60 * mu_vals(60) + ...
           L61 * mu_vals(61) + ...
           L62 * mu_vals(62) + ...
           L63 * mu_vals(63) + ...
           L64 * mu_vals(64);

        Kk = L_interp;
        tic
        dX = dXm + Kk * ([Qin; Qout] - [dQ1; dQ3]);
        x1 = dX(1);
        x2 = dX(2);
        x3 = dX(3);
        x4 = dX(4);
        x5 = dX(5);
        x6 = dX(6);
        x7 = dX(7);

        tic
    else
        x1 = x1i;
        x2 = x2i;
        x3 = x3i;
        x4 = x4i;
        x5 = x5i;
        x6 = x6i;
        x7 = x7i;

        % Inicialización de vectores
        dX = [x1; x2; x3; x4; x5; x6; x7];
        mus = 0;
        mu_vals = zeros(64, 1);  % Inicializamos todos los `mu_i` a 0
    end


for l = 1:64
    mu_Leak2{l}(n) = mu_vals(l); % Almacena el valor
end

    tic

    Q1Fuga2Observador(n)=dX(1);
    H2Fuga2Observador(n)=dX(2);
    Q2Fuga2Observador(n)=dX(3);
    H3Fuga2Observador(n)=dX(4);
    Q3Fuga2Observador(n)=dX(5);
    z2Fuga2Observador(n)=dX(6);
    lamda2Fuga2Observador(n)=dX(7);
            errorz2(n) = abs(33.5-(z2Fuga2Observador(n)*(Lr/Leq)));
            if enclave == 1
            errorcuadraticoz2 = errorcuadraticoz2 + (33.47-(z2Fuga2Observador(n)*(Lr/Leq)))^2;
            MAE = MAE + abs(33.47-(z2Fuga2Observador(n)*(Lr/Leq)));

            end
       n=n+1;
end
tic
RMSEz2 = sqrt(errorcuadraticoz2 / (130100-70218));
MAEz2 = 1/(130100-70218) * MAE;


tic

tic
Ts = 0.01; % tiempo de muestreo
t = 0:Ts:(length(QinFiltrado)-1)*Ts;



QoutF1 = Q2Fuga1Observador((1:70000));
QoutF2 = Q3Fuga2Observador(70001:end);

QoutTotalObs = [QoutF1, QoutF2];

QinF1 = Q1Fuga1Observador((1:70000));
QinF2 = Q1Fuga2Observador(70001:end);

QinTotalObs = [QinF1, QinF2];






tic
% %% GRAFICO LA EVOLUCION DE MU'S
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
% figure('Name','MUS EVOLUTION    ','NumberTitle','off','color', [1 1 1])
% axes1 = axes('FontSize',20);
% %% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[650 1300]); 
% ylim(axes1,[-0.01 0.22]); 
% box(axes1,'on');
% grid(axes1,'on');
% hold(axes1,'all');
% % Genera colores RGB
% colors = lines(64);  % Genera 64 colores distintos
% % Plotea cada mu_Leak2 en un color diferente usando un ciclo for
% for l = 1:64
%     % Plotea cada mu_Leak2 con su color correspondiente
%     plot(t, mu_Leak2{l}, 'Color', colors(l, :)); hold on;
%     tic
% end
% grid on;
% xlabel('(s)','FontName','Times New Roman','FontSize', 20);
% print(gcf, fullfile(output_folder, 'musEvolution.eps'), '-depsc');
% tic

%% GRAFICO LA EVOLUCION DE MU'S (64 curvas con tamaño fijo)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

figure('Name','MUS EVOLUTION','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);

% Límites de ejes
xlim(axes1,[650 1300]); 
ylim(axes1,[-0.01 0.23]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Colores
colors = lines(64);

% Graficar cada curva
for l = 1:64
    y = mu_Leak2{l}(:)'; % asegurar vector fila
    if length(y) ~= length(t)
        warning('Curva %d con longitud diferente: len(y)=%d, len(t)=%d', ...
                l, length(y), length(t));
        continue
    end
    plot(t, y, 'Color', colors(l,:), 'LineWidth', 0.5); 
    hold on;
end

% Etiquetas
xlabel('(s)','FontName','Times New Roman','FontSize', 20);

% ==== Ajuste de tamaño fijo ====
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 20, 15]); % [x y ancho alto]
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0, 0, 20, 15]);

% Exportar EPS con mismo tamaño
%print(gcf, fullfile(output_folder, 'musEvolutionL2.eps'), '-depsc','-vector');

% tic
% 
% %% Configuración básica
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
% set(groot, 'defaultLegendInterpreter','latex');
% 
% fig = figure('Name','MUS EVOLUTION','NumberTitle','off','color', [1 1 1]);
% ax = axes(fig,'FontSize',20);
% xlim(ax,[650 1300]); 
% ylim(ax,[-0.01 0.22]); 
% box(ax,'on'); grid(ax,'on'); hold(ax,'all');
% 
% xlabel('(s)','FontName','Times New Roman','FontSize', 20);
% 
% colors = lines(64);

% %% Graficar curvas rasterizadas
% set(fig,'Renderer','opengl'); % rasteriza las curvas
% for l = 1:64
%     y = mu_Leak2{l}(:)'; 
%     plot(ax, t, y, 'Color', colors(l,:), 'LineWidth', 0.5);
% end
% 
% %% Ajuste tamaño
% set(fig, 'Units', 'centimeters', 'Position', [2, 2, 20, 15]);
% set(fig, 'PaperUnits', 'centimeters');
% set(fig, 'PaperPosition', [0, 0, 20, 15]);

% %% Exportar EPS color
% % Las curvas quedan raster dentro del EPS, el resto vectorial
% exportgraphics(ax, fullfile(output_folder, 'musEvolutionL2BUENA.eps'), ...
%     'ContentType','none', ... % rasteriza el contenido (curvas)
%     'BackgroundColor','none', ...
%     'Resolution',800); % resolución ajustable

tic


%% GRÁFICA DE MUS, LEAK 1, LEAK2, LEAK3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','SUMAMUS','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Configuración de los límites y ticks del eje X
xlim(axes1,[0 1300]); 
ylim(axes1,[-0.1 1.1]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

musL1G = [musL1,ones(1,60100)];
muG = [mu];

% Colores en tonos de gris
plot(t, musL1G, 'color', [0, 0.45, 0.74], 'linewidth', 1.5), hold on  % Gris oscuro
plot(t, muG, 'color', [0.85, 0.33, 0.10], 'linewidth', 1.5), hold on  % Gris intermedio
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$\sum \psi_{Leak_1}(t)$';
leyenda2 = '$\sum \psi_{Leak_2}(t)$';
%leyenda = legend(leyenda1,leyenda2,'Location','southeast');
%print(gcf,fullfile(output_folder, 'MUSL2.eps'), '-depsc');

tic


%% GRÁFICA DE HIN Y HOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Hin&Hout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1300]); 
ylim(axes1,[8 22]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,HinFiltrado, 'color', [0, 0.45, 0.74],'linewidth',0.5), hold on
plot(t,HoutFiltrado,'color', [0.85, 0.33, 0.10] ,'linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$H_{in}(t)$';
leyenda2 = '$H_{out}(t)$';
%print(gcf, fullfile(output_folder, 'presionesEntradaYSalidaLPV.eps'), '-depsc');



%% GRÁFICA DE QIN Y QOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Qin&Qout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1300]); 
ylim(axes1,[0.0086 0.0097]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,QinFiltrado, 'color', [0, 0.45, 0.74],'linewidth',0.5), hold on
plot(t,QoutFiltrado,'color', [0.85, 0.33, 0.10],'linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^3/s)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$Q_{in}(t)$';
leyenda2 = '$Q_{out}(t)$';
%print(gcf,fullfile(output_folder, 'FlujoEntradaYSalida.eps'), '-depsc');


%% GRÁFICA DE QinFiltrado-QoutFiltrado
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','FLOWSDIFFERENCE','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1300]); 
ylim(axes1,[-0.0002 0.0010]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
plot(t, QinFiltrado-QoutFiltrado,'color', [0, 0.45, 0.74], 'linewidth', 0.5), hold on
plot(t,QinTotalObs-QoutTotalObs,'--','color',[0.85, 0.33, 0.10],'linewidth',1), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^{3}/s)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$|Q_{in}-Q_{out}|$';
%print(gcf,fullfile(output_folder, 'QRESTA.eps'), '-depsc');



%% GRÁFICA DE LAM1,LAM2,LAM3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','magnitudleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
Ultvalor = sum(lamda1Fuga1Observador(10000:end))/length(lamda1Fuga1Observador(10000:end))
extendido = ones(1,60100)*Ultvalor;
Lam1Ext = [lamda1Fuga1Observador,extendido];

Lam2Ext = [lamda2Fuga2Observador];
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1300]); 
ylim(axes1,[0 1.5e-04]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
% Línea azul punteada
plot(t, Lam1Ext, 'Color', [0, 0.45, 0.74], 'LineWidth', 0.5), hold on

% Línea roja punteada
plot(t, Lam2Ext, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.5), hold on
 
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^{5/2}/s)','FontName','Times New Roman','FontSize', 20);
%print(gcf,fullfile(output_folder, 'lambdasmagnitudes.eps'), '-depsc');







%% GRÁFICA DE z1,z2,z3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1300]); 
ylim(axes1,[10 40]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
z1Real = ones(1,130100)*17;
z2Real = ones(1,130100)*33.5;
z1Obs = mean(z1Fuga1Observador(11000:68000));
exten = ones(1,60100)*z1Obs;
z1Obs = [z1Fuga1Observador,exten]*(Lr/Leq);

z2Obs = [z2Fuga2Observador]*(Lr/Leq);
plot(t, z1Real, '--','color',[0.2 0.2 0.2], 'linewidth', 1.5), hold on
plot(t, z1Obs,'color', [0, 0.45, 0.74], 'linewidth', 1), hold on
plot(t, z2Real, '--', 'color',[0.4 0.4 0.4],'linewidth', 1.5), hold on
plot(t, z2Obs,  'color', [0.85, 0.33, 0.10], 'linewidth', 1), hold on
% plot(t,z2Real,'k--','linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$z_{1}(t)$';
leyenda2 =  '$\hat{z}_{1}(t)$';
leyenda3 =  '$z_{2}(t)$';
leyenda4 =  '$\hat{z}_{2}(t)$';
%print(gcf,fullfile(output_folder, 'leakspositions.eps'), '-depsc');


    tic






















lambda2prom = sum(lamda2Fuga2Observador(1,(70100:end)))/length(lamda2Fuga2Observador(1,(70100:end)));
%z2prom = sum(z2Fuga2Observador(1,(72500:end)))/length(z2Fuga2Observador(1,(72500:end)));
z2prom = sum(z2Fuga2Observador(1,(70100:end)))/length(z2Fuga2Observador(1,(70100:end)));
tic
z2 = z2prom;
lam2= lambda2prom;
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = Datos.data(:,1)*1e-03;
QinFiltrado = Datos.data(:,2)*1e-03;
HinFiltrado = Datos.data(:,5);
HoutFiltrado = Datos.data(:,6);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DATOS 
%%                       THIRD LEAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = QoutFiltrado((1:end),1);
QinFiltrado = QinFiltrado((1:end),1);
HoutFiltrado = HoutFiltrado((1:end),1);
HinFiltrado = HinFiltrado((1:end),1);

iteraciones = length(QinFiltrado);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Saco Temperatura Promedio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TemperaturaFiltrada = Datos.data(:, 9);
Temp_av = sum(TemperaturaFiltrada((1:end),1))/length(TemperaturaFiltrada(1:end));
% Parametro = ((Temp_av-55)/30.277);
% 
% YoungsModulus = (34.247*Parametro.^5 + 94.875*Parametro.^4 - 323.75*Parametro.^3 + 799.98*Parametro.^2 - 2549.7*Parametro + 3368.6)*g*1e4;
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Hago Tablas Para Interpolar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DATOS DEL LIBRO PAG 288 1 BAR (COMPRENSIBILIDAD)
Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
Comprensibilidad = [0.50881e-9 0.49167e-9 0.4779e-9 0.46691e-9 0.45825e-9 0.45157e-9 0.4466e-9 0.44314e-9 0.441e-9 0.44005e-9 0.4402e-9 0.44341e-9 0.45014e-9 0.4601e-9 0.47311e-9];
%% m^2/N
%% DATOS DEL LIBRO PAG 296 (\emph{v} Viscocidad cinemática)
Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
v = [1.7920 1.5182 1.3063 1.1386 1.0034 0.89266 0.80070 0.72344 0.65785 0.60166 0.55313 0.474 0.41273 0.36433 0.32547];
v = v*1e-06; %% 
%% m^2/s

%% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (vapor–liquid phase)
Temperatura1 = [273.160 278 284 288 294 298 304 308 314 318 324 336 346 356 366];
Temperatura1 = Temperatura1-(ones*273.15);
Densidad1 = [999.793 999.919 999.575 999.079 997.983 997.042 995.346 994.042  991.848  990.235 987.610  981.671  976.086 969.972 963.363];
%% kg/m^3
%% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (one-single phase)
Temperatura2 = [273.15 275 280 285 295 300 305 310 315 325 335 340 345 355 365];
Temperatura2 = Temperatura2-(ones*273.15);
Densidad2 = [999.843 999.937 999.910  999.517 997.807  996.556 995.075 993.383  991.496  987.187  982.233  979.535  976.699  970.628 964.057];
%% kg/m^3

for i=1:1:length(Comprensibilidad)
MEA(i)=1/Comprensibilidad(i);
Visabs(i)=v(i)*Densidad1(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Interpolo los valores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Viscosidad_av = interp1(Temperatura,  v,Temp_av, 'spline','extrap'); % Interpolación si es necesario
Comprensibilidad_av = interp1(Temperatura, Comprensibilidad, Temp_av, 'spline','extrap'); % Interpolación si es necesario
Densidad_av = interp1(Temperatura1, Densidad1, Temp_av, 'spline','extrap'); % Interpolación si es necesario
MEA_av = interp1(Temperatura,  MEA,Temp_av, 'spline','extrap'); % Interpolación si es necesario

z0=63;
x1 = Q1Fuga2Observador(end);
x2 = 17;
x3 = Q2Fuga2Observador(end);
x4 = 15;
x5 = Q2Fuga2Observador(end);
x6 = 11;
x7 = Q3Fuga2Observador(end);
x8 = z0;
x9 = 0;
tic

x1i = x1;
x2i = x2;
x3i = x3;
x4i = x4;
x5i = x5;
x6i = x6;
x7i = x7;
x8i = x8;
x9i = x9;

%% METO ECUACIONES PARA LMI
P       = sdpvar(9,9);
Y       = sdpvar(9,9,'symmetric');
W1  = sdpvar(2,9);
W2  = sdpvar(2,9);
W3  = sdpvar(2,9);
W4  = sdpvar(2,9);
W5  = sdpvar(2,9);
W6  = sdpvar(2,9);
W7  = sdpvar(2,9);
W8  = sdpvar(2,9);
W9  = sdpvar(2,9);
W10 = sdpvar(2,9);
W11 = sdpvar(2,9);
W12 = sdpvar(2,9);
W13 = sdpvar(2,9);
W14 = sdpvar(2,9);
W15 = sdpvar(2,9);
W16 = sdpvar(2,9);
W17 = sdpvar(2,9);
W18 = sdpvar(2,9);
W19 = sdpvar(2,9);
W20 = sdpvar(2,9);
W21 = sdpvar(2,9);
W22 = sdpvar(2,9);
W23 = sdpvar(2,9);
W24 = sdpvar(2,9);
W25 = sdpvar(2,9);
W26 = sdpvar(2,9);
W27 = sdpvar(2,9);
W28 = sdpvar(2,9);
W29 = sdpvar(2,9);
W30 = sdpvar(2,9);
W31 = sdpvar(2,9);
W32 = sdpvar(2,9);
W33 = sdpvar(2,9);
W34 = sdpvar(2,9);
W35 = sdpvar(2,9);
W36 = sdpvar(2,9);
W37 = sdpvar(2,9);
W38 = sdpvar(2,9);
W39 = sdpvar(2,9);
W40 = sdpvar(2,9);
W41 = sdpvar(2,9);
W42 = sdpvar(2,9);
W43 = sdpvar(2,9);
W44 = sdpvar(2,9);
W45 = sdpvar(2,9);
W46 = sdpvar(2,9);
W47 = sdpvar(2,9);
W48 = sdpvar(2,9);
W49 = sdpvar(2,9);
W50 = sdpvar(2,9);
W51 = sdpvar(2,9);
W52 = sdpvar(2,9);
W53 = sdpvar(2,9);
W54 = sdpvar(2,9);
W55 = sdpvar(2,9);
W56 = sdpvar(2,9);
W57 = sdpvar(2,9);
W58 = sdpvar(2,9);
W59 = sdpvar(2,9);
W60 = sdpvar(2,9);
W61 = sdpvar(2,9);
W62 = sdpvar(2,9);
W63 = sdpvar(2,9);
W64 = sdpvar(2,9);
W65 = sdpvar(2,9);
W66 = sdpvar(2,9);
W67 = sdpvar(2,9);
W68 = sdpvar(2,9);
W69 = sdpvar(2,9);
W70 = sdpvar(2,9);
W71 = sdpvar(2,9);
W72 = sdpvar(2,9);
W73 = sdpvar(2,9);
W74 = sdpvar(2,9);
W75 = sdpvar(2,9);
W76 = sdpvar(2,9);
W77 = sdpvar(2,9);
W78 = sdpvar(2,9);
W79 = sdpvar(2,9);
W80 = sdpvar(2,9);
W81 = sdpvar(2,9);
W82 = sdpvar(2,9);
W83 = sdpvar(2,9);
W84 = sdpvar(2,9);
W85 = sdpvar(2,9);
W86 = sdpvar(2,9);
W87 = sdpvar(2,9);
W88 = sdpvar(2,9);
W89 = sdpvar(2,9);
W90 = sdpvar(2,9);
W91 = sdpvar(2,9);
W92 = sdpvar(2,9);
W93 = sdpvar(2,9);
W94 = sdpvar(2,9);
W95 = sdpvar(2,9);
W96 = sdpvar(2,9);
W97 = sdpvar(2,9);
W98 = sdpvar(2,9);
W99 = sdpvar(2,9);
W100 = sdpvar(2,9);
W101 = sdpvar(2,9);
W102 = sdpvar(2,9);
W103 = sdpvar(2,9);
W104 = sdpvar(2,9);
W105 = sdpvar(2,9);
W106 = sdpvar(2,9);
W107 = sdpvar(2,9);
W108 = sdpvar(2,9);
W109 = sdpvar(2,9);
W110 = sdpvar(2,9);
W111 = sdpvar(2,9);
W112 = sdpvar(2,9);
W113 = sdpvar(2,9);
W114 = sdpvar(2,9);
W115 = sdpvar(2,9);
W116 = sdpvar(2,9);
W117 = sdpvar(2,9);
W118 = sdpvar(2,9);
W119 = sdpvar(2,9);
W120 = sdpvar(2,9);
W121 = sdpvar(2,9);
W122 = sdpvar(2,9);
W123 = sdpvar(2,9);
W124 = sdpvar(2,9);
W125 = sdpvar(2,9);
W126 = sdpvar(2,9);
W127 = sdpvar(2,9);
W128 = sdpvar(2,9);
W129 = sdpvar(2,9);
W130 = sdpvar(2,9);
W131 = sdpvar(2,9);
W132 = sdpvar(2,9);
W133 = sdpvar(2,9);
W134 = sdpvar(2,9);
W135 = sdpvar(2,9);
W136 = sdpvar(2,9);
W137 = sdpvar(2,9);
W138 = sdpvar(2,9);
W139 = sdpvar(2,9);
W140 = sdpvar(2,9);
W141 = sdpvar(2,9);
W142 = sdpvar(2,9);
W143 = sdpvar(2,9);
W144 = sdpvar(2,9);
W145 = sdpvar(2,9);
W146 = sdpvar(2,9);
W147 = sdpvar(2,9);
W148 = sdpvar(2,9);
W149 = sdpvar(2,9);
W150 = sdpvar(2,9);
W151 = sdpvar(2,9);
W152 = sdpvar(2,9);
W153 = sdpvar(2,9);
W154 = sdpvar(2,9);
W155 = sdpvar(2,9);
W156 = sdpvar(2,9);
W157 = sdpvar(2,9);
W158 = sdpvar(2,9);
W159 = sdpvar(2,9);
W160 = sdpvar(2,9);
W161 = sdpvar(2,9);
W162 = sdpvar(2,9);
W163 = sdpvar(2,9);
W164 = sdpvar(2,9);
W165 = sdpvar(2,9);
W166 = sdpvar(2,9);
W167 = sdpvar(2,9);
W168 = sdpvar(2,9);
W169 = sdpvar(2,9);
W170 = sdpvar(2,9);
W171 = sdpvar(2,9);
W172 = sdpvar(2,9);
W173 = sdpvar(2,9);
W174 = sdpvar(2,9);
W175 = sdpvar(2,9);
W176 = sdpvar(2,9);
W177 = sdpvar(2,9);
W178 = sdpvar(2,9);
W179 = sdpvar(2,9);
W180 = sdpvar(2,9);
W181 = sdpvar(2,9);
W182 = sdpvar(2,9);
W183 = sdpvar(2,9);
W184 = sdpvar(2,9);
W185 = sdpvar(2,9);
W186 = sdpvar(2,9);
W187 = sdpvar(2,9);
W188 = sdpvar(2,9);
W189 = sdpvar(2,9);
W190 = sdpvar(2,9);
W191 = sdpvar(2,9);
W192 = sdpvar(2,9);
W193 = sdpvar(2,9);
W194 = sdpvar(2,9);
W195 = sdpvar(2,9);
W196 = sdpvar(2,9);
W197 = sdpvar(2,9);
W198 = sdpvar(2,9);
W199 = sdpvar(2,9);
W200 = sdpvar(2,9);
W201 = sdpvar(2,9);
W202 = sdpvar(2,9);
W203 = sdpvar(2,9);
W204 = sdpvar(2,9);
W205 = sdpvar(2,9);
W206 = sdpvar(2,9);
W207 = sdpvar(2,9);
W208 = sdpvar(2,9);
W209 = sdpvar(2,9);
W210 = sdpvar(2,9);
W211 = sdpvar(2,9);
W212 = sdpvar(2,9);
W213 = sdpvar(2,9);
W214 = sdpvar(2,9);
W215 = sdpvar(2,9);
W216 = sdpvar(2,9);
W217 = sdpvar(2,9);
W218 = sdpvar(2,9);
W219 = sdpvar(2,9);
W220 = sdpvar(2,9);
W221 = sdpvar(2,9);
W222 = sdpvar(2,9);
W223 = sdpvar(2,9);
W224 = sdpvar(2,9);
W225 = sdpvar(2,9);
W226 = sdpvar(2,9);
W227 = sdpvar(2,9);
W228 = sdpvar(2,9);
W229 = sdpvar(2,9);
W230 = sdpvar(2,9);
W231 = sdpvar(2,9);
W232 = sdpvar(2,9);
W233 = sdpvar(2,9);
W234 = sdpvar(2,9);
W235 = sdpvar(2,9);
W236 = sdpvar(2,9);
W237 = sdpvar(2,9);
W238 = sdpvar(2,9);
W239 = sdpvar(2,9);
W240 = sdpvar(2,9);
W241 = sdpvar(2,9);
W242 = sdpvar(2,9);
W243 = sdpvar(2,9);
W244 = sdpvar(2,9);
W245 = sdpvar(2,9);
W246 = sdpvar(2,9);
W247 = sdpvar(2,9);
W248 = sdpvar(2,9);
W249 = sdpvar(2,9);
W250 = sdpvar(2,9);
W251 = sdpvar(2,9);
W252 = sdpvar(2,9);
W253 = sdpvar(2,9);
W254 = sdpvar(2,9);
W255 = sdpvar(2,9);
W256 = sdpvar(2,9);


%% EN ESTA PRUEBA LO ESTOY CALANDO CON EL QINFILTRADO DE TODO EL EXPEIRMENTO
%% **** EXPERIMENTAR SOLO CON LA MUESTRA DE 20000:END

%% FUNCIONA PARA 2000 Y 20000
x1max=max(QinFiltrado(1:end))+0.0001; %% (Plus Tiny error)
x1min=min(QinFiltrado(1:end))-0.0001; %% (Minus Tiny error)

x2max=19.5;
x2min=11.5;

tic


x3max=x1max-(lam1*sqrt(x2max))+(0.0002); %% (Plus Tiny error)
x3min= x1min-(lam1*sqrt(x2min))-(0.0001); %% (Minus Tiny error)

x4max=16.5;
x4min=10.5;

x5max=x3max-(lam2*sqrt(x4max))+(0.0002); %% (Plus Tiny error)
x5min=x3min-(lam2*sqrt(x4min))-(0.0001); %% (Minus Tiny error)

x6max=15.5;
x6min=9.5;

x7max = max(QoutFiltrado(1:end))+0.0001; %% (Plus Tiny error)
x7min = min(QoutFiltrado(1:end))-0.0001; %% (Minus Tiny error)

%% ESTE PARÁMETRO TIENE QUE HACER UN BARRIDO DE 20 EN 20, SIRVE PARA TENER MÁS PRECISIÓN Y UNA CONVERGENCIA MÁS RÁPIDA.
x8max = 80;
x8min = 60;
tic
Cd=[1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 1 0 0];
Dd=0;
Combinatoria = matrizVerdad(8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Saco Rein,Reout para f1 y f2, miu1, miu2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qintermediate = Q0in_av-(lam1*sqrt(H2Fuga1Observador(end)));
Qintermediate = (Q0in_av+Q0out_av)/2;

Rein=((Q0in_av/A)*D)/(Viscosidad_av);
Reinter = ((Qintermediate/A)*D)/(Viscosidad_av)
Reout=((Q0in_av/A)*D)/(Viscosidad_av);


% f1=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Rein)))^-1)^2;
f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
% f2=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Reout)))^-1)^2;
f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;
f3=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reinter) ) ).^2;
f4=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reinter) ) ).^2;

miu1=f1/(2*D*A);
miu2=f2/(2*D*A);
miu3=f3/(2*D*A);
miu4=f4/(2*D*A);

b=sqrt((K/ro)/(1+K*D/(e*E)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Leq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V1=2*(H0in_av-H0out_av)*D*g*A^2/(f1*Q0in_av^2);
V2=2*(H0in_av-H0out_av)*D*g*A^2/(f2*Q0out_av^2);
V3=2*(H0in_av-H0out_av)*D*g*A^2/(f3*Qintermediate^2);
V4=2*(H0in_av-H0out_av)*D*g*A^2/(f4*Qintermediate^2);
Leq=(V1+V2+V3+V4)/4;
dt=0.01;
Combinatoria = matrizVerdad(8);


%[adis,bdis,cdis,ddis,matricesA,matricesB] = evaluarF3(z1,lam1,z2,lam2,x1max,x1min,x2max,x2min,x3max,x3min,x4max,x4min,x5max,x5min,x6max,x6min,x7max,x7min,x8max,x8min,Combinatoria,8,Cd,Dd,dt);
[adis,bdis,cdis,ddis,matricesA,matricesB] = evaluarF3(D,Leq,b,f1,f2,f3,f4,z1,lam1,z2,lam2,x1max,x1min,x2max,x2min,x3max,x3min,x4max,x4min,x5max,x5min,x6max,x6min,x7max,x7min,x8max,x8min,Combinatoria,8,Cd,Dd,dt);
Ad1 = adis{1,1}; 
Ad2 = adis{2,1}; 
Ad3 = adis{3,1}; 
Ad4 = adis{4,1}; 
Ad5 = adis{5,1}; 
Ad6 = adis{6,1}; 
Ad7 = adis{7,1}; 
Ad8 = adis{8,1}; 
Ad9 = adis{9,1}; 
Ad10 = adis{10,1}; 
Ad11 = adis{11,1}; 
Ad12 = adis{12,1}; 
Ad13 = adis{13,1}; 
Ad14 = adis{14,1}; 
Ad15 = adis{15,1}; 
Ad16 = adis{16,1}; 
Ad17 = adis{17,1}; 
Ad18 = adis{18,1}; 
Ad19 = adis{19,1}; 
Ad20 = adis{20,1}; 
Ad21 = adis{21,1}; 
Ad22 = adis{22,1}; 
Ad23 = adis{23,1}; 
Ad24 = adis{24,1}; 
Ad25 = adis{25,1}; 
Ad26 = adis{26,1}; 
Ad27 = adis{27,1}; 
Ad28 = adis{28,1}; 
Ad29 = adis{29,1}; 
Ad30 = adis{30,1}; 
Ad31 = adis{31,1}; 
Ad32 = adis{32,1}; 
Ad33 = adis{33,1}; 
Ad34 = adis{34,1}; 
Ad35 = adis{35,1}; 
Ad36 = adis{36,1}; 
Ad37 = adis{37,1}; 
Ad38 = adis{38,1}; 
Ad39 = adis{39,1}; 
Ad40 = adis{40,1}; 
Ad41 = adis{41,1}; 
Ad42 = adis{42,1}; 
Ad43 = adis{43,1}; 
Ad44 = adis{44,1}; 
Ad45 = adis{45,1}; 
Ad46 = adis{46,1}; 
Ad47 = adis{47,1}; 
Ad48 = adis{48,1}; 
Ad49 = adis{49,1}; 
Ad50 = adis{50,1}; 
Ad51 = adis{51,1}; 
Ad52 = adis{52,1}; 
Ad53 = adis{53,1}; 
Ad54 = adis{54,1}; 
Ad55 = adis{55,1}; 
Ad56 = adis{56,1}; 
Ad57 = adis{57,1}; 
Ad58 = adis{58,1}; 
Ad59 = adis{59,1}; 
Ad60 = adis{60,1}; 
Ad61 = adis{61,1}; 
Ad62 = adis{62,1};
Ad63 = adis{63,1};
Ad64 = adis{64,1};
Ad65 = adis{65,1};
Ad66 = adis{66,1};
Ad67 = adis{67,1};
Ad68 = adis{68,1};
Ad69 = adis{69,1};
Ad70 = adis{70,1};
Ad71 = adis{71,1};
Ad72 = adis{72,1};
Ad73 = adis{73,1};
Ad74 = adis{74,1};
Ad75 = adis{75,1};
Ad76 = adis{76,1};
Ad77 = adis{77,1};
Ad78 = adis{78,1};
Ad79 = adis{79,1};
Ad80 = adis{80,1};
Ad81 = adis{81,1};
Ad82 = adis{82,1};
Ad83 = adis{83,1};
Ad84 = adis{84,1};
Ad85 = adis{85,1};
Ad86 = adis{86,1};
Ad87 = adis{87,1};
Ad88 = adis{88,1};
Ad89 = adis{89,1};
Ad90 = adis{90,1};
Ad91 = adis{91,1};
Ad92 = adis{92,1};
Ad93 = adis{93,1};
Ad94 = adis{94,1};
Ad95 = adis{95,1};
Ad96 = adis{96,1};
Ad97 = adis{97,1};
Ad98 = adis{98,1};
Ad99 = adis{99,1};
Ad100 = adis{100,1};
Ad101 = adis{101,1};
Ad102 = adis{102,1};
Ad103 = adis{103,1};
Ad104 = adis{104,1};
Ad105 = adis{105,1};
Ad106 = adis{106,1};
Ad107 = adis{107,1};
Ad108 = adis{108,1};
Ad109 = adis{109,1};
Ad110 = adis{110,1};
Ad111 = adis{111,1};
Ad112 = adis{112,1};
Ad113 = adis{113,1};
Ad114 = adis{114,1};
Ad115 = adis{115,1};
Ad116 = adis{116,1};
Ad117 = adis{117,1};
Ad118 = adis{118,1};
Ad119 = adis{119,1};
Ad120 = adis{120,1};
Ad121 = adis{121,1};
Ad122 = adis{122,1};
Ad123 = adis{123,1};
Ad124 = adis{124,1};
Ad125 = adis{125,1};
Ad126 = adis{126,1};
Ad127 = adis{127,1};
Ad128 = adis{128,1};
Ad129 = adis{129,1};
Ad130 = adis{130,1};
Ad131 = adis{131,1};
Ad132 = adis{132,1};
Ad133 = adis{133,1};
Ad134 = adis{134,1};
Ad135 = adis{135,1};
Ad136 = adis{136,1};
Ad137 = adis{137,1};
Ad138 = adis{138,1};
Ad139 = adis{139,1};
Ad140 = adis{140,1};
Ad141 = adis{141,1};
Ad142 = adis{142,1};
Ad143 = adis{143,1};
Ad144 = adis{144,1};
Ad145 = adis{145,1};
Ad146 = adis{146,1};
Ad147 = adis{147,1};
Ad148 = adis{148,1};
Ad149 = adis{149,1};
Ad150 = adis{150,1};
Ad151 = adis{151,1};
Ad152 = adis{152,1};
Ad153 = adis{153,1};
Ad154 = adis{154,1};
Ad155 = adis{155,1};
Ad156 = adis{156,1};
Ad157 = adis{157,1};
Ad158 = adis{158,1};
Ad159 = adis{159,1};
Ad160 = adis{160,1};
Ad161 = adis{161,1};
Ad162 = adis{162,1};
Ad163 = adis{163,1};
Ad164 = adis{164,1};
Ad165 = adis{165,1};
Ad166 = adis{166,1};
Ad167 = adis{167,1};
Ad168 = adis{168,1};
Ad169 = adis{169,1};
Ad170 = adis{170,1};
Ad171 = adis{171,1};
Ad172 = adis{172,1};
Ad173 = adis{173,1};
Ad174 = adis{174,1};
Ad175 = adis{175,1};
Ad176 = adis{176,1};
Ad177 = adis{177,1};
Ad178 = adis{178,1};
Ad179 = adis{179,1};
Ad180 = adis{180,1};
Ad181 = adis{181,1};
Ad182 = adis{182,1};
Ad183 = adis{183,1};
Ad184 = adis{184,1};
Ad185 = adis{185,1};
Ad186 = adis{186,1};
Ad187 = adis{187,1};
Ad188 = adis{188,1};
Ad189 = adis{189,1};
Ad190 = adis{190,1};
Ad191 = adis{191,1};
Ad192 = adis{192,1};
Ad193 = adis{193,1};
Ad194 = adis{194,1};
Ad195 = adis{195,1};
Ad196 = adis{196,1};
Ad197 = adis{197,1};
Ad198 = adis{198,1};
Ad199 = adis{199,1};
Ad200 = adis{200,1};
Ad201 = adis{201,1};
Ad202 = adis{202,1};
Ad203 = adis{203,1};
Ad204 = adis{204,1};
Ad205 = adis{205,1};
Ad206 = adis{206,1};
Ad207 = adis{207,1};
Ad208 = adis{208,1};
Ad209 = adis{209,1};
Ad210 = adis{210,1};
Ad211 = adis{211,1};
Ad212 = adis{212,1};
Ad213 = adis{213,1};
Ad214 = adis{214,1};
Ad215 = adis{215,1};
Ad216 = adis{216,1};
Ad217 = adis{217,1};
Ad218 = adis{218,1};
Ad219 = adis{219,1};
Ad220 = adis{220,1};
Ad221 = adis{221,1};
Ad222 = adis{222,1};
Ad223 = adis{223,1};
Ad224 = adis{224,1};
Ad225 = adis{225,1};
Ad226 = adis{226,1};
Ad227 = adis{227,1};
Ad228 = adis{228,1};
Ad229 = adis{229,1};
Ad230 = adis{230,1};
Ad231 = adis{231,1};
Ad232 = adis{232,1};
Ad233 = adis{233,1};
Ad234 = adis{234,1};
Ad235 = adis{235,1};
Ad236 = adis{236,1};
Ad237 = adis{237,1};
Ad238 = adis{238,1};
Ad239 = adis{239,1};
Ad240 = adis{240,1};
Ad241 = adis{241,1};
Ad242 = adis{242,1};
Ad243 = adis{243,1};
Ad244 = adis{244,1};
Ad245 = adis{245,1};
Ad246 = adis{246,1};
Ad247 = adis{247,1};
Ad248 = adis{248,1};
Ad249 = adis{249,1};
Ad250 = adis{250,1};
Ad251 = adis{251,1};
Ad252 = adis{252,1};
Ad253 = adis{253,1};
Ad254 = adis{254,1};
Ad255 = adis{255,1};
Ad256 = adis{256,1};

C=Cd;

tic
Qk=[1*10^-7 0 0 0 0 0 0 0 0;
    0 1*10^-4 0 0 0 0 0 0 0;
    0 0 1*10^-7 0 0 0 0 0 0;
    0 0 0 1*10^-4 0 0 0 0 0;
    0 0 0 0 1*10^-7 0 0 0 0;
    0 0 0 0 0 1*10^-4 0 0 0;
    0 0 0 0 0 0 1*10^-7 0 0;
    0 0 0 0 0 0 0 10000 0;
    0 0 0 0 0 0 0 0 1*10^-6];

% Qk=[1*10^-7 0 0 0 0 0 0 0 0;
%     0 1*10^-6 0 0 0 0 0 0 0;
%     0 0 1*10^-7 0 0 0 0 0 0;
%     0 0 0 1*10^-6 0 0 0 0 0;
%     0 0 0 0 1*10^-7 0 0 0 0;
%     0 0 0 0 0 1*10^-6 0 0 0;
%     0 0 0 0 0 0 1*10^-7 0 0;
%     0 0 0 0 0 0 0 1000 0;
%     0 0 0 0 0 0 0 0 0.5*10^-7];


Rk=[1*10^-5 0;0 1*10^-5]; %funciona bien para muestro de 100 Hz

R_obs=Rk;
H = Qk.^0.5;
tic
sdpvar gamma_LQR_obs;

Trace_cond  = [gamma_LQR_obs*eye(9) eye(9); eye(9) Y];

% Presupongo que todas las matrices Y, C, H, y R_obs están definidas

% Presupongo que todas las matrices Y, C, H, y R_obs están definidas

H1  = [ -Y , Y*Ad1 - W1'*C , Y*H' , W1';
        Ad1'*Y - C'*W1 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W1 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H2  = [ -Y , Y*Ad2 - W2'*C , Y*H' , W2';
        Ad2'*Y - C'*W2 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W2 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H3  = [ -Y , Y*Ad3 - W3'*C , Y*H' , W3';
        Ad3'*Y - C'*W3 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W3 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H4  = [ -Y , Y*Ad4 - W4'*C , Y*H' , W4';
        Ad4'*Y - C'*W4 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W4 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H5  = [ -Y , Y*Ad5 - W5'*C , Y*H' , W5';
        Ad5'*Y - C'*W5 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W5 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H6  = [ -Y , Y*Ad6 - W6'*C , Y*H' , W6';
        Ad6'*Y - C'*W6 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W6 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H7  = [ -Y , Y*Ad7 - W7'*C , Y*H' , W7';
        Ad7'*Y - C'*W7 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W7 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H8  = [ -Y , Y*Ad8 - W8'*C , Y*H' , W8';
        Ad8'*Y - C'*W8 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W8 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H9  = [ -Y , Y*Ad9 - W9'*C , Y*H' , W9';
        Ad9'*Y - C'*W9 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W9 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H10 = [ -Y , Y*Ad10 - W10'*C , Y*H' , W10';
        Ad10'*Y - C'*W10 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W10 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H11 = [ -Y , Y*Ad11 - W11'*C , Y*H' , W11';
        Ad11'*Y - C'*W11 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W11 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H12 = [ -Y , Y*Ad12 - W12'*C , Y*H' , W12';
        Ad12'*Y - C'*W12 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W12 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H13 = [ -Y , Y*Ad13 - W13'*C , Y*H' , W13';
        Ad13'*Y - C'*W13 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W13 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H14 = [ -Y , Y*Ad14 - W14'*C , Y*H' , W14';
        Ad14'*Y - C'*W14 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W14 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H15 = [ -Y , Y*Ad15 - W15'*C , Y*H' , W15';
        Ad15'*Y - C'*W15 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W15 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H16 = [ -Y , Y*Ad16 - W16'*C , Y*H' , W16';
        Ad16'*Y - C'*W16 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W16 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H17 = [ -Y , Y*Ad17 - W17'*C , Y*H' , W17';
        Ad17'*Y - C'*W17 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W17 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H18 = [ -Y , Y*Ad18 - W18'*C , Y*H' , W18';
        Ad18'*Y - C'*W18 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W18 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H19 = [ -Y , Y*Ad19 - W19'*C , Y*H' , W19';
        Ad19'*Y - C'*W19 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W19 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H20 = [ -Y , Y*Ad20 - W20'*C , Y*H' , W20';
        Ad20'*Y - C'*W20 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W20 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H21 = [ -Y , Y*Ad21 - W21'*C , Y*H' , W21';
        Ad21'*Y - C'*W21 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W21 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H22 = [ -Y , Y*Ad22 - W22'*C , Y*H' , W22';
        Ad22'*Y - C'*W22 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W22 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H23 = [ -Y , Y*Ad23 - W23'*C , Y*H' , W23';
        Ad23'*Y - C'*W23 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W23 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H24 = [ -Y , Y*Ad24 - W24'*C , Y*H' , W24';
        Ad24'*Y - C'*W24 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W24 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H25 = [ -Y , Y*Ad25 - W25'*C , Y*H' , W25';
        Ad25'*Y - C'*W25 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W25 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H26 = [ -Y , Y*Ad26 - W26'*C , Y*H' , W26';
        Ad26'*Y - C'*W26 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W26 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H27 = [ -Y , Y*Ad27 - W27'*C , Y*H' , W27';
        Ad27'*Y - C'*W27 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W27 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H28 = [ -Y , Y*Ad28 - W28'*C , Y*H' , W28';
        Ad28'*Y - C'*W28 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W28 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H29 = [ -Y , Y*Ad29 - W29'*C , Y*H' , W29';
        Ad29'*Y - C'*W29 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W29 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H30 = [ -Y , Y*Ad30 - W30'*C , Y*H' , W30';
        Ad30'*Y - C'*W30 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W30 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H31 = [ -Y , Y*Ad31 - W31'*C , Y*H' , W31';
        Ad31'*Y - C'*W31 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W31 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H32 = [ -Y , Y*Ad32 - W32'*C , Y*H' , W32';
        Ad32'*Y - C'*W32 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W32 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H33 = [ -Y , Y*Ad33 - W33'*C , Y*H' , W33';
        Ad33'*Y - C'*W33 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W33 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H34 = [ -Y , Y*Ad34 - W34'*C , Y*H' , W34';
        Ad34'*Y - C'*W34 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W34 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H35 = [ -Y , Y*Ad35 - W35'*C , Y*H' , W35';
        Ad35'*Y - C'*W35 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W35 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H36 = [ -Y , Y*Ad36 - W36'*C , Y*H' , W36';
        Ad36'*Y - C'*W36 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W36 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H37 = [ -Y , Y*Ad37 - W37'*C , Y*H' , W37';
        Ad37'*Y - C'*W37 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W37 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H38 = [ -Y , Y*Ad38 - W38'*C , Y*H' , W38';
        Ad38'*Y - C'*W38 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W38 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H39 = [ -Y , Y*Ad39 - W39'*C , Y*H' , W39';
        Ad39'*Y - C'*W39 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W39 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H40 = [ -Y , Y*Ad40 - W40'*C , Y*H' , W40';
        Ad40'*Y - C'*W40 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W40 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H41 = [ -Y , Y*Ad41 - W41'*C , Y*H' , W41';
        Ad41'*Y - C'*W41 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W41 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H42 = [ -Y , Y*Ad42 - W42'*C , Y*H' , W42';
        Ad42'*Y - C'*W42 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W42 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H43 = [ -Y , Y*Ad43 - W43'*C , Y*H' , W43';
        Ad43'*Y - C'*W43 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W43 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H44 = [ -Y , Y*Ad44 - W44'*C , Y*H' , W44';
        Ad44'*Y - C'*W44 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W44 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H45 = [ -Y , Y*Ad45 - W45'*C , Y*H' , W45';
        Ad45'*Y - C'*W45 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W45 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H46 = [ -Y , Y*Ad46 - W46'*C , Y*H' , W46';
        Ad46'*Y - C'*W46 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W46 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H47 = [ -Y , Y*Ad47 - W47'*C , Y*H' , W47';
        Ad47'*Y - C'*W47 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W47 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H48 = [ -Y , Y*Ad48 - W48'*C , Y*H' , W48';
        Ad48'*Y - C'*W48 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W48 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H49 = [ -Y , Y*Ad49 - W49'*C , Y*H' , W49';
        Ad49'*Y - C'*W49 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W49 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H50 = [ -Y , Y*Ad50 - W50'*C , Y*H' , W50';
        Ad50'*Y - C'*W50 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W50 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H51 = [ -Y , Y*Ad51 - W51'*C , Y*H' , W51';
        Ad51'*Y - C'*W51 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W51 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H52 = [ -Y , Y*Ad52 - W52'*C , Y*H' , W52';
        Ad52'*Y - C'*W52 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W52 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H53 = [ -Y , Y*Ad53 - W53'*C , Y*H' , W53';
        Ad53'*Y - C'*W53 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W53 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H54 = [ -Y , Y*Ad54 - W54'*C , Y*H' , W54';
        Ad54'*Y - C'*W54 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W54 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H55 = [ -Y , Y*Ad55 - W55'*C , Y*H' , W55';
        Ad55'*Y - C'*W55 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W55 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H56 = [ -Y , Y*Ad56 - W56'*C , Y*H' , W56';
        Ad56'*Y - C'*W56 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W56 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H57 = [ -Y , Y*Ad57 - W57'*C , Y*H' , W57';
        Ad57'*Y - C'*W57 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W57 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H58 = [ -Y , Y*Ad58 - W58'*C , Y*H' , W58';
        Ad58'*Y - C'*W58 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W58 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H59 = [ -Y , Y*Ad59 - W59'*C , Y*H' , W59';
        Ad59'*Y - C'*W59 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W59 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H60 = [ -Y , Y*Ad60 - W60'*C , Y*H' , W60';
        Ad60'*Y - C'*W60 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W60 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H61 = [ -Y , Y*Ad61 - W61'*C , Y*H' , W61';
        Ad61'*Y - C'*W61 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W61 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H62 = [ -Y , Y*Ad62 - W62'*C , Y*H' , W62';
        Ad62'*Y - C'*W62 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W62 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H63 = [ -Y , Y*Ad63 - W63'*C , Y*H' , W63';
        Ad63'*Y - C'*W63 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W63 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H64 = [ -Y , Y*Ad64 - W64'*C , Y*H' , W64';
        Ad64'*Y - C'*W64 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W64 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H65 = [ -Y , Y*Ad65 - W65'*C , Y*H' , W65';
        Ad65'*Y - C'*W65 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W65 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H66 = [ -Y , Y*Ad66 - W66'*C , Y*H' , W66';
        Ad66'*Y - C'*W66 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W66 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H67 = [ -Y , Y*Ad67 - W67'*C , Y*H' , W67';
        Ad67'*Y - C'*W67 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W67 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H68 = [ -Y , Y*Ad68 - W68'*C , Y*H' , W68';
        Ad68'*Y - C'*W68 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W68 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H69 = [ -Y , Y*Ad69 - W69'*C , Y*H' , W69';
        Ad69'*Y - C'*W69 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W69 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H70 = [ -Y , Y*Ad70 - W70'*C , Y*H' , W70';
        Ad70'*Y - C'*W70 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W70 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H71 = [ -Y , Y*Ad71 - W71'*C , Y*H' , W71';
        Ad71'*Y - C'*W71 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W71 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H72 = [ -Y , Y*Ad72 - W72'*C , Y*H' , W72';
        Ad72'*Y - C'*W72 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W72 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H73 = [ -Y , Y*Ad73 - W73'*C , Y*H' , W73';
        Ad73'*Y - C'*W73 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W73 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H74 = [ -Y , Y*Ad74 - W74'*C , Y*H' , W74';
        Ad74'*Y - C'*W74 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W74 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H75 = [ -Y , Y*Ad75 - W75'*C , Y*H' , W75';
        Ad75'*Y - C'*W75 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W75 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H76 = [ -Y , Y*Ad76 - W76'*C , Y*H' , W76';
        Ad76'*Y - C'*W76 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W76 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H77 = [ -Y , Y*Ad77 - W77'*C , Y*H' , W77';
        Ad77'*Y - C'*W77 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W77 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H78 = [ -Y , Y*Ad78 - W78'*C , Y*H' , W78';
        Ad78'*Y - C'*W78 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W78 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H79 = [ -Y , Y*Ad79 - W79'*C , Y*H' , W79';
        Ad79'*Y - C'*W79 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W79 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H80 = [ -Y , Y*Ad80 - W80'*C , Y*H' , W80';
        Ad80'*Y - C'*W80 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W80 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H81 = [ -Y , Y*Ad81 - W81'*C , Y*H' , W81';
        Ad81'*Y - C'*W81 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W81 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H82 = [ -Y , Y*Ad82 - W82'*C , Y*H' , W82';
        Ad82'*Y - C'*W82 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W82 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H83 = [ -Y , Y*Ad83 - W83'*C , Y*H' , W83';
        Ad83'*Y - C'*W83 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W83 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H84 = [ -Y , Y*Ad84 - W84'*C , Y*H' , W84';
        Ad84'*Y - C'*W84 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W84 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H85 = [ -Y , Y*Ad85 - W85'*C , Y*H' , W85';
        Ad85'*Y - C'*W85 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W85 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H86 = [ -Y , Y*Ad86 - W86'*C , Y*H' , W86';
        Ad86'*Y - C'*W86 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W86 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H87 = [ -Y , Y*Ad87 - W87'*C , Y*H' , W87';
        Ad87'*Y - C'*W87 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W87 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H88 = [ -Y , Y*Ad88 - W88'*C , Y*H' , W88';
        Ad88'*Y - C'*W88 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W88 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H89 = [ -Y , Y*Ad89 - W89'*C , Y*H' , W89';
        Ad89'*Y - C'*W89 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W89 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H90 = [ -Y , Y*Ad90 - W90'*C , Y*H' , W90';
        Ad90'*Y - C'*W90 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W90 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H91 = [ -Y , Y*Ad91 - W91'*C , Y*H' , W91';
        Ad91'*Y - C'*W91 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W91 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H92 = [ -Y , Y*Ad92 - W92'*C , Y*H' , W92';
        Ad92'*Y - C'*W92 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W92 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H93 = [ -Y , Y*Ad93 - W93'*C , Y*H' , W93';
        Ad93'*Y - C'*W93 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W93 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H94 = [ -Y , Y*Ad94 - W94'*C , Y*H' , W94';
        Ad94'*Y - C'*W94 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W94 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H95 = [ -Y , Y*Ad95 - W95'*C , Y*H' , W95';
        Ad95'*Y - C'*W95 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W95 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H96 = [ -Y , Y*Ad96 - W96'*C , Y*H' , W96';
        Ad96'*Y - C'*W96 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W96 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H97 = [ -Y , Y*Ad97 - W97'*C , Y*H' , W97';
        Ad97'*Y - C'*W97 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W97 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H98 = [ -Y , Y*Ad98 - W98'*C , Y*H' , W98';
        Ad98'*Y - C'*W98 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W98 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H99 = [ -Y , Y*Ad99 - W99'*C , Y*H' , W99';
        Ad99'*Y - C'*W99 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W99 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H100 = [ -Y , Y*Ad100 - W100'*C , Y*H' , W100';
        Ad100'*Y - C'*W100 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W100 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H101 = [ -Y , Y*Ad101 - W101'*C , Y*H' , W101';
        Ad101'*Y - C'*W101 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W101 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H102 = [ -Y , Y*Ad102 - W102'*C , Y*H' , W102';
        Ad102'*Y - C'*W102 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W102 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H103 = [ -Y , Y*Ad103 - W103'*C , Y*H' , W103';
        Ad103'*Y - C'*W103 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W103 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H104 = [ -Y , Y*Ad104 - W104'*C , Y*H' , W104';
        Ad104'*Y - C'*W104 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W104 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H105 = [ -Y , Y*Ad105 - W105'*C , Y*H' , W105';
        Ad105'*Y - C'*W105 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W105 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H106 = [ -Y , Y*Ad106 - W106'*C , Y*H' , W106';
        Ad106'*Y - C'*W106 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W106 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H107 = [ -Y , Y*Ad107 - W107'*C , Y*H' , W107';
        Ad107'*Y - C'*W107 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W107 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H108 = [ -Y , Y*Ad108 - W108'*C , Y*H' , W108';
        Ad108'*Y - C'*W108 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W108 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H109 = [ -Y , Y*Ad109 - W109'*C , Y*H' , W109';
        Ad109'*Y - C'*W109 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W109 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H110 = [ -Y , Y*Ad110 - W110'*C , Y*H' , W110';
        Ad110'*Y - C'*W110 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W110 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H111 = [ -Y , Y*Ad111 - W111'*C , Y*H' , W111';
        Ad111'*Y - C'*W111 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W111 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H112 = [ -Y , Y*Ad112 - W112'*C , Y*H' , W112';
        Ad112'*Y - C'*W112 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W112 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H113 = [ -Y , Y*Ad113 - W113'*C , Y*H' , W113';
        Ad113'*Y - C'*W113 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W113 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H114 = [ -Y , Y*Ad114 - W114'*C , Y*H' , W114';
        Ad114'*Y - C'*W114 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W114 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H115 = [ -Y , Y*Ad115 - W115'*C , Y*H' , W115';
        Ad115'*Y - C'*W115 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W115 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H116 = [ -Y , Y*Ad116 - W116'*C , Y*H' , W116';
        Ad116'*Y - C'*W116 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W116 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H117 = [ -Y , Y*Ad117 - W117'*C , Y*H' , W117';
        Ad117'*Y - C'*W117 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W117 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H118 = [ -Y , Y*Ad118 - W118'*C , Y*H' , W118';
        Ad118'*Y - C'*W118 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W118 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H119 = [ -Y , Y*Ad119 - W119'*C , Y*H' , W119';
        Ad119'*Y - C'*W119 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W119 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H120 = [ -Y , Y*Ad120 - W120'*C , Y*H' , W120';
        Ad120'*Y - C'*W120 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W120 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H121 = [ -Y , Y*Ad121 - W121'*C , Y*H' , W121';
        Ad121'*Y - C'*W121 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W121 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H122 = [ -Y , Y*Ad122 - W122'*C , Y*H' , W122';
        Ad122'*Y - C'*W122 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W122 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H123 = [ -Y , Y*Ad123 - W123'*C , Y*H' , W123';
        Ad123'*Y - C'*W123 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W123 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H124 = [ -Y , Y*Ad124 - W124'*C , Y*H' , W124';
        Ad124'*Y - C'*W124 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W124 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H125 = [ -Y , Y*Ad125 - W125'*C , Y*H' , W125';
        Ad125'*Y - C'*W125 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W125 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H126 = [ -Y , Y*Ad126 - W126'*C , Y*H' , W126';
        Ad126'*Y - C'*W126 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W126 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H127 = [ -Y , Y*Ad127 - W127'*C , Y*H' , W127';
        Ad127'*Y - C'*W127 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W127 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H128 = [ -Y , Y*Ad128 - W128'*C , Y*H' , W128';
        Ad128'*Y - C'*W128 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W128 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H129 = [ -Y , Y*Ad129 - W129'*C , Y*H' , W129';
        Ad129'*Y - C'*W129 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W129 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H130 = [ -Y , Y*Ad130 - W130'*C , Y*H' , W130';
        Ad130'*Y - C'*W130 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W130 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H131 = [ -Y , Y*Ad131 - W131'*C , Y*H' , W131';
        Ad131'*Y - C'*W131 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W131 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H132 = [ -Y , Y*Ad132 - W132'*C , Y*H' , W132';
        Ad132'*Y - C'*W132 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W132 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H133 = [ -Y , Y*Ad133 - W133'*C , Y*H' , W133';
        Ad133'*Y - C'*W133 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W133 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H134 = [ -Y , Y*Ad134 - W134'*C , Y*H' , W134';
        Ad134'*Y - C'*W134 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W134 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H135 = [ -Y , Y*Ad135 - W135'*C , Y*H' , W135';
        Ad135'*Y - C'*W135 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W135 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H136 = [ -Y , Y*Ad136 - W136'*C , Y*H' , W136';
        Ad136'*Y - C'*W136 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W136 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H137 = [ -Y , Y*Ad137 - W137'*C , Y*H' , W137';
        Ad137'*Y - C'*W137 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W137 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H138 = [ -Y , Y*Ad138 - W138'*C , Y*H' , W138';
        Ad138'*Y - C'*W138 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W138 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H139 = [ -Y , Y*Ad139 - W139'*C , Y*H' , W139';
        Ad139'*Y - C'*W139 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W139 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H140 = [ -Y , Y*Ad140 - W140'*C , Y*H' , W140';
        Ad140'*Y - C'*W140 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W140 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H141 = [ -Y , Y*Ad141 - W141'*C , Y*H' , W141';
        Ad141'*Y - C'*W141 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W141 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H142 = [ -Y , Y*Ad142 - W142'*C , Y*H' , W142';
        Ad142'*Y - C'*W142 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W142 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H143 = [ -Y , Y*Ad143 - W143'*C , Y*H' , W143';
        Ad143'*Y - C'*W143 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W143 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H144 = [ -Y , Y*Ad144 - W144'*C , Y*H' , W144';
        Ad144'*Y - C'*W144 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W144 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H145 = [ -Y , Y*Ad145 - W145'*C , Y*H' , W145';
        Ad145'*Y - C'*W145 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W145 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H146 = [ -Y , Y*Ad146 - W146'*C , Y*H' , W146';
        Ad146'*Y - C'*W146 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W146 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H147 = [ -Y , Y*Ad147 - W147'*C , Y*H' , W147';
        Ad147'*Y - C'*W147 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W147 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H148 = [ -Y , Y*Ad148 - W148'*C , Y*H' , W148';
        Ad148'*Y - C'*W148 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W148 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H149 = [ -Y , Y*Ad149 - W149'*C , Y*H' , W149';
        Ad149'*Y - C'*W149 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W149 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H150 = [ -Y , Y*Ad150 - W150'*C , Y*H' , W150';
        Ad150'*Y - C'*W150 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W150 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H151 = [ -Y , Y*Ad151 - W151'*C , Y*H' , W151';
        Ad151'*Y - C'*W151 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W151 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H152 = [ -Y , Y*Ad152 - W152'*C , Y*H' , W152';
        Ad152'*Y - C'*W152 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W152 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H153 = [ -Y , Y*Ad153 - W153'*C , Y*H' , W153';
        Ad153'*Y - C'*W153 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W153 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H154 = [ -Y , Y*Ad154 - W154'*C , Y*H' , W154';
        Ad154'*Y - C'*W154 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W154 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H155 = [ -Y , Y*Ad155 - W155'*C , Y*H' , W155';
        Ad155'*Y - C'*W155 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W155 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H156 = [ -Y , Y*Ad156 - W156'*C , Y*H' , W156';
        Ad156'*Y - C'*W156 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W156 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H157 = [ -Y , Y*Ad157 - W157'*C , Y*H' , W157';
        Ad157'*Y - C'*W157 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W157 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H158 = [ -Y , Y*Ad158 - W158'*C , Y*H' , W158';
        Ad158'*Y - C'*W158 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W158 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H159 = [ -Y , Y*Ad159 - W159'*C , Y*H' , W159';
        Ad159'*Y - C'*W159 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W159 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H160 = [ -Y , Y*Ad160 - W160'*C , Y*H' , W160';
        Ad160'*Y - C'*W160 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W160 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H161 = [ -Y , Y*Ad161 - W161'*C , Y*H' , W161';
        Ad161'*Y - C'*W161 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W161 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H162 = [ -Y , Y*Ad162 - W162'*C , Y*H' , W162';
        Ad162'*Y - C'*W162 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W162 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H163 = [ -Y , Y*Ad163 - W163'*C , Y*H' , W163';
        Ad163'*Y - C'*W163 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W163 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H164 = [ -Y , Y*Ad164 - W164'*C , Y*H' , W164';
        Ad164'*Y - C'*W164 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W164 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H165 = [ -Y , Y*Ad165 - W165'*C , Y*H' , W165';
        Ad165'*Y - C'*W165 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W165 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H166 = [ -Y , Y*Ad166 - W166'*C , Y*H' , W166';
        Ad166'*Y - C'*W166 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W166 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H167 = [ -Y , Y*Ad167 - W167'*C , Y*H' , W167';
        Ad167'*Y - C'*W167 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W167 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H168 = [ -Y , Y*Ad168 - W168'*C , Y*H' , W168';
        Ad168'*Y - C'*W168 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W168 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H169 = [ -Y , Y*Ad169 - W169'*C , Y*H' , W169';
        Ad169'*Y - C'*W169 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W169 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H170 = [ -Y , Y*Ad170 - W170'*C , Y*H' , W170';
        Ad170'*Y - C'*W170 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W170 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H171 = [ -Y , Y*Ad171 - W171'*C , Y*H' , W171';
        Ad171'*Y - C'*W171 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W171 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H172 = [ -Y , Y*Ad172 - W172'*C , Y*H' , W172';
        Ad172'*Y - C'*W172 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W172 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H173 = [ -Y , Y*Ad173 - W173'*C , Y*H' , W173';
        Ad173'*Y - C'*W173 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W173 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H174 = [ -Y , Y*Ad174 - W174'*C , Y*H' , W174';
        Ad174'*Y - C'*W174 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W174 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H175 = [ -Y , Y*Ad175 - W175'*C , Y*H' , W175';
        Ad175'*Y - C'*W175 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W175 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H176 = [ -Y , Y*Ad176 - W176'*C , Y*H' , W176';
        Ad176'*Y - C'*W176 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W176 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H177 = [ -Y , Y*Ad177 - W177'*C , Y*H' , W177';
        Ad177'*Y - C'*W177 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W177 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H178 = [ -Y , Y*Ad178 - W178'*C , Y*H' , W178';
        Ad178'*Y - C'*W178 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W178 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H179 = [ -Y , Y*Ad179 - W179'*C , Y*H' , W179';
        Ad179'*Y - C'*W179 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W179 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H180 = [ -Y , Y*Ad180 - W180'*C , Y*H' , W180';
        Ad180'*Y - C'*W180 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W180 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H181 = [ -Y , Y*Ad181 - W181'*C , Y*H' , W181';
        Ad181'*Y - C'*W181 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W181 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H182 = [ -Y , Y*Ad182 - W182'*C , Y*H' , W182';
        Ad182'*Y - C'*W182 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W182 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H183 = [ -Y , Y*Ad183 - W183'*C , Y*H' , W183';
        Ad183'*Y - C'*W183 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W183 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H184 = [ -Y , Y*Ad184 - W184'*C , Y*H' , W184';
        Ad184'*Y - C'*W184 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W184 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H185 = [ -Y , Y*Ad185 - W185'*C , Y*H' , W185';
        Ad185'*Y - C'*W185 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W185 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H186 = [ -Y , Y*Ad186 - W186'*C , Y*H' , W186';
        Ad186'*Y - C'*W186 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W186 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H187 = [ -Y , Y*Ad187 - W187'*C , Y*H' , W187';
        Ad187'*Y - C'*W187 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W187 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H188 = [ -Y , Y*Ad188 - W188'*C , Y*H' , W188';
        Ad188'*Y - C'*W188 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W188 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H189 = [ -Y , Y*Ad189 - W189'*C , Y*H' , W189';
        Ad189'*Y - C'*W189 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W189 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H190 = [ -Y , Y*Ad190 - W190'*C , Y*H' , W190';
        Ad190'*Y - C'*W190 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W190 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H191 = [ -Y , Y*Ad191 - W191'*C , Y*H' , W191';
        Ad191'*Y - C'*W191 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W191 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H192 = [ -Y , Y*Ad192 - W192'*C , Y*H' , W192';
        Ad192'*Y - C'*W192 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W192 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H193 = [ -Y , Y*Ad193 - W193'*C , Y*H' , W193';
        Ad193'*Y - C'*W193 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W193 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H194 = [ -Y , Y*Ad194 - W194'*C , Y*H' , W194';
        Ad194'*Y - C'*W194 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W194 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H195 = [ -Y , Y*Ad195 - W195'*C , Y*H' , W195';
        Ad195'*Y - C'*W195 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W195 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H196 = [ -Y , Y*Ad196 - W196'*C , Y*H' , W196';
        Ad196'*Y - C'*W196 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W196 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H197 = [ -Y , Y*Ad197 - W197'*C , Y*H' , W197';
        Ad197'*Y - C'*W197 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W197 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H198 = [ -Y , Y*Ad198 - W198'*C , Y*H' , W198';
        Ad198'*Y - C'*W198 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W198 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H199 = [ -Y , Y*Ad199 - W199'*C , Y*H' , W199';
        Ad199'*Y - C'*W199 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W199 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H200 = [ -Y , Y*Ad200 - W200'*C , Y*H' , W200';
        Ad200'*Y - C'*W200 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W200 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H201 = [ -Y , Y*Ad201 - W201'*C , Y*H' , W201';
        Ad201'*Y - C'*W201 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W201 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H202 = [ -Y , Y*Ad202 - W202'*C , Y*H' , W202';
        Ad202'*Y - C'*W202 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W202 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H203 = [ -Y , Y*Ad203 - W203'*C , Y*H' , W203';
        Ad203'*Y - C'*W203 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W203 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H204 = [ -Y , Y*Ad204 - W204'*C , Y*H' , W204';
        Ad204'*Y - C'*W204 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W204 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H205 = [ -Y , Y*Ad205 - W205'*C , Y*H' , W205';
        Ad205'*Y - C'*W205 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W205 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H206 = [ -Y , Y*Ad206 - W206'*C , Y*H' , W206';
        Ad206'*Y - C'*W206 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W206 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H207 = [ -Y , Y*Ad207 - W207'*C , Y*H' , W207';
        Ad207'*Y - C'*W207 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W207 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H208 = [ -Y , Y*Ad208 - W208'*C , Y*H' , W208';
        Ad208'*Y - C'*W208 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W208 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H209 = [ -Y , Y*Ad209 - W209'*C , Y*H' , W209';
        Ad209'*Y - C'*W209 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W209 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H210 = [ -Y , Y*Ad210 - W210'*C , Y*H' , W210';
        Ad210'*Y - C'*W210 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W210 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H211 = [ -Y , Y*Ad211 - W211'*C , Y*H' , W211';
        Ad211'*Y - C'*W211 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W211 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H212 = [ -Y , Y*Ad212 - W212'*C , Y*H' , W212';
        Ad212'*Y - C'*W212 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W212 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H213 = [ -Y , Y*Ad213 - W213'*C , Y*H' , W213';
        Ad213'*Y - C'*W213 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W213 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H214 = [ -Y , Y*Ad214 - W214'*C , Y*H' , W214';
        Ad214'*Y - C'*W214 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W214 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H215 = [ -Y , Y*Ad215 - W215'*C , Y*H' , W215';
        Ad215'*Y - C'*W215 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W215 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H216 = [ -Y , Y*Ad216 - W216'*C , Y*H' , W216';
        Ad216'*Y - C'*W216 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W216 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H217 = [ -Y , Y*Ad217 - W217'*C , Y*H' , W217';
        Ad217'*Y - C'*W217 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W217 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H218 = [ -Y , Y*Ad218 - W218'*C , Y*H' , W218';
        Ad218'*Y - C'*W218 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W218 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H219 = [ -Y , Y*Ad219 - W219'*C , Y*H' , W219';
        Ad219'*Y - C'*W219 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W219 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H220 = [ -Y , Y*Ad220 - W220'*C , Y*H' , W220';
        Ad220'*Y - C'*W220 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W220 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H221 = [ -Y , Y*Ad221 - W221'*C , Y*H' , W221';
        Ad221'*Y - C'*W221 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W221 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H222 = [ -Y , Y*Ad222 - W222'*C , Y*H' , W222';
        Ad222'*Y - C'*W222 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W222 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H223 = [ -Y , Y*Ad223 - W223'*C , Y*H' , W223';
        Ad223'*Y - C'*W223 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W223 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H224 = [ -Y , Y*Ad224 - W224'*C , Y*H' , W224';
        Ad224'*Y - C'*W224 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W224 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H225 = [ -Y , Y*Ad225 - W225'*C , Y*H' , W225';
        Ad225'*Y - C'*W225 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W225 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H226 = [ -Y , Y*Ad226 - W226'*C , Y*H' , W226';
        Ad226'*Y - C'*W226 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W226 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H227 = [ -Y , Y*Ad227 - W227'*C , Y*H' , W227';
        Ad227'*Y - C'*W227 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W227 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H228 = [ -Y , Y*Ad228 - W228'*C , Y*H' , W228';
        Ad228'*Y - C'*W228 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W228 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H229 = [ -Y , Y*Ad229 - W229'*C , Y*H' , W229';
        Ad229'*Y - C'*W229 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W229 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H230 = [ -Y , Y*Ad230 - W230'*C , Y*H' , W230';
        Ad230'*Y - C'*W230 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W230 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H231 = [ -Y , Y*Ad231 - W231'*C , Y*H' , W231';
        Ad231'*Y - C'*W231 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W231 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H232 = [ -Y , Y*Ad232 - W232'*C , Y*H' , W232';
        Ad232'*Y - C'*W232 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W232 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H233 = [ -Y , Y*Ad233 - W233'*C , Y*H' , W233';
        Ad233'*Y - C'*W233 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W233 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H234 = [ -Y , Y*Ad234 - W234'*C , Y*H' , W234';
        Ad234'*Y - C'*W234 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W234 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H235 = [ -Y , Y*Ad235 - W235'*C , Y*H' , W235';
        Ad235'*Y - C'*W235 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W235 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H236 = [ -Y , Y*Ad236 - W236'*C , Y*H' , W236';
        Ad236'*Y - C'*W236 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W236 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H237 = [ -Y , Y*Ad237 - W237'*C , Y*H' , W237';
        Ad237'*Y - C'*W237 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W237 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H238 = [ -Y , Y*Ad238 - W238'*C , Y*H' , W238';
        Ad238'*Y - C'*W238 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W238 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H239 = [ -Y , Y*Ad239 - W239'*C , Y*H' , W239';
        Ad239'*Y - C'*W239 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W239 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H240 = [ -Y , Y*Ad240 - W240'*C , Y*H' , W240';
        Ad240'*Y - C'*W240 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W240 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H241 = [ -Y , Y*Ad241 - W241'*C , Y*H' , W241';
        Ad241'*Y - C'*W241 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W241 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H242 = [ -Y , Y*Ad242 - W242'*C , Y*H' , W242';
        Ad242'*Y - C'*W242 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W242 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H243 = [ -Y , Y*Ad243 - W243'*C , Y*H' , W243';
        Ad243'*Y - C'*W243 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W243 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H244 = [ -Y , Y*Ad244 - W244'*C , Y*H' , W244';
        Ad244'*Y - C'*W244 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W244 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H245 = [ -Y , Y*Ad245 - W245'*C , Y*H' , W245';
        Ad245'*Y - C'*W245 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W245 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H246 = [ -Y , Y*Ad246 - W246'*C , Y*H' , W246';
        Ad246'*Y - C'*W246 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W246 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H247 = [ -Y , Y*Ad247 - W247'*C , Y*H' , W247';
        Ad247'*Y - C'*W247 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W247 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H248 = [ -Y , Y*Ad248 - W248'*C , Y*H' , W248';
        Ad248'*Y - C'*W248 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W248 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H249 = [ -Y , Y*Ad249 - W249'*C , Y*H' , W249';
        Ad249'*Y - C'*W249 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W249 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H250 = [ -Y , Y*Ad250 - W250'*C , Y*H' , W250';
        Ad250'*Y - C'*W250 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W250 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H251 = [ -Y , Y*Ad251 - W251'*C , Y*H' , W251';
        Ad251'*Y - C'*W251 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W251 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H252 = [ -Y , Y*Ad252 - W252'*C , Y*H' , W252';
        Ad252'*Y - C'*W252 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W252 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H253 = [ -Y , Y*Ad253 - W253'*C , Y*H' , W253';
        Ad253'*Y - C'*W253 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W253 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H254 = [ -Y , Y*Ad254 - W254'*C , Y*H' , W254';
        Ad254'*Y - C'*W254 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W254 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H255 = [ -Y , Y*Ad255 - W255'*C , Y*H' , W255';
        Ad255'*Y - C'*W255 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W255 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];

H256 = [ -Y , Y*Ad256 - W256'*C , Y*H' , W256';
        Ad256'*Y - C'*W256 , -Y , zeros(9,9) , zeros(9,2);
        H*Y , zeros(9,9) , -eye(9) , zeros(9,2);
        W256 , zeros(2,9) , zeros(2,9) , -inv(R_obs)];



tic
ops= sdpsettings('solver','mosek','verbose',0);

[opt_info] = optimize([[H1<=0] + [H2<=0] + [H3<=0] + [H4<=0] + [H5<=0] + [H6<=0] + [H7<=0] + [H8<=0] + [H9<=0] + [H10<=0] + [H11<=0] + [H12<=0] + [H13<=0] + [H14<=0] + [H15<=0] + [H16<=0] + ...
[H17<=0] + [H18<=0] + [H19<=0] + [H20<=0] + [H21<=0] + [H22<=0] + [H23<=0] + [H24<=0] + [H25<=0] + [H26<=0] + [H27<=0] + [H28<=0] + [H29<=0] + [H30<=0] + [H31<=0] + [H32<=0] + ...
[H33<=0] + [H34<=0] + [H35<=0] + [H36<=0] + [H37<=0] + [H38<=0] + [H39<=0] + [H40<=0] + [H41<=0] + [H42<=0] + [H43<=0] + [H44<=0] + [H45<=0] + [H46<=0] + [H47<=0] + [H48<=0] + ...
[H49<=0] + [H50<=0] + [H51<=0] + [H52<=0] + [H53<=0] + [H54<=0] + [H55<=0] + [H56<=0] + [H57<=0] + [H58<=0] + [H59<=0] + [H60<=0] + [H61<=0] + [H62<=0] + [H63<=0] + [H64<=0] + ...
[H65<=0] + [H66<=0] + [H67<=0] + [H68<=0] + [H69<=0] + [H70<=0] + [H71<=0] + [H72<=0] + [H73<=0] + [H74<=0] + [H75<=0] + [H76<=0] + [H77<=0] + [H78<=0] + [H79<=0] + [H80<=0] + ...
[H81<=0] + [H82<=0] + [H83<=0] + [H84<=0] + [H85<=0] + [H86<=0] + [H87<=0] + [H88<=0] + [H89<=0] + [H90<=0] + [H91<=0] + [H92<=0] + [H93<=0] + [H94<=0] + [H95<=0] + [H96<=0] + ...
[H97<=0] + [H98<=0] + [H99<=0] + [H100<=0] + [H101<=0] + [H102<=0] + [H103<=0] + [H104<=0] + [H105<=0] + [H106<=0] + [H107<=0] + [H108<=0] + [H109<=0] + [H110<=0] + [H111<=0] + [H112<=0] + ...
[H113<=0] + [H114<=0] + [H115<=0] + [H116<=0] + [H117<=0] + [H118<=0] + [H119<=0] + [H120<=0] + [H121<=0] + [H122<=0] + [H123<=0] + [H124<=0] + [H125<=0] + [H126<=0] + [H127<=0] + [H128<=0] + ...
[H129<=0] + [H130<=0] + [H131<=0] + [H132<=0] + [H133<=0] + [H134<=0] + [H135<=0] + [H136<=0] + [H137<=0] + [H138<=0] + [H139<=0] + [H140<=0] + [H141<=0] + [H142<=0] + [H143<=0] + [H144<=0] + ...
[H145<=0] + [H146<=0] + [H147<=0] + [H148<=0] + [H149<=0] + [H150<=0] + [H151<=0] + [H152<=0] + [H153<=0] + [H154<=0] + [H155<=0] + [H156<=0] + [H157<=0] + [H158<=0] + [H159<=0] + [H160<=0] + ...
[H161<=0] + [H162<=0] + [H163<=0] + [H164<=0] + [H165<=0] + [H166<=0] + [H167<=0] + [H168<=0] + [H169<=0] + [H170<=0] + [H171<=0] + [H172<=0] + [H173<=0] + [H174<=0] + [H175<=0] + [H176<=0] + ...
[H177<=0] + [H178<=0] + [H179<=0] + [H180<=0] + [H181<=0] + [H182<=0] + [H183<=0] + [H184<=0] + [H185<=0] + [H186<=0] + [H187<=0] + [H188<=0] + [H189<=0] + [H190<=0] + [H191<=0] + [H192<=0] + ...
[H193<=0] + [H194<=0] + [H195<=0] + [H196<=0] + [H197<=0] + [H198<=0] + [H199<=0] + [H200<=0] + [H201<=0] + [H202<=0] + [H203<=0] + [H204<=0] + [H205<=0] + [H206<=0] + [H207<=0] + [H208<=0] + ...
[H209<=0] + [H210<=0] + [H211<=0] + [H212<=0] + [H213<=0] + [H214<=0] + [H215<=0] + [H216<=0] + [H217<=0] + [H218<=0] + [H219<=0] + [H220<=0] + [H221<=0] + [H222<=0] + [H223<=0] + [H224<=0] + ...
[H225<=0] + [H226<=0] + [H227<=0] + [H228<=0] + [H229<=0] + [H230<=0] + [H231<=0] + [H232<=0] + [H233<=0] + [H234<=0] + [H235<=0] + [H236<=0] + [H237<=0] + [H238<=0] + [H239<=0] + [H240<=0] + ...
[H241<=0] + [H242<=0] + [H243<=0] + [H244<=0] + [H245<=0] + [H246<=0] + [H247<=0] + [H248<=0] + [H249<=0] + [H250<=0] + [H251<=0] + [H252<=0] + [H253<=0] + [H254<=0] + [H255<=0] + [H256<=0] + ...
[Y>=0] + [Trace_cond>=0] ], [gamma_LQR_obs],ops);

tic
tic
Y_sol  = value(Y);%  Y solution -> WARNING: LQR solves for -Y!!
W1_sol  = value(W1); % W1 solution
W2_sol  = value(W2); % W2 solution
W3_sol  = value(W3); % W3 solution
W4_sol  = value(W4); % W4 solution
W5_sol  = value(W5); % W5 solution
W6_sol  = value(W6); % W6 solution
W7_sol  = value(W7); % W7 solution
W8_sol  = value(W8); % W8 solution
W9_sol  = value(W9); % W9 solution
W10_sol = value(W10); % W10 solution
W11_sol = value(W11); % W11 solution
W12_sol = value(W12); % W12 solution
W13_sol = value(W13); % W13 solution
W14_sol = value(W14); % W14 solution
W15_sol = value(W15); % W15 solution
W16_sol = value(W16); % W16 solution
W17_sol = value(W17); % W17 solution
W18_sol = value(W18); % W18 solution
W19_sol = value(W19); % W19 solution
W20_sol = value(W20); % W20 solution
W21_sol = value(W21); % W21 solution
W22_sol = value(W22); % W22 solution
W23_sol = value(W23); % W23 solution
W24_sol = value(W24); % W24 solution
W25_sol = value(W25); % W25 solution
W26_sol = value(W26); % W26 solution
W27_sol = value(W27); % W27 solution
W28_sol = value(W28); % W28 solution
W29_sol = value(W29); % W29 solution
W30_sol = value(W30); % W30 solution
W31_sol = value(W31); % W31 solution
W32_sol = value(W32); % W32 solution
W33_sol = value(W33); % W33 solution
W34_sol = value(W34); % W34 solution
W35_sol = value(W35); % W35 solution
W36_sol = value(W36); % W36 solution
W37_sol = value(W37); % W37 solution
W38_sol = value(W38); % W38 solution
W39_sol = value(W39); % W39 solution
W40_sol = value(W40); % W40 solution
W41_sol = value(W41); % W41 solution
W42_sol = value(W42); % W42 solution
W43_sol = value(W43); % W43 solution
W44_sol = value(W44); % W44 solution
W45_sol = value(W45); % W45 solution
W46_sol = value(W46); % W46 solution
W47_sol = value(W47); % W47 solution
W48_sol = value(W48); % W48 solution
W49_sol = value(W49); % W49 solution
W50_sol = value(W50); % W50 solution
W51_sol = value(W51); % W51 solution
W52_sol = value(W52); % W52 solution
W53_sol = value(W53); % W53 solution
W54_sol = value(W54); % W54 solution
W55_sol = value(W55); % W55 solution
W56_sol = value(W56); % W56 solution
W57_sol = value(W57); % W57 solution
W58_sol = value(W58); % W58 solution
W59_sol = value(W59); % W59 solution
W60_sol = value(W60); % W60 solution
W61_sol = value(W61); % W61 solution
W62_sol = value(W62); % W62 solution
W63_sol = value(W63); % W63 solution
W64_sol = value(W64); % W64 solution
W65_sol = value(W65); % W65 solution
W66_sol = value(W66); % W66 solution
W67_sol = value(W67); % W67 solution
W68_sol = value(W68); % W68 solution
W69_sol = value(W69); % W69 solution
W70_sol = value(W70); % W70 solution
W71_sol = value(W71); % W71 solution
W72_sol = value(W72); % W72 solution
W73_sol = value(W73); % W73 solution
W74_sol = value(W74); % W74 solution
W75_sol = value(W75); % W75 solution
W76_sol = value(W76); % W76 solution
W77_sol = value(W77); % W77 solution
W78_sol = value(W78); % W78 solution
W79_sol = value(W79); % W79 solution
W80_sol = value(W80); % W80 solution
W81_sol = value(W81); % W81 solution
W82_sol = value(W82); % W82 solution
W83_sol = value(W83); % W83 solution
W84_sol = value(W84); % W84 solution
W85_sol = value(W85); % W85 solution
W86_sol = value(W86); % W86 solution
W87_sol = value(W87); % W87 solution
W88_sol = value(W88); % W88 solution
W89_sol = value(W89); % W89 solution
W90_sol = value(W90); % W90 solution
W91_sol = value(W91); % W91 solution
W92_sol = value(W92); % W92 solution
W93_sol = value(W93); % W93 solution
W94_sol = value(W94); % W94 solution
W95_sol = value(W95); % W95 solution
W96_sol = value(W96); % W96 solution
W97_sol = value(W97); % W97 solution
W98_sol = value(W98); % W98 solution
W99_sol = value(W99); % W99 solution
W100_sol = value(W100); % W100 solution
W101_sol = value(W101); % W101 solution
W102_sol = value(W102); % W102 solution
W103_sol = value(W103); % W103 solution
W104_sol = value(W104); % W104 solution
W105_sol = value(W105); % W105 solution
W106_sol = value(W106); % W106 solution
W107_sol = value(W107); % W107 solution
W108_sol = value(W108); % W108 solution
W109_sol = value(W109); % W109 solution
W110_sol = value(W110); % W110 solution
W111_sol = value(W111); % W111 solution
W112_sol = value(W112); % W112 solution
W113_sol = value(W113); % W113 solution
W114_sol = value(W114); % W114 solution
W115_sol = value(W115); % W115 solution
W116_sol = value(W116); % W116 solution
W117_sol = value(W117); % W117 solution
W118_sol = value(W118); % W118 solution
W119_sol = value(W119); % W119 solution
W120_sol = value(W120); % W120 solution
W121_sol = value(W121); % W121 solution
W122_sol = value(W122); % W122 solution
W123_sol = value(W123); % W123 solution
W124_sol = value(W124); % W124 solution
W125_sol = value(W125); % W125 solution
W126_sol = value(W126); % W126 solution
W127_sol = value(W127); % W127 solution
W128_sol = value(W128); % W128 solution
W129_sol = value(W129); % W129 solution
W130_sol = value(W130); % W130 solution
W131_sol = value(W131); % W131 solution
W132_sol = value(W132); % W132 solution
W133_sol = value(W133); % W133 solution
W134_sol = value(W134); % W134 solution
W135_sol = value(W135); % W135 solution
W136_sol = value(W136); % W136 solution
W137_sol = value(W137); % W137 solution
W138_sol = value(W138); % W138 solution
W139_sol = value(W139); % W139 solution
W140_sol = value(W140); % W140 solution
W141_sol = value(W141); % W141 solution
W142_sol = value(W142); % W142 solution
W143_sol = value(W143); % W143 solution
W144_sol = value(W144); % W144 solution
W145_sol = value(W145); % W145 solution
W146_sol = value(W146); % W146 solution
W147_sol = value(W147); % W147 solution
W148_sol = value(W148); % W148 solution
W149_sol = value(W149); % W149 solution
W150_sol = value(W150); % W150 solution
W151_sol = value(W151); % W151 solution
W152_sol = value(W152); % W152 solution
W153_sol = value(W153); % W153 solution
W154_sol = value(W154); % W154 solution
W155_sol = value(W155); % W155 solution
W156_sol = value(W156); % W156 solution
W157_sol = value(W157); % W157 solution
W158_sol = value(W158); % W158 solution
W159_sol = value(W159); % W159 solution
W160_sol = value(W160); % W160 solution
W161_sol = value(W161); % W161 solution
W162_sol = value(W162); % W162 solution
W163_sol = value(W163); % W163 solution
W164_sol = value(W164); % W164 solution
W165_sol = value(W165); % W165 solution
W166_sol = value(W166); % W166 solution
W167_sol = value(W167); % W167 solution
W168_sol = value(W168); % W168 solution
W169_sol = value(W169); % W169 solution
W170_sol = value(W170); % W170 solution
W171_sol = value(W171); % W171 solution
W172_sol = value(W172); % W172 solution
W173_sol = value(W173); % W173 solution
W174_sol = value(W174); % W174 solution
W175_sol = value(W175); % W175 solution
W176_sol = value(W176); % W176 solution
W177_sol = value(W177); % W177 solution
W178_sol = value(W178); % W178 solution
W179_sol = value(W179); % W179 solution
W180_sol = value(W180); % W180 solution
W181_sol = value(W181); % W181 solution
W182_sol = value(W182); % W182 solution
W183_sol = value(W183); % W183 solution
W184_sol = value(W184); % W184 solution
W185_sol = value(W185); % W185 solution
W186_sol = value(W186); % W186 solution
W187_sol = value(W187); % W187 solution
W188_sol = value(W188); % W188 solution
W189_sol = value(W189); % W189 solution
W190_sol = value(W190); % W190 solution
W191_sol = value(W191); % W191 solution
W192_sol = value(W192); % W192 solution
W193_sol = value(W193); % W193 solution
W194_sol = value(W194); % W194 solution
W195_sol = value(W195); % W195 solution
W196_sol = value(W196); % W196 solution
W197_sol = value(W197); % W197 solution
W198_sol = value(W198); % W198 solution
W199_sol = value(W199); % W199 solution
W200_sol = value(W200); % W200 solution
W201_sol = value(W201); % W201 solution
W202_sol = value(W202); % W202 solution
W203_sol = value(W203); % W203 solution
W204_sol = value(W204); % W204 solution
W205_sol = value(W205); % W205 solution
W206_sol = value(W206); % W206 solution
W207_sol = value(W207); % W207 solution
W208_sol = value(W208); % W208 solution
W209_sol = value(W209); % W209 solution
W210_sol = value(W210); % W210 solution
W211_sol = value(W211); % W211 solution
W212_sol = value(W212); % W212 solution
W213_sol = value(W213); % W213 solution
W214_sol = value(W214); % W214 solution
W215_sol = value(W215); % W215 solution
W216_sol = value(W216); % W216 solution
W217_sol = value(W217); % W217 solution
W218_sol = value(W218); % W218 solution
W219_sol = value(W219); % W219 solution
W220_sol = value(W220); % W220 solution
W221_sol = value(W221); % W221 solution
W222_sol = value(W222); % W222 solution
W223_sol = value(W223); % W223 solution
W224_sol = value(W224); % W224 solution
W225_sol = value(W225); % W225 solution
W226_sol = value(W226); % W226 solution
W227_sol = value(W227); % W227 solution
W228_sol = value(W228); % W228 solution
W229_sol = value(W229); % W229 solution
W230_sol = value(W230); % W230 solution
W231_sol = value(W231); % W231 solution
W232_sol = value(W232); % W232 solution
W233_sol = value(W233); % W233 solution
W234_sol = value(W234); % W234 solution
W235_sol = value(W235); % W235 solution
W236_sol = value(W236); % W236 solution
W237_sol = value(W237); % W237 solution
W238_sol = value(W238); % W238 solution
W239_sol = value(W239); % W239 solution
W240_sol = value(W240); % W240 solution
W241_sol = value(W241); % W241 solution
W242_sol = value(W242); % W242 solution
W243_sol = value(W243); % W243 solution
W244_sol = value(W244); % W244 solution
W245_sol = value(W245); % W245 solution
W246_sol = value(W246); % W246 solution
W247_sol = value(W247); % W247 solution
W248_sol = value(W248); % W248 solution
W249_sol = value(W249); % W249 solution
W250_sol = value(W250); % W250 solution
W251_sol = value(W251); % W251 solution
W252_sol = value(W252); % W252 solution
W253_sol = value(W253); % W253 solution
W254_sol = value(W254); % W254 solution
W255_sol = value(W255); % W255 solution
W256_sol = value(W256); % W256 solution

L1  = (inv(Y_sol))*W1_sol.'; 
L2  = (inv(Y_sol))*W2_sol.'; 
L3  = (inv(Y_sol))*W3_sol.'; 
L4  = (inv(Y_sol))*W4_sol.'; 
L5  = (inv(Y_sol))*W5_sol.'; 
L6  = (inv(Y_sol))*W6_sol.'; 
L7  = (inv(Y_sol))*W7_sol.'; 
L8  = (inv(Y_sol))*W8_sol.'; 
L9  = (inv(Y_sol))*W9_sol.'; 
L10 = (inv(Y_sol))*W10_sol.'; 
L11 = (inv(Y_sol))*W11_sol.'; 
L12 = (inv(Y_sol))*W12_sol.'; 
L13 = (inv(Y_sol))*W13_sol.'; 
L14 = (inv(Y_sol))*W14_sol.'; 
L15 = (inv(Y_sol))*W15_sol.'; 
L16 = (inv(Y_sol))*W16_sol.'; 
L17 = (inv(Y_sol))*W17_sol.'; 
L18 = (inv(Y_sol))*W18_sol.'; 
L19 = (inv(Y_sol))*W19_sol.'; 
L20 = (inv(Y_sol))*W20_sol.'; 
L21 = (inv(Y_sol))*W21_sol.'; 
L22 = (inv(Y_sol))*W22_sol.'; 
L23 = (inv(Y_sol))*W23_sol.'; 
L24 = (inv(Y_sol))*W24_sol.'; 
L25 = (inv(Y_sol))*W25_sol.'; 
L26 = (inv(Y_sol))*W26_sol.'; 
L27 = (inv(Y_sol))*W27_sol.'; 
L28 = (inv(Y_sol))*W28_sol.'; 
L29 = (inv(Y_sol))*W29_sol.'; 
L30 = (inv(Y_sol))*W30_sol.'; 
L31 = (inv(Y_sol))*W31_sol.'; 
L32 = (inv(Y_sol))*W32_sol.'; 
L33 = (inv(Y_sol))*W33_sol.'; 
L34 = (inv(Y_sol))*W34_sol.'; 
L35 = (inv(Y_sol))*W35_sol.'; 
L36 = (inv(Y_sol))*W36_sol.'; 
L37 = (inv(Y_sol))*W37_sol.'; 
L38 = (inv(Y_sol))*W38_sol.'; 
L39 = (inv(Y_sol))*W39_sol.'; 
L40 = (inv(Y_sol))*W40_sol.'; 
L41 = (inv(Y_sol))*W41_sol.'; 
L42 = (inv(Y_sol))*W42_sol.'; 
L43 = (inv(Y_sol))*W43_sol.'; 
L44 = (inv(Y_sol))*W44_sol.'; 
L45 = (inv(Y_sol))*W45_sol.'; 
L46 = (inv(Y_sol))*W46_sol.'; 
L47 = (inv(Y_sol))*W47_sol.'; 
L48 = (inv(Y_sol))*W48_sol.'; 
L49 = (inv(Y_sol))*W49_sol.'; 
L50 = (inv(Y_sol))*W50_sol.'; 
L51 = (inv(Y_sol))*W51_sol.'; 
L52 = (inv(Y_sol))*W52_sol.'; 
L53 = (inv(Y_sol))*W53_sol.'; 
L54 = (inv(Y_sol))*W54_sol.'; 
L55 = (inv(Y_sol))*W55_sol.'; 
L56 = (inv(Y_sol))*W56_sol.'; 
L57 = (inv(Y_sol))*W57_sol.'; 
L58 = (inv(Y_sol))*W58_sol.'; 
L59 = (inv(Y_sol))*W59_sol.'; 
L60 = (inv(Y_sol))*W60_sol.'; 
L61 = (inv(Y_sol))*W61_sol.'; 
L62 = (inv(Y_sol))*W62_sol.'; 
L63 = (inv(Y_sol))*W63_sol.'; 
L64 = (inv(Y_sol))*W64_sol.'; 
L65 = (inv(Y_sol))*W65_sol.'; 
L66 = (inv(Y_sol))*W66_sol.'; 
L67 = (inv(Y_sol))*W67_sol.'; 
L68 = (inv(Y_sol))*W68_sol.'; 
L69 = (inv(Y_sol))*W69_sol.'; 
L70 = (inv(Y_sol))*W70_sol.'; 
L71 = (inv(Y_sol))*W71_sol.'; 
L72 = (inv(Y_sol))*W72_sol.'; 
L73 = (inv(Y_sol))*W73_sol.'; 
L74 = (inv(Y_sol))*W74_sol.'; 
L75 = (inv(Y_sol))*W75_sol.'; 
L76 = (inv(Y_sol))*W76_sol.'; 
L77 = (inv(Y_sol))*W77_sol.'; 
L78 = (inv(Y_sol))*W78_sol.'; 
L79 = (inv(Y_sol))*W79_sol.'; 
L80 = (inv(Y_sol))*W80_sol.'; 
L81 = (inv(Y_sol))*W81_sol.'; 
L82 = (inv(Y_sol))*W82_sol.'; 
L83 = (inv(Y_sol))*W83_sol.'; 
L84 = (inv(Y_sol))*W84_sol.'; 
L85 = (inv(Y_sol))*W85_sol.'; 
L86 = (inv(Y_sol))*W86_sol.'; 
L87 = (inv(Y_sol))*W87_sol.'; 
L88 = (inv(Y_sol))*W88_sol.'; 
L89 = (inv(Y_sol))*W89_sol.'; 
L90 = (inv(Y_sol))*W90_sol.'; 
L91 = (inv(Y_sol))*W91_sol.'; 
L92 = (inv(Y_sol))*W92_sol.'; 
L93 = (inv(Y_sol))*W93_sol.'; 
L94 = (inv(Y_sol))*W94_sol.'; 
L95 = (inv(Y_sol))*W95_sol.'; 
L96 = (inv(Y_sol))*W96_sol.'; 
L97 = (inv(Y_sol))*W97_sol.'; 
L98 = (inv(Y_sol))*W98_sol.'; 
L99 = (inv(Y_sol))*W99_sol.'; 
L100 = (inv(Y_sol))*W100_sol.'; 
L101 = (inv(Y_sol))*W101_sol.'; 
L102 = (inv(Y_sol))*W102_sol.'; 
L103 = (inv(Y_sol))*W103_sol.'; 
L104 = (inv(Y_sol))*W104_sol.'; 
L105 = (inv(Y_sol))*W105_sol.'; 
L106 = (inv(Y_sol))*W106_sol.'; 
L107 = (inv(Y_sol))*W107_sol.'; 
L108 = (inv(Y_sol))*W108_sol.'; 
L109 = (inv(Y_sol))*W109_sol.'; 
L110 = (inv(Y_sol))*W110_sol.'; 
L111 = (inv(Y_sol))*W111_sol.'; 
L112 = (inv(Y_sol))*W112_sol.'; 
L113 = (inv(Y_sol))*W113_sol.'; 
L114 = (inv(Y_sol))*W114_sol.'; 
L115 = (inv(Y_sol))*W115_sol.'; 
L116 = (inv(Y_sol))*W116_sol.'; 
L117 = (inv(Y_sol))*W117_sol.'; 
L118 = (inv(Y_sol))*W118_sol.'; 
L119 = (inv(Y_sol))*W119_sol.'; 
L120 = (inv(Y_sol))*W120_sol.'; 
L121 = (inv(Y_sol))*W121_sol.'; 
L122 = (inv(Y_sol))*W122_sol.'; 
L123 = (inv(Y_sol))*W123_sol.'; 
L124 = (inv(Y_sol))*W124_sol.'; 
L125 = (inv(Y_sol))*W125_sol.'; 
L126 = (inv(Y_sol))*W126_sol.'; 
L127 = (inv(Y_sol))*W127_sol.'; 
L128 = (inv(Y_sol))*W128_sol.'; 
L129 = (inv(Y_sol))*W129_sol.'; 
L130 = (inv(Y_sol))*W130_sol.'; 
L131 = (inv(Y_sol))*W131_sol.'; 
L132 = (inv(Y_sol))*W132_sol.'; 
L133 = (inv(Y_sol))*W133_sol.'; 
L134 = (inv(Y_sol))*W134_sol.'; 
L135 = (inv(Y_sol))*W135_sol.'; 
L136 = (inv(Y_sol))*W136_sol.'; 
L137 = (inv(Y_sol))*W137_sol.'; 
L138 = (inv(Y_sol))*W138_sol.'; 
L139 = (inv(Y_sol))*W139_sol.'; 
L140 = (inv(Y_sol))*W140_sol.'; 
L141 = (inv(Y_sol))*W141_sol.'; 
L142 = (inv(Y_sol))*W142_sol.'; 
L143 = (inv(Y_sol))*W143_sol.'; 
L144 = (inv(Y_sol))*W144_sol.'; 
L145 = (inv(Y_sol))*W145_sol.'; 
L146 = (inv(Y_sol))*W146_sol.'; 
L147 = (inv(Y_sol))*W147_sol.'; 
L148 = (inv(Y_sol))*W148_sol.'; 
L149 = (inv(Y_sol))*W149_sol.'; 
L150 = (inv(Y_sol))*W150_sol.'; 
L151 = (inv(Y_sol))*W151_sol.'; 
L152 = (inv(Y_sol))*W152_sol.'; 
L153 = (inv(Y_sol))*W153_sol.'; 
L154 = (inv(Y_sol))*W154_sol.'; 
L155 = (inv(Y_sol))*W155_sol.'; 
L156 = (inv(Y_sol))*W156_sol.'; 
L157 = (inv(Y_sol))*W157_sol.'; 
L158 = (inv(Y_sol))*W158_sol.'; 
L159 = (inv(Y_sol))*W159_sol.'; 
L160 = (inv(Y_sol))*W160_sol.'; 
L161 = (inv(Y_sol))*W161_sol.'; 
L162 = (inv(Y_sol))*W162_sol.'; 
L163 = (inv(Y_sol))*W163_sol.'; 
L164 = (inv(Y_sol))*W164_sol.'; 
L165 = (inv(Y_sol))*W165_sol.'; 
L166 = (inv(Y_sol))*W166_sol.'; 
L167 = (inv(Y_sol))*W167_sol.'; 
L168 = (inv(Y_sol))*W168_sol.'; 
L169 = (inv(Y_sol))*W169_sol.'; 
L170 = (inv(Y_sol))*W170_sol.'; 
L171 = (inv(Y_sol))*W171_sol.'; 
L172 = (inv(Y_sol))*W172_sol.'; 
L173 = (inv(Y_sol))*W173_sol.'; 
L174 = (inv(Y_sol))*W174_sol.'; 
L175 = (inv(Y_sol))*W175_sol.'; 
L176 = (inv(Y_sol))*W176_sol.'; 
L177 = (inv(Y_sol))*W177_sol.'; 
L178 = (inv(Y_sol))*W178_sol.'; 
L179 = (inv(Y_sol))*W179_sol.'; 
L180 = (inv(Y_sol))*W180_sol.'; 
L181 = (inv(Y_sol))*W181_sol.'; 
L182 = (inv(Y_sol))*W182_sol.'; 
L183 = (inv(Y_sol))*W183_sol.'; 
L184 = (inv(Y_sol))*W184_sol.'; 
L185 = (inv(Y_sol))*W185_sol.'; 
L186 = (inv(Y_sol))*W186_sol.'; 
L187 = (inv(Y_sol))*W187_sol.'; 
L188 = (inv(Y_sol))*W188_sol.'; 
L189 = (inv(Y_sol))*W189_sol.'; 
L190 = (inv(Y_sol))*W190_sol.'; 
L191 = (inv(Y_sol))*W191_sol.'; 
L192 = (inv(Y_sol))*W192_sol.'; 
L193 = (inv(Y_sol))*W193_sol.'; 
L194 = (inv(Y_sol))*W194_sol.'; 
L195 = (inv(Y_sol))*W195_sol.'; 
L196 = (inv(Y_sol))*W196_sol.'; 
L197 = (inv(Y_sol))*W197_sol.'; 
L198 = (inv(Y_sol))*W198_sol.'; 
L199 = (inv(Y_sol))*W199_sol.'; 
L200 = (inv(Y_sol))*W200_sol.'; 
L201 = (inv(Y_sol))*W201_sol.'; 
L202 = (inv(Y_sol))*W202_sol.'; 
L203 = (inv(Y_sol))*W203_sol.'; 
L204 = (inv(Y_sol))*W204_sol.'; 
L205 = (inv(Y_sol))*W205_sol.'; 
L206 = (inv(Y_sol))*W206_sol.'; 
L207 = (inv(Y_sol))*W207_sol.'; 
L208 = (inv(Y_sol))*W208_sol.'; 
L209 = (inv(Y_sol))*W209_sol.'; 
L210 = (inv(Y_sol))*W210_sol.'; 
L211 = (inv(Y_sol))*W211_sol.'; 
L212 = (inv(Y_sol))*W212_sol.'; 
L213 = (inv(Y_sol))*W213_sol.'; 
L214 = (inv(Y_sol))*W214_sol.'; 
L215 = (inv(Y_sol))*W215_sol.'; 
L216 = (inv(Y_sol))*W216_sol.'; 
L217 = (inv(Y_sol))*W217_sol.'; 
L218 = (inv(Y_sol))*W218_sol.'; 
L219 = (inv(Y_sol))*W219_sol.'; 
L220 = (inv(Y_sol))*W220_sol.'; 
L221 = (inv(Y_sol))*W221_sol.'; 
L222 = (inv(Y_sol))*W222_sol.'; 
L223 = (inv(Y_sol))*W223_sol.'; 
L224 = (inv(Y_sol))*W224_sol.'; 
L225 = (inv(Y_sol))*W225_sol.'; 
L226 = (inv(Y_sol))*W226_sol.'; 
L227 = (inv(Y_sol))*W227_sol.'; 
L228 = (inv(Y_sol))*W228_sol.'; 
L229 = (inv(Y_sol))*W229_sol.'; 
L230 = (inv(Y_sol))*W230_sol.'; 
L231 = (inv(Y_sol))*W231_sol.'; 
L232 = (inv(Y_sol))*W232_sol.'; 
L233 = (inv(Y_sol))*W233_sol.'; 
L234 = (inv(Y_sol))*W234_sol.'; 
L235 = (inv(Y_sol))*W235_sol.'; 
L236 = (inv(Y_sol))*W236_sol.'; 
L237 = (inv(Y_sol))*W237_sol.'; 
L238 = (inv(Y_sol))*W238_sol.'; 
L239 = (inv(Y_sol))*W239_sol.'; 
L240 = (inv(Y_sol))*W240_sol.'; 
L241 = (inv(Y_sol))*W241_sol.'; 
L242 = (inv(Y_sol))*W242_sol.'; 
L243 = (inv(Y_sol))*W243_sol.'; 
L244 = (inv(Y_sol))*W244_sol.'; 
L245 = (inv(Y_sol))*W245_sol.'; 
L246 = (inv(Y_sol))*W246_sol.'; 
L247 = (inv(Y_sol))*W247_sol.'; 
L248 = (inv(Y_sol))*W248_sol.'; 
L249 = (inv(Y_sol))*W249_sol.'; 
L250 = (inv(Y_sol))*W250_sol.'; 
L251 = (inv(Y_sol))*W251_sol.'; 
L252 = (inv(Y_sol))*W252_sol.'; 
L253 = (inv(Y_sol))*W253_sol.'; 
L254 = (inv(Y_sol))*W254_sol.'; 
L255 = (inv(Y_sol))*W255_sol.'; 
L256 = (inv(Y_sol))*W256_sol.'; 



polos1_obs  = eig(Ad1 - L1*C);
polos2_obs  = eig(Ad2 - L2*C);
polos3_obs  = eig(Ad3 - L3*C);
polos4_obs  = eig(Ad4 - L4*C);
polos5_obs  = eig(Ad5 - L5*C);
polos6_obs  = eig(Ad6 - L6*C);
polos7_obs  = eig(Ad7 - L7*C);
polos8_obs  = eig(Ad8 - L8*C);
polos9_obs  = eig(Ad9 - L9*C);
polos10_obs = eig(Ad10 - L10*C);
polos11_obs = eig(Ad11 - L11*C);
polos12_obs = eig(Ad12 - L12*C);
polos13_obs = eig(Ad13 - L13*C);
polos14_obs = eig(Ad14 - L14*C);
polos15_obs = eig(Ad15 - L15*C);
polos16_obs = eig(Ad16 - L16*C);
polos17_obs  = eig(Ad17 - L17*C);
polos18_obs  = eig(Ad18 - L18*C);
polos19_obs  = eig(Ad19 - L19*C);
polos20_obs  = eig(Ad20 - L20*C);
polos21_obs  = eig(Ad21 - L21*C);
polos22_obs  = eig(Ad22 - L22*C);
polos23_obs  = eig(Ad23 - L23*C);
polos24_obs  = eig(Ad24 - L24*C);
polos25_obs  = eig(Ad25 - L25*C);
polos26_obs  = eig(Ad26 - L26*C);
polos27_obs  = eig(Ad27 - L27*C);
polos28_obs  = eig(Ad28 - L28*C);
polos29_obs  = eig(Ad29 - L29*C);
polos30_obs  = eig(Ad30 - L30*C);
polos31_obs  = eig(Ad31 - L31*C);
polos32_obs  = eig(Ad32 - L32*C);
polos33_obs  = eig(Ad33 - L33*C);
polos34_obs  = eig(Ad34 - L34*C);
polos35_obs  = eig(Ad35 - L35*C);
polos36_obs  = eig(Ad36 - L36*C);
polos37_obs  = eig(Ad37 - L37*C);
polos38_obs  = eig(Ad38 - L38*C);
polos39_obs  = eig(Ad39 - L39*C);
polos40_obs  = eig(Ad40 - L40*C);
polos41_obs  = eig(Ad41 - L41*C);
polos42_obs  = eig(Ad42 - L42*C);
polos43_obs  = eig(Ad43 - L43*C);
polos44_obs  = eig(Ad44 - L44*C);
polos45_obs  = eig(Ad45 - L45*C);
polos46_obs  = eig(Ad46 - L46*C);
polos47_obs  = eig(Ad47 - L47*C);
polos48_obs  = eig(Ad48 - L48*C);
polos49_obs  = eig(Ad49 - L49*C);
polos50_obs  = eig(Ad50 - L50*C);
polos51_obs  = eig(Ad51 - L51*C);
polos52_obs  = eig(Ad52 - L52*C);
polos53_obs  = eig(Ad53 - L53*C);
polos54_obs  = eig(Ad54 - L54*C);
polos55_obs  = eig(Ad55 - L55*C);
polos56_obs  = eig(Ad56 - L56*C);
polos57_obs  = eig(Ad57 - L57*C);
polos58_obs  = eig(Ad58 - L58*C);
polos59_obs  = eig(Ad59 - L59*C);
polos60_obs  = eig(Ad60 - L60*C);
polos61_obs  = eig(Ad61 - L61*C);
polos62_obs  = eig(Ad62 - L62*C);
polos63_obs  = eig(Ad63 - L63*C);
polos64_obs  = eig(Ad64 - L64*C);
polos65_obs  = eig(Ad65 - L65*C);
polos66_obs  = eig(Ad66 - L66*C);
polos67_obs  = eig(Ad67 - L67*C);
polos68_obs  = eig(Ad68 - L68*C);
polos69_obs  = eig(Ad69 - L69*C);
polos70_obs  = eig(Ad70 - L70*C);
polos71_obs  = eig(Ad71 - L71*C);
polos72_obs  = eig(Ad72 - L72*C);
polos73_obs  = eig(Ad73 - L73*C);
polos74_obs  = eig(Ad74 - L74*C);
polos75_obs  = eig(Ad75 - L75*C);
polos76_obs  = eig(Ad76 - L76*C);
polos77_obs  = eig(Ad77 - L77*C);
polos78_obs  = eig(Ad78 - L78*C);
polos79_obs  = eig(Ad79 - L79*C);
polos80_obs  = eig(Ad80 - L80*C);
polos81_obs  = eig(Ad81 - L81*C);
polos82_obs  = eig(Ad82 - L82*C);
polos83_obs  = eig(Ad83 - L83*C);
polos84_obs  = eig(Ad84 - L84*C);
polos85_obs  = eig(Ad85 - L85*C);
polos86_obs  = eig(Ad86 - L86*C);
polos87_obs  = eig(Ad87 - L87*C);
polos88_obs  = eig(Ad88 - L88*C);
polos89_obs  = eig(Ad89 - L89*C);
polos90_obs  = eig(Ad90 - L90*C);
polos91_obs  = eig(Ad91 - L91*C);
polos92_obs  = eig(Ad92 - L92*C);
polos93_obs  = eig(Ad93 - L93*C);
polos94_obs  = eig(Ad94 - L94*C);
polos95_obs  = eig(Ad95 - L95*C);
polos96_obs  = eig(Ad96 - L96*C);
polos97_obs  = eig(Ad97 - L97*C);
polos98_obs  = eig(Ad98 - L98*C);
polos99_obs  = eig(Ad99 - L99*C);
polos100_obs = eig(Ad100 - L100*C);
polos101_obs = eig(Ad101 - L101*C);
polos102_obs = eig(Ad102 - L102*C);
polos103_obs = eig(Ad103 - L103*C);
polos104_obs = eig(Ad104 - L104*C);
polos105_obs = eig(Ad105 - L105*C);
polos106_obs = eig(Ad106 - L106*C);
polos107_obs = eig(Ad107 - L107*C);
polos108_obs = eig(Ad108 - L108*C);
polos109_obs = eig(Ad109 - L109*C);
polos110_obs = eig(Ad110 - L110*C);
polos111_obs = eig(Ad111 - L111*C);
polos112_obs = eig(Ad112 - L112*C);
polos113_obs = eig(Ad113 - L113*C);
polos114_obs = eig(Ad114 - L114*C);
polos115_obs = eig(Ad115 - L115*C);
polos116_obs = eig(Ad116 - L116*C);
polos117_obs = eig(Ad117 - L117*C);
polos118_obs = eig(Ad118 - L118*C);
polos119_obs = eig(Ad119 - L119*C);
polos120_obs = eig(Ad120 - L120*C);
polos121_obs = eig(Ad121 - L121*C);
polos122_obs = eig(Ad122 - L122*C);
polos123_obs = eig(Ad123 - L123*C);
polos124_obs = eig(Ad124 - L124*C);
polos125_obs = eig(Ad125 - L125*C);
polos126_obs = eig(Ad126 - L126*C);
polos127_obs = eig(Ad127 - L127*C);
polos128_obs = eig(Ad128 - L128*C);
polos129_obs = eig(Ad129 - L129*C);
polos130_obs = eig(Ad130 - L130*C);
polos131_obs = eig(Ad131 - L131*C);
polos132_obs = eig(Ad132 - L132*C);
polos133_obs = eig(Ad133 - L133*C);
polos134_obs = eig(Ad134 - L134*C);
polos135_obs = eig(Ad135 - L135*C);
polos136_obs = eig(Ad136 - L136*C);
polos137_obs = eig(Ad137 - L137*C);
polos138_obs = eig(Ad138 - L138*C);
polos139_obs = eig(Ad139 - L139*C);
polos140_obs = eig(Ad140 - L140*C);
polos141_obs = eig(Ad141 - L141*C);
polos142_obs = eig(Ad142 - L142*C);
polos143_obs = eig(Ad143 - L143*C);
polos144_obs = eig(Ad144 - L144*C);
polos145_obs = eig(Ad145 - L145*C);
polos146_obs = eig(Ad146 - L146*C);
polos147_obs = eig(Ad147 - L147*C);
polos148_obs = eig(Ad148 - L148*C);
polos149_obs = eig(Ad149 - L149*C);
polos150_obs = eig(Ad150 - L150*C);
polos151_obs = eig(Ad151 - L151*C);
polos152_obs = eig(Ad152 - L152*C);
polos153_obs = eig(Ad153 - L153*C);
polos154_obs = eig(Ad154 - L154*C);
polos155_obs = eig(Ad155 - L155*C);
polos156_obs = eig(Ad156 - L156*C);
polos157_obs = eig(Ad157 - L157*C);
polos158_obs = eig(Ad158 - L158*C);
polos159_obs = eig(Ad159 - L159*C);
polos160_obs = eig(Ad160 - L160*C);
polos161_obs = eig(Ad161 - L161*C);
polos162_obs = eig(Ad162 - L162*C);
polos163_obs = eig(Ad163 - L163*C);
polos164_obs = eig(Ad164 - L164*C);
polos165_obs = eig(Ad165 - L165*C);
polos166_obs = eig(Ad166 - L166*C);
polos167_obs = eig(Ad167 - L167*C);
polos168_obs = eig(Ad168 - L168*C);
polos169_obs = eig(Ad169 - L169*C);
polos170_obs = eig(Ad170 - L170*C);
polos171_obs = eig(Ad171 - L171*C);
polos172_obs = eig(Ad172 - L172*C);
polos173_obs = eig(Ad173 - L173*C);
polos174_obs = eig(Ad174 - L174*C);
polos175_obs = eig(Ad175 - L175*C);
polos176_obs = eig(Ad176 - L176*C);
polos177_obs = eig(Ad177 - L177*C);
polos178_obs = eig(Ad178 - L178*C);
polos179_obs = eig(Ad179 - L179*C);
polos180_obs = eig(Ad180 - L180*C);
polos181_obs = eig(Ad181 - L181*C);
polos182_obs = eig(Ad182 - L182*C);
polos183_obs = eig(Ad183 - L183*C);
polos184_obs = eig(Ad184 - L184*C);
polos185_obs = eig(Ad185 - L185*C);
polos186_obs = eig(Ad186 - L186*C);
polos187_obs = eig(Ad187 - L187*C);
polos188_obs = eig(Ad188 - L188*C);
polos189_obs = eig(Ad189 - L189*C);
polos190_obs = eig(Ad190 - L190*C);
polos191_obs = eig(Ad191 - L191*C);
polos192_obs = eig(Ad192 - L192*C);
polos193_obs = eig(Ad193 - L193*C);
polos194_obs = eig(Ad194 - L194*C);
polos195_obs = eig(Ad195 - L195*C);
polos196_obs = eig(Ad196 - L196*C);
polos197_obs = eig(Ad197 - L197*C);
polos198_obs = eig(Ad198 - L198*C);
polos199_obs = eig(Ad199 - L199*C);
polos200_obs = eig(Ad200 - L200*C);
polos201_obs = eig(Ad201 - L201*C);
polos202_obs = eig(Ad202 - L202*C);
polos203_obs = eig(Ad203 - L203*C);
polos204_obs = eig(Ad204 - L204*C);
polos205_obs = eig(Ad205 - L205*C);
polos206_obs = eig(Ad206 - L206*C);
polos207_obs = eig(Ad207 - L207*C);
polos208_obs = eig(Ad208 - L208*C);
polos209_obs = eig(Ad209 - L209*C);
polos210_obs = eig(Ad210 - L210*C);
polos211_obs = eig(Ad211 - L211*C);
polos212_obs = eig(Ad212 - L212*C);
polos213_obs = eig(Ad213 - L213*C);
polos214_obs = eig(Ad214 - L214*C);
polos215_obs = eig(Ad215 - L215*C);
polos216_obs = eig(Ad216 - L216*C);
polos217_obs = eig(Ad217 - L217*C);
polos218_obs = eig(Ad218 - L218*C);
polos219_obs = eig(Ad219 - L219*C);
polos220_obs = eig(Ad220 - L220*C);
polos221_obs = eig(Ad221 - L221*C);
polos222_obs = eig(Ad222 - L222*C);
polos223_obs = eig(Ad223 - L223*C);
polos224_obs = eig(Ad224 - L224*C);
polos225_obs = eig(Ad225 - L225*C);
polos226_obs = eig(Ad226 - L226*C);
polos227_obs = eig(Ad227 - L227*C);
polos228_obs = eig(Ad228 - L228*C);
polos229_obs = eig(Ad229 - L229*C);
polos230_obs = eig(Ad230 - L230*C);
polos231_obs = eig(Ad231 - L231*C);
polos232_obs = eig(Ad232 - L232*C);
polos233_obs = eig(Ad233 - L233*C);
polos234_obs = eig(Ad234 - L234*C);
polos235_obs = eig(Ad235 - L235*C);
polos236_obs = eig(Ad236 - L236*C);
polos237_obs = eig(Ad237 - L237*C);
polos238_obs = eig(Ad238 - L238*C);
polos239_obs = eig(Ad239 - L239*C);
polos240_obs = eig(Ad240 - L240*C);
polos241_obs = eig(Ad241 - L241*C);
polos242_obs = eig(Ad242 - L242*C);
polos243_obs = eig(Ad243 - L243*C);
polos244_obs = eig(Ad244 - L244*C);
polos245_obs = eig(Ad245 - L245*C);
polos246_obs = eig(Ad246 - L246*C);
polos247_obs = eig(Ad247 - L247*C);
polos248_obs = eig(Ad248 - L248*C);
polos249_obs = eig(Ad249 - L249*C);
polos250_obs = eig(Ad250 - L250*C);
polos251_obs = eig(Ad251 - L251*C);
polos252_obs = eig(Ad252 - L252*C);
polos253_obs = eig(Ad253 - L253*C);
polos254_obs = eig(Ad254 - L254*C);
polos255_obs = eig(Ad255 - L255*C);
polos256_obs = eig(Ad256 - L256*C);

% Vector de valores propios (polos)
% Graficar el círculo unitario
theta = linspace(0, 2*pi, 100);
x_circulo = cos(theta); % Parte real
y_circulo = sin(theta); % Parte imaginaria

figure;
hold on;
plot(x_circulo, y_circulo, 'r--', 'LineWidth', 2); % Círculo unitario
for i = 1:256
    poloeval = eval(['polos' num2str(i) '_obs'])
    plot(real(poloeval), ...
         imag(poloeval), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Valores propios
end
hold off

% Configuración del gráfico
axis equal;
xlabel('Parte real');
ylabel('Parte imaginaria');
title('Valores propios y el círculo unitario');
legend('Círculo unitario', 'Valores propios');
grid on;

% %% CHECO POLOS
% for i = 1:256
%     % Acceder a las matrices L y W usando eval para construir los nombres de las variables
%     Lsol = eval(['L' num2str(i)]);  % Accede a L1FUGA2, L2FUGA2, ..., L64FUGA2
%     % Calcular los polos del observador
%     Polos_Obs{i,1} = eig(adis{i,1} - Lsol * Cd);
% end
% %% CHECO QUE ESTEN DENTRO DEL CIRCULO UNITARIO
% % for i=1:1:64
%     CeldaPolo = abs(real(Polos_Obs{i,1}))
%     tic
% 
%     for j=1:1:7
%         if CeldaPolo<=1
%             fprintf('Polo dentro del círculo unitario \n');
%         else
%             fprintf('Hay un polo fuera del círculo unitario en la iteración %d, posición %d \n', i,j)
%         end
%     end
%     tic
% 
% end


tic
n=1
enclave=0;

errorcuadraticoz3=0;
MAE = 0;
while(n<=iteraciones)
tic
    Qin = QinFiltrado(n,1);
    Qout = QoutFiltrado(n,1);
    Hin = HinFiltrado(n,1);
    Hout = HoutFiltrado(n,1);

    if (Qin-Qout>=10e-04)  && enclave ==0 && n>=120000
        enclave =1;
        tic
    end

    if enclave ==1
        Q1max=max(QinFiltrado)+0.0001;
        Q1min=min(QinFiltrado)-0.0001;

        H2max=19.5;
        H2min=11.5;


        Q2max=(Q1max-(lam1*sqrt(H2max)))+0.0002;
        Q2min= Q1min-(lam1*sqrt(H2min))-0.0001;



        H3max=16.5;
        H3min=10.5;


        Q3max = (Q2max-(lam2*sqrt(H3max)))+0.0002;
        Q3min = Q2min-(lam2*sqrt(H3min))-0.0001;


        H4max = 15.5;
        H4min = 9.5;

        Q4max = max(QoutFiltrado)+0.0001;
        Q4min = min(QoutFiltrado)-0.0001;

        z3max = 80;
        z3min = 60;

        Q1Sch = (Q1max - x1) / (Q1max - Q1min);
        H2Sch = (H2max - x2) / (H2max - H2min);
        Q2Sch = (Q2max - x3) / (Q2max - Q2min);
        H3Sch = (H3max - x4) / (H3max - H3min);
        Q3Sch = (Q3max - x5) / (Q3max - Q3min);
        H4Sch = (H4max - x6) / (H4max - H4min);
        Q4Sch = (Q4max - x7) / (Q4max - Q4min);
        z3Sch = (z3max - x8) / (z3max - z3min);
        tic

        if n>=50073 
            tic
        end
        if Q1Sch>1 || Q1Sch <=0 || H2Sch >1 || H2Sch <=0 || Q2Sch>1 || Q2Sch<=0 || H3Sch>1 || H3Sch<=0 ||Q3Sch>1 || Q3Sch<=0 || n>=50014
            tic
        end
        Q1 = x1;
        H2 = x2;
        Q2 = x3;
        H3 = x4;
        Q3 = x5;
        H4 = x6;
        Q4 = x7;
        z3 = x8;
        lam = x9;
        Hk = [1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0];

        A = Ar;
        promediolongitud = Leq;
                    promediolongitud = Leq;
            Temp_av = TemperaturaFiltrada(n);
            Viscosidad_av = interp1(Temperatura,  v,Temp_av, 'spline','extrap'); % Interpolación si es necesario
            Comprensibilidad_av = interp1(Temperatura, Comprensibilidad, Temp_av, 'spline','extrap'); % Interpolación si es necesario
            Densidad_av = interp1(Temperatura1, Densidad1, Temp_av, 'spline','extrap'); % Interpolación si es necesario
            MEA_av = interp1(Temperatura,  MEA,Temp_av, 'spline','extrap'); % Interpolación si es necesario
            E=E;
            ro = Densidad_av;
            D = D;
            e=e;
            K=MEA_av;
            Rein=((Qin/A)*D)/(Viscosidad_av);
            Reout=((Qout/A)*D)/(Viscosidad_av);
% % 
% % 
% % 
% % % f1=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Rein)))^-1)^2;
 f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
% % % f2=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Reout)))^-1)^2;
 f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;
miu1=f1/(2*D*A);
miu2=f2/(2*D*A);
miu3=f2/(2*D*A);
miu4=f2/(2*D*A);

% Varvis1(n)=miu1;
% Varvis2(n)=miu2;
 %E=(34.247*((Temp_av-55)/30.277)^5 + 94.874*((Temp_av-55)/30.277)^4 - 323.75*((Temp_av-55)/30.277)^3 + 799.98*((Temp_av-55)/30.277)^2 - 2549.7*((Temp_av-55)/30.277) + 3368.6)*g*1e4;

            b=sqrt((K/ro)/(1+K*D/(e*E)));
        dz1 = (z1);
        dz2 = (z2-z1);
        dz3 = (z3-z2);
        dz4 = (L-z3);
        lambda1 = lam1;
        lambda2 = lam2;

        %% HAGO MIS Q ESTIMADAS (Qi gorro j+1)
        Q11=Q1+dt*((-((g*A)*(H2-Hin))/(dz1))-miu1*Q1^2);
        Q21=Q2+dt*((-((g*A)*(H3-H2))/(dz2))-miu2*Q2^2);
        Q31=Q3+dt*((-((g*A)*(H4-H3))/(dz3))-miu3*Q3^2);
        Q41=Q4+dt*((-((g*A)*(Hout-H4))/(dz4))-miu4*Q4^2);

        %% HAGO MIS H ESTIMADAS (Hi gorro j+1)
        H21=H2+dt*(-(b^2)*(Q2-Q1+lam1*sqrt(abs(H2)))/(g*A*dz1));
        H31=H3+dt*(-(b^2)*(Q3-Q2+lam2*sqrt(abs(H3)))/(g*A*dz2));
        H41=H4+dt*(-(b^2)*(Q4-Q3+lam*sqrt(abs(H4)))/(g*A*dz3));

        %%Construyo la primera linea
        phi1=-miu1*Q11;
        phi2=((g*Ar)*(H21-1))/(z1);
        phi3=-(g*Ar*H21^2)/(z1*z3);
        tic
        
        %%Construyo la segunda linea
        phi4=(b^2/(g*Ar*(z1)))*(1 - (lambda1*sqrt(abs(H21)))/(2*Q11));
        phi5=(-b^2/(g*Ar*(z1)))*(1 + (lambda1*sqrt(abs(H21)))/(2*Q21));
        tic
        %%Construyo la tercera linea
        phi6=(g*Ar)/(z2-z1);
        phi7=-miu2*Q21;
        phi8=-(phi6);
        tic
        %%construyo la cuarta linea
        phi9=(b^2/(g*Ar*(z2-z1)))*(1 - (lambda2*sqrt(abs(H31)))/(2*Q21));
        phi10=(-b^2/(g*Ar*(z2-z1)))*(1 + (lambda2*sqrt(abs(H31)))/(2*Q31));
        tic
        %%Construyo la quinta linea
        phi11=(g*Ar)/(z3-z2);
        phi12=-miu3*Q31;
        phi13=-phi11;
        tic
        
        %%Construyo la sexta linea
        phi14=b^2/(g*Ar*(z3-z2));
        phi15=-phi14;
        phi16=-phi14*sqrt(abs(H41));
        
        %%Construyo la septima linea
        phi17=(g*Ar)/(L-z3);
        phi18=-(miu4*Q41);
        


        Aact1 = [phi1 phi2 0 0 0 0 0 phi3 0;
        phi4 0 phi5 0 0 0 0 0 0;
        0 phi6 phi7 phi8 0 0 0 0 0;
        0 0 phi9 0 phi10 0 0 0 0;
        0 0 0 phi11 phi12 phi13 0 0 0;
        0 0 0 0 phi14 0 phi15 0 phi16;
        0 0 0 0 0 phi17 phi18 0 0;
        zeros(2,9)];

        xact1 = [Q11, H21, Q21, H31, Q31, H41, Q41, z3, lam]' * (dt / 2);
        Bact1=[(g*Ar)/z1, 0;
        0, 0;
        0, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z3);
        0, 0;
        0, 0];
        u = [Hin * (dt / 2); Hout * (dt / 2)];
        Ax1 = Aact1 * xact1;
        Bu1 = Bact1 * u;

        %%Construyo la primera linea
        phi1=-miu1*Q1;
        phi2=((g*Ar)*(H2-1))/(z1);
        phi3=-(g*Ar*H2^2)/(z1*z3);
        tic
        
        %%Construyo la segunda linea
        phi4=(b^2/(g*Ar*(z1)))*(1 - (lambda1*sqrt(abs(H2)))/(2*Q1));
        phi5=(-b^2/(g*Ar*(z1)))*(1 + (lambda1*sqrt(abs(H2)))/(2*Q2));
        tic
        %%Construyo la tercera linea
        phi6=(g*Ar)/(z2-z1);
        phi7=-miu2*Q2;
        phi8=-(phi6);
        tic
        %%construyo la cuarta linea
        phi9=(b^2/(g*Ar*(z2-z1)))*(1 - (lambda2*sqrt(abs(H3)))/(2*Q2));
        phi10=(-b^2/(g*Ar*(z2-z1)))*(1 + (lambda2*sqrt(abs(H3)))/(2*Q3));
        tic
        %%Construyo la quinta linea
        phi11=(g*Ar)/(z3-z2);
        phi12=-miu3*Q3;
        phi13=-phi11;
        tic
        
        %%Construyo la sexta linea
        phi14=b^2/(g*Ar*(z3-z2));
        phi15=-phi14;
        phi16=-phi14*sqrt(abs(H4));
        
        %%Construyo la septima linea
        phi17=(g*Ar)/(L-z3);
        phi18=-(miu4*Q4);

        Aact2 = [phi1 phi2 0 0 0 0 0 phi3 0;
        phi4 0 phi5 0 0 0 0 0 0;
        0 phi6 phi7 phi8 0 0 0 0 0;
        0 0 phi9 0 phi10 0 0 0 0;
        0 0 0 phi11 phi12 phi13 0 0 0;
        0 0 0 0 phi14 0 phi15 0 phi16;
        0 0 0 0 0 phi17 phi18 0 0;
        zeros(2,9)];

        xact2 = [Q1, H2, Q2, H3, Q3, H4, Q4, z3, lam]' * (dt / 2);
        Bact2=[(g*Ar)/z1, 0;
        0, 0;
        0, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z3);
        0, 0;
        0, 0];
        u = [Hin * (dt / 2); Hout * (dt / 2)];
        Ax2 = Aact2 * xact2;
        Bu2 = Bact2 * u;

        xx = [Q1; H2; Q2; H3; Q3; H4; Q4; z3; lam];
        dXmat = (Ax1 + Ax2) + (Bu1 + Bu2) + xx;
        dQ1 = dXmat(1);
        dH2 = dXmat(2);
        dQ2 = dXmat(3);
        dH3 = dXmat(4);
        dQ3 = dXmat(5);
        dH4 = dXmat(6);
        dQ4 = dXmat(7);
        dXm = [dQ1; dH2; dQ2; dH3; dQ3; dH4; dQ4; z3; lam];

        %% Cálculo de los coeficientes de pertenencia `mu`
        Q1SchMax= (1 - Q1Sch);
        Q1SchMin= Q1Sch;

        H2SchMax= (1 - H2Sch);
        H2SchMin= H2Sch;

        Q2SchMax= (1 - Q2Sch);
        Q2SchMin= Q2Sch;

        H3SchMax= (1 - H3Sch);
        H3SchMin= H3Sch;

        Q3SchMax= (1 - Q3Sch);
        Q3SchMin= Q3Sch;

        H4SchMax = (1-H4Sch);
        H4SchMin = H4Sch;

        Q4SchMax = (1-Q4Sch);
        Q4SchMin = Q4Sch;

        z3SchMax= (1-z3Sch);
        z3SchMin = z3Sch;

        mu_vals = obtenermusF3(Q1SchMax, Q1SchMin, H2SchMax, H2SchMin, Q2SchMax, Q2SchMin,H3SchMax,H3SchMin,Q3SchMax,Q3SchMin,H4SchMax,H4SchMin,Q4SchMax,Q4SchMin, z3SchMax, z3SchMin, Combinatoria, 8);
        tic

        muOBS3(n) = sum(mu_vals);
        tic
        L_interp=0;
        L_interp = L1 * mu_vals(1) + ...
           L2 * mu_vals(2) + ...
           L3 * mu_vals(3) + ...
           L4 * mu_vals(4) + ...
           L5 * mu_vals(5) + ...
           L6 * mu_vals(6) + ...
           L7 * mu_vals(7) + ...
           L8 * mu_vals(8) + ...
           L9 * mu_vals(9) + ...
           L10 * mu_vals(10) + ...
           L11 * mu_vals(11) + ...
           L12 * mu_vals(12) + ...
           L13 * mu_vals(13) + ...
           L14 * mu_vals(14) + ...
           L15 * mu_vals(15) + ...
           L16 * mu_vals(16) + ...
           L17 * mu_vals(17) + ...
           L18 * mu_vals(18) + ...
           L19 * mu_vals(19) + ...
           L20 * mu_vals(20) + ...
           L21 * mu_vals(21) + ...
           L22 * mu_vals(22) + ...
           L23 * mu_vals(23) + ...
           L24 * mu_vals(24) + ...
           L25 * mu_vals(25) + ...
           L26 * mu_vals(26) + ...
           L27 * mu_vals(27) + ...
           L28 * mu_vals(28) + ...
           L29 * mu_vals(29) + ...
           L30 * mu_vals(30) + ...
           L31 * mu_vals(31) + ...
           L32 * mu_vals(32) + ...
           L33 * mu_vals(33) + ...
           L34 * mu_vals(34) + ...
           L35 * mu_vals(35) + ...
           L36 * mu_vals(36) + ...
           L37 * mu_vals(37) + ...
           L38 * mu_vals(38) + ...
           L39 * mu_vals(39) + ...
           L40 * mu_vals(40) + ...
           L41 * mu_vals(41) + ...
           L42 * mu_vals(42) + ...
           L43 * mu_vals(43) + ...
           L44 * mu_vals(44) + ...
           L45 * mu_vals(45) + ...
           L46 * mu_vals(46) + ...
           L47 * mu_vals(47) + ...
           L48 * mu_vals(48) + ...
           L49 * mu_vals(49) + ...
           L50 * mu_vals(50) + ...
           L51 * mu_vals(51) + ...
           L52 * mu_vals(52) + ...
           L53 * mu_vals(53) + ...
           L54 * mu_vals(54) + ...
           L55 * mu_vals(55) + ...
           L56 * mu_vals(56) + ...
           L57 * mu_vals(57) + ...
           L58 * mu_vals(58) + ...
           L59 * mu_vals(59) + ...
           L60 * mu_vals(60) + ...
           L61 * mu_vals(61) + ...
            L62  * mu_vals(62)  + ...
            L63  * mu_vals(63)  + ...
            L64  * mu_vals(64)  + ...
            L65  * mu_vals(65)  + ...
            L66  * mu_vals(66)  + ...
            L67  * mu_vals(67)  + ...
            L68  * mu_vals(68)  + ...
            L69  * mu_vals(69)  + ...
            L70  * mu_vals(70)  + ...
            L71  * mu_vals(71)  + ...
            L72  * mu_vals(72)  + ...
            L73  * mu_vals(73)  + ...
            L74  * mu_vals(74)  + ...
            L75  * mu_vals(75)  + ...
            L76  * mu_vals(76)  + ...
            L77  * mu_vals(77)  + ...
            L78  * mu_vals(78)  + ...
            L79  * mu_vals(79)  + ...
            L80  * mu_vals(80)  + ...
            L81  * mu_vals(81)  + ...
            L82  * mu_vals(82)  + ...
            L83  * mu_vals(83)  + ...
            L84  * mu_vals(84)  + ...
            L85  * mu_vals(85)  + ...
            L86  * mu_vals(86)  + ...
            L87  * mu_vals(87)  + ...
            L88  * mu_vals(88)  + ...
            L89  * mu_vals(89)  + ...
            L90  * mu_vals(90)  + ...
            L91  * mu_vals(91)  + ...
            L92  * mu_vals(92)  + ...
            L93  * mu_vals(93)  + ...
            L94  * mu_vals(94)  + ...
            L95  * mu_vals(95)  + ...
            L96  * mu_vals(96)  + ...
            L97  * mu_vals(97)  + ...
            L98  * mu_vals(98)  + ...
            L99  * mu_vals(99)  + ...
            L100 * mu_vals(100) + ...
            L101 * mu_vals(101) + ...
            L102 * mu_vals(102) + ...
            L103 * mu_vals(103) + ...
            L104 * mu_vals(104) + ...
            L105 * mu_vals(105) + ...
            L106 * mu_vals(106) + ...
            L107 * mu_vals(107) + ...
            L108 * mu_vals(108) + ...
            L109 * mu_vals(109) + ...
            L110 * mu_vals(110) + ...
            L111 * mu_vals(111) + ...
            L112 * mu_vals(112) + ...
            L113 * mu_vals(113) + ...
            L114 * mu_vals(114) + ...
            L115 * mu_vals(115) + ...
            L116 * mu_vals(116) + ...
            L117 * mu_vals(117) + ...
            L118 * mu_vals(118) + ...
            L119 * mu_vals(119) + ...
            L120 * mu_vals(120) + ...
            L121 * mu_vals(121) + ...
            L122 * mu_vals(122) + ...
            L123 * mu_vals(123) + ...
            L124 * mu_vals(124) + ...
            L125 * mu_vals(125) + ...
            L126 * mu_vals(126) + ...
            L127 * mu_vals(127) + ...
            L128 * mu_vals(128) + ...
            L129 * mu_vals(129) + ...
            L130 * mu_vals(130) + ...
            L131 * mu_vals(131) + ...
            L132 * mu_vals(132) + ...
            L133 * mu_vals(133) + ...
            L134 * mu_vals(134) + ...
            L135 * mu_vals(135) + ...
            L136 * mu_vals(136) + ...
            L137 * mu_vals(137) + ...
            L138 * mu_vals(138) + ...
            L139 * mu_vals(139) + ...
            L140 * mu_vals(140) + ...
            L141 * mu_vals(141) + ...
            L142 * mu_vals(142) + ...
            L143 * mu_vals(143) + ...
            L144 * mu_vals(144) + ...
            L145 * mu_vals(145) + ...
            L146 * mu_vals(146) + ...
            L147 * mu_vals(147) + ...
            L148 * mu_vals(148) + ...
            L149 * mu_vals(149) + ...
            L150 * mu_vals(150) + ...
            L151 * mu_vals(151) + ...
            L152 * mu_vals(152) + ...
            L153 * mu_vals(153) + ...
            L154 * mu_vals(154) + ...
            L155 * mu_vals(155) + ...
            L156 * mu_vals(156) + ...
            L157 * mu_vals(157) + ...
            L158 * mu_vals(158) + ...
            L159 * mu_vals(159) + ...
            L160 * mu_vals(160) + ...
            L161 * mu_vals(161) + ...
            L162 * mu_vals(162) + ...
            L163 * mu_vals(163) + ...
            L164 * mu_vals(164) + ...
            L165 * mu_vals(165) + ...
            L166 * mu_vals(166) + ...
            L167 * mu_vals(167) + ...
            L168 * mu_vals(168) + ...
            L169 * mu_vals(169) + ...
            L170 * mu_vals(170) + ...
            L171 * mu_vals(171) + ...
            L172 * mu_vals(172) + ...
            L173 * mu_vals(173) + ...
            L174 * mu_vals(174) + ...
            L175 * mu_vals(175) + ...
            L176 * mu_vals(176) + ...
            L177 * mu_vals(177) + ...
            L178 * mu_vals(178) + ...
            L179 * mu_vals(179) + ...
            L180 * mu_vals(180) + ...
            L181 * mu_vals(181) + ...
            L182 * mu_vals(182) + ...
            L183 * mu_vals(183) + ...
            L184 * mu_vals(184) + ...
            L185 * mu_vals(185) + ...
            L186 * mu_vals(186) + ...
            L187 * mu_vals(187) + ...
            L188 * mu_vals(188) + ...
            L189 * mu_vals(189) + ...
            L190 * mu_vals(190) + ...
            L191 * mu_vals(191) + ...
            L192 * mu_vals(192) + ...
            L193 * mu_vals(193) + ...
            L194 * mu_vals(194) + ...
            L195 * mu_vals(195) + ...
            L196 * mu_vals(196) + ...
            L197 * mu_vals(197) + ...
            L198 * mu_vals(198) + ...
            L199 * mu_vals(199) + ...
            L200 * mu_vals(200) + ...
            L201 * mu_vals(201) + ...
            L202 * mu_vals(202) + ...
            L203 * mu_vals(203) + ...
            L204 * mu_vals(204) + ...
            L205 * mu_vals(205) + ...
            L206 * mu_vals(206) + ...
            L207 * mu_vals(207) + ...
            L208 * mu_vals(208) + ...
            L209 * mu_vals(209) + ...
            L210 * mu_vals(210) + ...
            L211 * mu_vals(211) + ...
            L212 * mu_vals(212) + ...
            L213 * mu_vals(213) + ...
            L214 * mu_vals(214) + ...
            L215 * mu_vals(215) + ...
            L216 * mu_vals(216) + ...
            L217 * mu_vals(217) + ...
            L218 * mu_vals(218) + ...
            L219 * mu_vals(219) + ...
            L220 * mu_vals(220) + ...
            L221 * mu_vals(221) + ...
            L222 * mu_vals(222) + ...
            L223 * mu_vals(223) + ...
            L224 * mu_vals(224) + ...
            L225 * mu_vals(225) + ...
            L226 * mu_vals(226) + ...
            L227 * mu_vals(227) + ...
            L228 * mu_vals(228) + ...
            L229 * mu_vals(229) + ...
            L230 * mu_vals(230) + ...
            L231 * mu_vals(231) + ...
            L232 * mu_vals(232) + ...
            L233 * mu_vals(233) + ...
            L234 * mu_vals(234) + ...
            L235 * mu_vals(235) + ...
            L236 * mu_vals(236) + ...
            L237 * mu_vals(237) + ...
            L238 * mu_vals(238) + ...
            L239 * mu_vals(239) + ...
            L240 * mu_vals(240) + ...
            L241 * mu_vals(241) + ...
            L242 * mu_vals(242) + ...
            L243 * mu_vals(243) + ...
            L244 * mu_vals(244) + ...
            L245 * mu_vals(245) + ...
            L246 * mu_vals(246) + ...
            L247 * mu_vals(247) + ...
            L248 * mu_vals(248) + ...
            L249 * mu_vals(249) + ...
            L250 * mu_vals(250) + ...
            L251 * mu_vals(251) + ...
            L252 * mu_vals(252) + ...
            L253 * mu_vals(253) + ...
            L254 * mu_vals(254) + ...
            L255 * mu_vals(255) + ...
            L256 * mu_vals(256);


        Kk = L_interp;
        tic
        dX = dXm + Kk * ([Qin; Qout] - [dQ1; dQ4]);
        x1 = dX(1);
        x2 = dX(2);
        x3 = dX(3);
        x4 = dX(4);
        x5 = dX(5);
        x6 = dX(6);
        x7 = dX(7);
        x8 = dX(8);
        x9 = dX(9);

        tic
    else
        x1 = x1i;
        x2 = x2i;
        x3 = x3i;
        x4 = x4i;
        x5 = x5i;
        x6 = x6i;
        x7 = x7i;
        x8 = x8i;
        x9 = x9i;

        % Inicialización de vectores
        dX = [x1; x2; x3; x4; x5; x6; x7; x8; x9];
        mus = 0;
        mu_vals = zeros(256, 1);  % Inicializamos todos los `mu_i` a 0
    end


for l = 1:256
    mu_Leak2{l}(n) = mu_vals(l); % Almacena el valor
end

    tic

    Q1Fuga3Observador(n)=dX(1);
    H2Fuga3Observador(n)=dX(2);
    Q2Fuga3Observador(n)=dX(3);
    H3Fuga3Observador(n)=dX(4);
    Q3Fuga3Observador(n)=dX(5);
    H4Fuga3Observador(n)=dX(6);
    Q4Fuga3Observador(n)=dX(7);
    z3Fuga3Observador(n)=dX(8);
    lambda3Fuga3Observador(n)=dX(9);
            errorz3(n) = abs(49.9-(z3Fuga3Observador(n)*(Lr/Leq)));
            if enclave == 1
            errorcuadraticoz3 = errorcuadraticoz3 + (49.9-(z3Fuga3Observador(n)*(Lr/Leq)))^2;
            MAE = MAE + abs(49.9-(z3Fuga3Observador(n)*(Lr/Leq)));

            end
    n=n+1;
end
tic
RMSEz3 = sqrt(errorcuadraticoz3 / (190000-130150));
MAEz3 = 1/(190000-130150) * MAE;

tic
Ts = 0.01; % tiempo de muestreo
t = 0:Ts:(length(QinFiltrado)-1)*Ts;

H2prom = sum(H2Fuga1Observador(1,(12000:end)))/length(H2Fuga1Observador(1,(12000:end)));
H3prom = sum(H3Fuga2Observador(1,(70300:end)))/length(H3Fuga2Observador(1,(70300:end)));

Qleak1 = lamda1Fuga1Observador.*sqrt(H2Fuga1Observador);
Qleak1Total = Qleak1;

Qleak2 = lamda2Fuga2Observador.*sqrt(H3Fuga2Observador);
Qleak2Total = Qleak2(70001:130100)+lam1*sqrt(H2prom);


Qleak3 = lambda3Fuga3Observador.*sqrt(H4Fuga3Observador);
Qleak3Total = Qleak3(130101:end)+lam1*sqrt(H2prom)+lam2*sqrt(H3prom);

QLeakTotalObs = [Qleak1Total, Qleak2Total,Qleak3Total];
%QLeakTotalObs = QLeakTotalObs(1:190000);




% %% GRAFICO LA EVOLUCION DE MU'S
% set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
% figure('Name','MUS EVOLUTION    ','NumberTitle','off','color', [1 1 1])
% axes1 = axes('FontSize',20);
% %% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[1250 1900]); 
% ylim(axes1,[-0.001 0.05]); 
% box(axes1,'on');
% grid(axes1,'on');
% hold(axes1,'all');
% % Genera colores RGB
% colors = lines(256);  % Genera 64 colores distintos
% % Plotea cada mu_Leak2 en un color diferente usando un ciclo for
% for l = 1:256
%     % Plotea cada mu_Leak2 con su color correspondiente
%     plot(t, mu_Leak2{l}, 'Color', colors(l, :)); hold on;
%     tic
% end
% grid on;
% xlabel('(s)','FontName','Times New Roman','FontSize', 20);
% print(gcf, fullfile(output_folder, 'musEvolution.eps'), '-depsc');
% tic

%% GRAFICO LA EVOLUCION DE MU'S (256 curvas con tamaño fijo)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

figure('Name','MUS EVOLUTION','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);

% Límites de ejes
xlim(axes1,[1250 1900]); 
ylim(axes1,[-0.001 0.07]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Colores
colors = lines(256);

% Graficar cada curva
for l = 1:256
    y = mu_Leak2{l}(:)'; % asegurar vector fila
    if length(y) ~= length(t)
        warning('Curva %d con longitud diferente: len(y)=%d, len(t)=%d', ...
                l, length(y), length(t));
        continue
    end
    plot(t, y, 'Color', colors(l,:), 'LineWidth', 0.5); 
    hold on;
end

% Etiquetas
xlabel('(s)','FontName','Times New Roman','FontSize', 20);

% ==== Ajuste de tamaño fijo ====
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 20, 15]); % [x y ancho alto]
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0, 0, 20, 15]);

% Exportar EPS con mismo tamaño
%print(gcf, fullfile(output_folder, 'musEvolutionL3.eps'), '-depsc','-vector');




%% GRÁFICA DE MUS, LEAK 1, LEAK2, LEAK3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','SUMAMUS','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Configuración de los límites y ticks del eje X
xlim(axes1,[0 1900]); 
ylim(axes1,[-0.1 1.1]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

musL1 = [musL1,ones(1,120000)];
mu = [mu, ones(1,59900)];

% Colores en tonos de gris
plot(t, musL1, 'color', [0, 0.45, 0.74], 'linewidth', 1.5), hold on  % Gris oscuro
plot(t, mu, 'color', [0.85, 0.33, 0.10], 'linewidth', 1.5), hold on  % Gris intermedio
plot(t, muOBS3, 'color', [0.47, 0.67, 0.19], 'linewidth', 1.5), hold on  % Gris claro

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$\sum \psi_{Leak_1}(t)$';
leyenda2 = '$\sum \psi_{Leak_2}(t)$';
%leyenda = legend(leyenda1,leyenda2,'Location','southeast');
%print(gcf,fullfile(output_folder, 'MUS.eps'), '-depsc');

tic


%% GRÁFICA DE HIN Y HOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Hin&Hout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1900]); 
ylim(axes1,[8 22]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,HinFiltrado, 'color', [0, 0.45, 0.74],'linewidth',0.5), hold on
plot(t,HoutFiltrado,'color', [0.85, 0.33, 0.10] ,'linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$H_{in}(t)$';
leyenda2 = '$H_{out}(t)$';
%print(gcf, fullfile(output_folder, 'presionesEntradaYSalidaLPV.eps'), '-depsc');



%% GRÁFICA DE QIN Y QOUT
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Qin&Qout','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1900]); 
ylim(axes1,[0.0083 0.01]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
plot(t,QinFiltrado, 'color', [0, 0.45, 0.74],'linewidth',0.5), hold on
plot(t,QoutFiltrado,'color', [0.85, 0.33, 0.10],'linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^3/s)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$Q_{in}(t)$';
leyenda2 = '$Q_{out}(t)$';
%print(gcf,fullfile(output_folder, 'FlujoEntradaYSalida.eps'), '-depsc');


%% GRÁFICA DE QinFiltrado-QoutFiltrado
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','FLOWSDIFFERENCE','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1900]); 
ylim(axes1,[-0.0002 0.0016]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
plot(t, QinFiltrado-QoutFiltrado,'color', [0, 0.45, 0.74], 'linewidth', 0.5), hold on
plot(t,QLeakTotalObs,'--','color',[0.85, 0.33, 0.10],'linewidth',1), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^{3}/s)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$|Q_{in}-Q_{out}|$';
%print(gcf,fullfile(output_folder, 'QRESTA.eps'), '-depsc');




%% GRÁFICA DE LAM1,LAM2,LAM3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','magnitudleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
Ultvalor = sum(lamda1Fuga1Observador(10000:end))/length(lamda1Fuga1Observador(10000:end))
extendido = ones(1,120000)*Ultvalor;
Lam1Ext = [lamda1Fuga1Observador,extendido];
Ultvalor = sum(lamda2Fuga2Observador(72000:end))/length(lamda2Fuga2Observador(72000:end))
extendido = ones(1,59900)*Ultvalor;
Lam2Ext = [lamda2Fuga2Observador,extendido];
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1900]); 
ylim(axes1,[0 2e-04]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
% Línea azul punteada
plot(t, Lam1Ext, 'Color', [0, 0.45, 0.74], 'LineWidth', 0.5), hold on

% Línea roja punteada
plot(t, Lam2Ext, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 0.5), hold on

% Línea verde con marcador 'x'
plot(t, lambda3Fuga3Observador, 'Color', [0.47, 0.67, 0.19], ...
     'LineWidth', 0.5), hold on 
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^{5/2}/s)','FontName','Times New Roman','FontSize', 20);
%print(gcf,fullfile(output_folder, 'lambdasmagnitudes.eps'), '-depsc');







%% GRÁFICA DE z1,z2,z3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1900]); 
ylim(axes1,[10 52]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
z1Real = ones(1,190000)*17;
z2Real = ones(1,190000)*33.5;
z3Real = ones(1,190000)*49.9;
z1Obs = mean(z1Fuga1Observador(10200:end));
exten = ones(1,120000)*z1Obs;
z1Obs = [z1Fuga1Observador,exten]*(Lr/Leq);

z2Obs = z2prom;
exten = ones(1,59900)*z2Obs;
z2Obs = [z2Fuga2Observador,exten]*(Lr/Leq);

z3Obs = z3Fuga3Observador*(Lr/Leq);
plot(t, z1Real, '--','color',[0.2 0.2 0.2], 'linewidth', 2), hold on
plot(t, z1Obs,'color', [0, 0.45, 0.74], 'linewidth', 1), hold on
plot(t, z2Real, '--', 'color',[0.4 0.4 0.4],'linewidth', 2), hold on
plot(t, z2Obs,  'color', [0.85, 0.33, 0.10], 'linewidth', 1), hold on
plot(t, z3Real, '--', 'color',[0.6 0.6 0.6],'linewidth', 2), hold on
plot(t, z3Obs,  'color', [0.47, 0.67, 0.19], 'linewidth', 1), hold on

%z1Obs_filtrado = movmean(z1Obs, 2000, 'Endpoints', 'shrink');  % este es el más común
%z2Obs_filtrado = movmean(z2Obs, 2000, 'Endpoints', 'shrink');  % este es el más común
%z3Obs_filtrado = movmean(z3Obs, 2000, 'Endpoints', 'shrink');  % este es el más común

%plot(t, z1Obs_filtrado,'color',  [0, 0.45, 0.74], 'linewidth', 1.5), hold on
%plot(t, z2Obs_filtrado,'color',  [0.85, 0.33, 0.10], 'linewidth', 1.5), hold on
%plot(t, z3Obs_filtrado,'color',  [0.47, 0.67, 0.19], 'linewidth', 1.5), hold on

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$z_{1}(t)$';
leyenda2 =  '$\hat{z}_{1}(t)$';
leyenda3 =  '$z_{2}(t)$';
leyenda4 =  '$\hat{z}_{2}(t)$';
%print(gcf,fullfile(output_folder, 'leakspositions.eps'), '-depsc');
tic

%% GRÁFICA DE ERRORES
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Errores','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Configuración de los límites y ticks del eje X
xlim(axes1,[0 1900]); 
ylim(axes1,[0 10]);
%yticks(0:1:10);   % <-- añade esta línea
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

errorz1mean = abs(z1Obs(end)-17.04);
exten = ones(1,120000)*errorz1mean;
errorz1G = [errorz1,exten];

errorz2mean = abs(z2Obs(end)-33.5);
exten = ones(1,59900)*errorz2mean;
errorz2G = [errorz2,exten];


% Colores en tonos de gris
plot(t, errorz1G, 'color', [0, 0.45, 0.74], 'linewidth', 1), hold on  % Gris oscuro
plot(t, errorz2G, 'color', [0.85, 0.33, 0.10], 'linewidth', 1), hold on  % Gris intermedio
plot(t, errorz3, 'color', [0.47, 0.67, 0.19], 'linewidth', 1), hold on  % Gris claro

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$\sum \psi_{Leak_1}(t)$';
leyenda2 = '$\sum \psi_{Leak_2}(t)$';
%leyenda = legend(leyenda1,leyenda2,'Location','southeast');
%print(gcf,fullfile(output_folder, 'Errores.eps'), '-depsc');
tic


















function Combinatoria = matrizVerdad(NoParametrosAdaptar)
numerobits=NoParametrosAdaptar;
tamano=2^numerobits;
matriz=zeros(tamano,numerobits);
n=0;
acumulador=0;
inicializador=0;

for columna=numerobits:-1:1
inicializador=(2^(n))+1;
incrementos = 2^(n+1);
    if inicializador<2 
    inicializador = 2;
    end
    for filas=incrementos:incrementos:tamano
        matriz(inicializador:filas,columna)=1;
        inicializador=inicializador+incrementos;
    end
n=n+1;
end

    Combinatoria = matriz;  % Asignar la matriz al argumento de salida
end


%%% FUNCIONES PARA FUGA2
function [adis,bdis,cdis,ddis,matricesA,matricesB] = evaluarF2(D,Leq,b,f1,f2,f3,z1,lam1,x1max,x1min,x2max,x2min,x3max,x3min,x4max,x4min,x5max,x5min,x6max,x6min,matrizcombinatoria,NoParametros,Cd,Dd,dt)


matricesA = cell(2^NoParametros,1);
matricesB = cell(2^NoParametros,1);
adis = cell(2^NoParametros, 1);
bdis = cell(2^NoParametros, 1);
cdis = cell(2^NoParametros, 1);
ddis = cell(2^NoParametros, 1);
SYS = cell(2^NoParametros, 1);
SYSd = cell(2^NoParametros, 1);

D=D;
Ar=pi*(D/2)^2;
L=Leq;
b=b;
f1=f1;
f2=f2;
f3=f3;
miu1=f1/(2*D*Ar);
miu2=f2/(2*D*Ar);
miu3=f3/(2*D*Ar);
g=9.78;

for i = 1:2^NoParametros

    
    % Leer el valor de cada bit en la fila i (todas las columnas de la fila i)
    fila = matrizcombinatoria(i, :);
    
    % Asignar x1 según el primer valor de la fila
    if fila(1) == 0
        x1 = x1min;
    else
        x1 = x1max;
    end
    
    % Asignar x2 según el segundo valor de la fila
    if fila(2) == 0
        x2 = x2min;
    else
        x2 = x2max;
    end
    
    % Asignar x3 según el tercer valor de la fila
    if fila(3) == 0
        x3 = x3min;
    else
        x3 = x3max;
    end
    
    % Asignar x4 según el cuarto valor de la fila
    if fila(4) == 0
        x4 = x4min;
    else
        x4 = x4max;
    end

    % Asignar x5 según el cuarto valor de la fila
    if fila(5) == 0
        x5 = x5min;
    else
        x5 = x5max;
    end

    % Asignar x6 según el cuarto valor de la fila
    if fila(6) == 0
        x6 = x6min;
    else
        x6 = x6max;
    end
    
phi1 = -miu1*x1;
phi2 = (g*Ar*(x2-1))/z1;
phi3 = (-g*Ar*x2^2)/(z1*x6);
phi4 = (b^2/(g*Ar*z1))*(1-((lam1*sqrt(abs(x2)))/(2*x1)));
phi5 = (-b^2/(g*Ar*z1))*(1+((lam1*sqrt(abs(x2)))/(2*x3)));
phi6 = (g*Ar)/(x6-z1);
phi7 = -miu2*x3;
phi8 = (-g*Ar)/(x6-z1);
phi9 = (b^2)/(g*Ar*(x6-z1));
phi10 = (-b^2)/(g*Ar*(x6-z1));
phi11 = (-b^2*sqrt(abs(x4)))/(g*Ar*(x6-z1));
phi12 = (g*Ar)/(L-x6);
phi13 = -miu3*x5;

matricesA{i} = [phi1 ,phi2, 0, 0 ,0 ,phi3, 0;
        phi4 ,0, phi5, 0, 0, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        0, 0, phi9 ,0 ,phi10, 0, phi11;
        0, 0, 0, phi12, phi13, 0, 0;
    zeros(2,7)];


        matricesB{i}=[(g*Ar)/z1, 0;
    0, 0;
    0, 0;
    0, 0;
    0, -(g*Ar)/(L-x6);
    0, 0;
    0, 0];
    
    dt = dt;
    % Crear el sistema continuo en espacio de estados
    SYS{i} = ss(matricesA{i}, matricesB{i}, Cd, Dd);

    % Discretizar el sistema
    SYSd{i} = c2d(SYS{i,1}, dt);

    % Extraer las matrices discretas
    [adis{i}, bdis{i}, cdis{i}, ddis{i}] = ssdata(SYSd{i,1});
    tic
end
    
end

function [adis,bdis,cdis,ddis,matricesA,matricesB] = evaluarF3(D,Leq,b,f1,f2,f3,f4,z1,lam1,z2,lam2,x1max,x1min,x2max,x2min,x3max,x3min,x4max,x4min,x5max,x5min,x6max,x6min,x7max,x7min,x8max,x8min,matrizcombinatoria,NoParametros,Cd,Dd,dt)


matricesA = cell(2^NoParametros,1);
matricesB = cell(2^NoParametros,1);
adis = cell(2^NoParametros, 1);
bdis = cell(2^NoParametros, 1);
cdis = cell(2^NoParametros, 1);
ddis = cell(2^NoParametros, 1);
SYS = cell(2^NoParametros, 1);
SYSd = cell(2^NoParametros, 1);

D=D;
Ar=pi*(D/2)^2;
L=Leq;
b=b;
f1=f1;
f2=f2;
f3=f3;
f4=f4;
miu1=f1/(2*D*Ar);
miu2=f2/(2*D*Ar);
miu3=f3/(2*D*Ar);
miu4=f4/(2*D*Ar);
g=9.78;
lambda1=lam1;
lambda2=lam2;

for i = 1:2^NoParametros

    
    % Leer el valor de cada bit en la fila i (todas las columnas de la fila i)
    fila = matrizcombinatoria(i, :);
    
    % Asignar x1 según el primer valor de la fila
    if fila(1) == 0
        x1 = x1min;
    else
        x1 = x1max;
    end
    
    % Asignar x2 según el segundo valor de la fila
    if fila(2) == 0
        x2 = x2min;
    else
        x2 = x2max;
    end
    
    % Asignar x3 según el tercer valor de la fila
    if fila(3) == 0
        x3 = x3min;
    else
        x3 = x3max;
    end
    
    % Asignar x4 según el cuarto valor de la fila
    if fila(4) == 0
        x4 = x4min;
    else
        x4 = x4max;
    end

    % Asignar x5 según el cuarto valor de la fila
    if fila(5) == 0
        x5 = x5min;
    else
        x5 = x5max;
    end

    % Asignar x6 según el cuarto valor de la fila
    if fila(6) == 0
        x6 = x6min;
    else
        x6 = x6max;
    end

    if fila(7) == 0
        x7 = x7min;
    else
        x7 = x7max;
    end

    if fila(8) == 0
        x8 = x8min;
    else
        x8 = x8max;
    end
    
%%Construyo la primera linea
phi1=-miu1*x1;
phi2=((g*Ar)*(x2-1))/(z1);
phi3=-(g*Ar*x2^2)/(z1*x8);
tic

%%Construyo la segunda linea
phi4=(b^2/(g*Ar*(z1)))*(1 - (lambda1*sqrt(abs(x2)))/(2*x1));
phi5=(-b^2/(g*Ar*(z1)))*(1 + (lambda1*sqrt(abs(x2)))/(2*x3));
    tic
%%Construyo la tercera linea
phi6=(g*Ar)/(z2-z1);
phi7=-miu2*x3;
phi8=-(phi6);
tic
%%construyo la cuarta linea
phi9=(b^2/(g*Ar*(z2-z1)))*(1 - (lambda2*sqrt(abs(x4)))/(2*x3));
phi10=(-b^2/(g*Ar*(z2-z1)))*(1 + (lambda2*sqrt(abs(x4)))/(2*x5));
tic
%%Construyo la quinta linea
phi11=(g*Ar)/(x8-z2);
phi12=-miu3*x5;
phi13=-phi11;
tic

%%Construyo la sexta linea
phi14=b^2/(g*Ar*(x8-z2));
phi15=-phi14;
phi16=-phi14*sqrt(abs(x6));

%%Construyo la septima linea
phi17=(g*Ar)/(L-x8);
phi18=-(miu4*x7);

matricesA{i} = [phi1 phi2 0 0 0 0 0 phi3 0;
    phi4 0 phi5 0 0 0 0 0 0;
    0 phi6 phi7 phi8 0 0 0 0 0;
    0 0 phi9 0 phi10 0 0 0 0;
    0 0 0 phi11 phi12 phi13 0 0 0;
    0 0 0 0 phi14 0 phi15 0 phi16;
    0 0 0 0 0 phi17 phi18 0 0;
    zeros(2,9)];


        matricesB{i}=[(g*Ar)/z1, 0;
    0, 0;
    0, 0;
    0, 0;
    0, 0;
    0, 0;
    0, -(g*Ar)/(L-x8);
    0, 0;
    0, 0];
    
    dt = dt;
    % Crear el sistema continuo en espacio de estados
    SYS{i} = ss(matricesA{i}, matricesB{i}, Cd, Dd);

    % Discretizar el sistema
    SYSd{i} = c2d(SYS{i,1}, dt);

    % Extraer las matrices discretas
    [adis{i}, bdis{i}, cdis{i}, ddis{i}] = ssdata(SYSd{i,1});
    tic
end
    
end

function mu_vals = obtenermusF2(Q1SchMax, Q1SchMin, H2SchMax, H2SchMin, Q2SchMax, Q2SchMin,H3SchMax,H3SchMin,Q3SchMax,Q3SchMin, z2SchMax, z2SchMin, MatrizCombinatoria, NoParametros)

    % Inicializamos el vector mu_vals con ceros
    mu_vals = zeros(2^NoParametros, 1);  % Se crea un vector de ceros con el tamaño adecuado

    for i = 1:(2^NoParametros)
        % Obtener la fila actual de la Matriz Combinatoria
        fila = MatrizCombinatoria(i, :);

        % Asignar los valores dependiendo del bit
        if fila(1) == 1
            Q1Sch = Q1SchMax;
        else
            Q1Sch = Q1SchMin;
        end

        if fila(2) == 1
            H2Sch = H2SchMax;
        else
            H2Sch = H2SchMin;
        end

        if fila(3) == 1
            Q2Sch = Q2SchMax;
        else
            Q2Sch = Q2SchMin;
        end

        if fila(4) == 1
            H3Sch = H3SchMax;
        else
            H3Sch = H3SchMin;
        end

        if fila(5) == 1
            Q3Sch = Q3SchMax;
        else
            Q3Sch = Q3SchMin;
        end

        if fila(6) == 1
            z2Sch = z2SchMax;
        else
            z2Sch = z2SchMin;
        end

        % Calcular y almacenar el valor de mu_vals(i)
        mu_vals(i) = Q1Sch * H2Sch * Q2Sch * H3Sch * Q3Sch * z2Sch;
    end

end

function mu_vals = obtenermusF3(Q1SchMax, Q1SchMin, H2SchMax, H2SchMin, Q2SchMax, Q2SchMin,H3SchMax,H3SchMin,Q3SchMax,Q3SchMin,H4SchMax,H4SchMin,Q4SchMax,Q4SchMin, z3SchMax, z3SchMin, MatrizCombinatoria, NoParametros)

    % Inicializamos el vector mu_vals con ceros
    mu_vals = zeros(2^NoParametros, 1);  % Se crea un vector de ceros con el tamaño adecuado

    for i = 1:(2^NoParametros)
        % Obtener la fila actual de la Matriz Combinatoria
        fila = MatrizCombinatoria(i, :);

        % Asignar los valores dependiendo del bit
        if fila(1) == 1
            Q1Sch = Q1SchMax;
        else
            Q1Sch = Q1SchMin;
        end

        if fila(2) == 1
            H2Sch = H2SchMax;
        else
            H2Sch = H2SchMin;
        end

        if fila(3) == 1
            Q2Sch = Q2SchMax;
        else
            Q2Sch = Q2SchMin;
        end

        if fila(4) == 1
            H3Sch = H3SchMax;
        else
            H3Sch = H3SchMin;
        end

        if fila(5) == 1
            Q3Sch = Q3SchMax;
        else
            Q3Sch = Q3SchMin;
        end

        if fila(6) == 1
            H4Sch = H4SchMax;
        else
            H4Sch = H4SchMin;
        end
        
        if fila(7) == 1
            Q4Sch = Q4SchMax;
        else
            Q4Sch = Q4SchMin;
        end
        
        if fila(8) == 1
            z3Sch = z3SchMax;
        else
            z3Sch = z3SchMin;
        end

        % Calcular y almacenar el valor de mu_vals(i)
        mu_vals(i) = Q1Sch * H2Sch * Q2Sch * H3Sch * Q3Sch * H4Sch * Q4Sch *z3Sch;
    end

end


