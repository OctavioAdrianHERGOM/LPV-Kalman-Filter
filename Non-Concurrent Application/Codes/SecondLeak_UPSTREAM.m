%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPLICACIÓN DE RESULTADOS CON DATOS REALES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LIMPIO TODO
clear all
clear classes
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
clc
%close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       YALMIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('C:\Users\adria\Desktop\DOCUMENTOS TESIS\solvers\YALMIP-master'));

%Example:
%addpath(genpath('C:\Users\adria\Desktop\DOCUMENTOS TESIS\solvers\YALMIP-master'));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       AÑADIR SOLVER/S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('C:\Program Files\Mosek'));

%Example:
%addpath(genpath('C:\Program Files\Mosek'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       IMPORTAR DATOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USE THE DATA BASE datos_reales_3.csv%%

Datos = importdata('\datos_reales_3.csv');

%Example
%Datos = importdata('C:\Users\adria\Desktop\COREA\Non-Concurrent Application\Datasets\datos_reales_3.csv');




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
%%                     SINGLE LEAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = QoutFiltrado((1:70000),1);
QinFiltrado = QinFiltrado((1:70000),1);
HoutFiltrado = HoutFiltrado((1:70000),1);
HinFiltrado = HinFiltrado((1:70000),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       VECTOR DE TIEMPO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;
fs = 1/dt;
muestras = length(QinFiltrado);
tiempo = [dt:dt:(muestras/fs)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       PARÁMETROS DE LA TUBERÍA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lr=68.147;                            %longitud física de la tubería (m)
g=9.81;                              %Fuerza de gravedad (m/s^2)
D=0.06271;%                          %Diámetro interno de la tubería (m)
A=pi*(D/2)^2;                        %Área interna del tubo (m^2)
e=0.013095;                           %Espesor de la pared del tubo (m)
eps=7e-6;                        %Rugosidad de la tubería (m)
%E = 8e+08 %% SACADA DE LA TESIS DE Erick Axel Padilla García PAG 36
%E=(34.247*((temp-55)/30.277)^5 + 94.874*((temp-55)/30.277)^4 - 323.75*((temp-55)/30.277)^3 + 799.98*((temp-55)/30.277)^2 - 2549.7*((temp-55)/30.277) + 3368.6)*g*1e4;    %módulo de elasticidad del tubo 1Pa=N/m^2=J/m^3=kg/m.s^2 [1kg/m^2=9.8Pa]
E=6.8646e+08; %% How can the temperature affect the performance of a classical pipeline model when plastic pipes are used?
%E=8.0e+08;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Saco Temperatura Promedio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TemperaturaFiltrada = Datos.data(:, 9);

%Temp_av = sum(TemperaturaFiltrada(10000,1))/length(TemperaturaFiltrada(10000,1));
Temp_av = sum(TemperaturaFiltrada((1:70000),1))/length(TemperaturaFiltrada((1:70000),1));
%Temp_av = sum(TemperaturaFiltrada)/length(TemperaturaFiltrada);

%Temp_av = TemperaturaFiltrada(1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Saco Qin0, Qout0, Hin0, Hout0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q0in_av= mean(QinFiltrado(1:8500,1));
Q0out_av= mean(QoutFiltrado(1:8500,1));
H0in_av= mean(HinFiltrado(1:8500,1));
H0out_av= mean(HoutFiltrado(1:8500,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Saco Rein,Reout para f1 y f2, miu1, miu2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rein=((Q0in_av/A)*D)/(Viscosidad_av);
Reout=((Q0in_av/A)*D)/(Viscosidad_av);


 %f1=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Rein)))^-1)^2;
f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
 %f2=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Reout)))^-1)^2;
f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;
miu1=f1/(2*D*A);
miu2=f2/(2*D*A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Leq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V1=2*(H0in_av-H0out_av)*D*g*A^2/(f1*Q0in_av^2);
V2=2*(H0in_av-H0out_av)*D*g*A^2/(f2*Q0out_av^2);
Leq=(V1+V2)/2;
z1 = 59;
%z1 = Leq*0.5;
H20=15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Declaro nuevamente E,ro,D,e,K para obtener b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=E;
ro = Densidad_av;
D = D;
e=e;
K=MEA_av;

b=sqrt((K/ro)/(1+K*D/(e*E)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Inicializo el observador
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X1=Q0in_av;
X2=15;
X3=Q0out_av;
X4=z1;
X5=0;

X1i=X1;
X2i=X2;
X3i=X3;
X4i=X4;%z1
X5i=X5;%0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Límite de los vértices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1max=max(QinFiltrado)+0.0001;
x1min=min(QinFiltrado)-0.0001;

x2max=17;
x2min=11;

x3max=max(QoutFiltrado)+0.0001;
x3min=min(QoutFiltrado)-0.0001;

x4max=85;
x4min=20;

Cd=[1 0 0 0 0;0 0 1 0 0];
Dd=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Ecuaciones LMI
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


%% ORIGINALES
Qk=[1*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-6 0 0;0 0 0 20000 0;0 0 0 0 1*10^-7];
Rk=[1*10^-5 0;0 1*10^-5]; %

%% CANDIDATA V1
Qk=[1*10^-7 0 0 0 0;0 0.001 0 0 0;0 0 1*10^-8 0 0;0 0 0 20000 0;0 0 0 0 1*10^-8];

%% CANDIDATA V2
Qk=[1*10^-7 0 0 0 0;0 0.001 0 0 0;0 0 1*10^-8 0 0;0 0 0 100 0;0 0 0 0 1*10^-9];

%% CANDIDATA V2
Qk=[1*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-8 0 0;0 0 0 10000 0;0 0 0 0 1*10^-9];


% % ORIGINALES
Qk=[5*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-6 0 0;0 0 0 10000 0;0 0 0 0 1*10^-9];
Rk=[1*10^-5 0;0 1*10^-5]; %

%% SIRVE PERO ES LENTA
%%Qk=[7.5*10^-5 0 0 0 0;0 0.01 0 0 0;0 0 1*10^-5 0 0;0 0 0 3000 0;0 0 0 0 1*10^-7];


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
[opt_info]  = optimize(Constraints,gamma_LQR_obs,ops);

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   L en los vértices
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
errorcuadraticoz1=0;
while (n <= muestras)
    % Filtrar entradas en cada iteración
    Qin = QinFiltrado(n, 1);
    Qout = QoutFiltrado(n, 1);
    Hin = HinFiltrado(n, 1);
    Hout = HoutFiltrado(n, 1);


    if (Qin - Qout) > 3e-04 && enclave ==0
        enclave = 1;
    end

    if enclave ==1
            Q1max = max(QinFiltrado)+0.001;
            Q1min = min(QinFiltrado);
            H2max = 17;
            H2min = 11;
            Q2max = max(QoutFiltrado);
            Q2min = min(QoutFiltrado)-0.0001;
            z1max = 85;
            z1min = 20;

            % Normalización (escalado) de las variables
            Q1Sch = (Q1max - x1) / (Q1max - Q1min);
            H2Sch = (H2max - x2) / (H2max - H2min);
            Q2Sch = (Q2max - x3) / (Q2max - Q2min);
            z1Sch = (z1max - x4) / (z1max - z1min);
            tic
            if Q1Sch <=0 || H2Sch <=0 || Q2Sch<=0 || z1Sch<=0
                tic
                tic
            end
            % Cálculo de estados del sistema
            tic
            Q1 = x1;
            H2 = x2;
            Q2 = x3;
            z1 = x4;
            lam = x5;
            Hk = [1 0 0 0 0; 0 0 1 0 0];

            % Parámetros y cálculos para dinámica del sistema
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
miu1=f1/(2*D*A);
miu2=f2/(2*D*A);
% Varvis1(n)=miu1;
% Varvis2(n)=miu2;
 %E=(34.247*((Temp_av-55)/30.277)^5 + 94.874*((Temp_av-55)/30.277)^4 - 323.75*((Temp_av-55)/30.277)^3 + 799.98*((Temp_av-55)/30.277)^2 - 2549.7*((Temp_av-55)/30.277) + 3368.6)*g*1e4;

            b=sqrt((K/ro)/(1+K*D/(e*E)));
            Q11=(-g*A*(H2-Hin)/(z1)-f1*Q1^2/(2*D*A))*dt+Q1;
            H21=(-b^2*(Q2-Q1+lam*sqrt(H2))/(g*A*z1))*dt+H2; 
            Q21=(-g*A*(Hout-H2)/(promediolongitud-z1)-f2*Q2^2/(2*D*A))*dt+Q2;
            tic
            %% Matrices para la dinámica del sistema
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

            %% Segunda parte de la dinámica del sistema
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

            %% Actualización de estados
            xx = [Q1; H2; Q2; z1; lam];
            dXmat = (Ax1 + Ax2) + (Bu1 + Bu2) + xx;
            dQ1 = dXmat(1);
            dH2 = dXmat(2);
            dQ2 = dXmat(3);
            dXm = [dQ1; dH2; dQ2; z1; lam];



            %% Cálculo de los coeficientes de pertenencia `mu`
            Q1SchMax= (1 - Q1Sch);
            Q1SchMin= Q1Sch;
            H2SchMax= (1 - H2Sch);
            H2SchMin= H2Sch;
            Q2SchMax= (1 - Q2Sch);
            Q2SchMin= Q2Sch;
            z1SchMax= (1-z1Sch);
            z1SchMin = z1Sch;

            %% FUNCIÓN PARA OBTENER LA COMBINATORIA DE MU'S
mu_1(n)= Q1Sch * H2Sch * Q2Sch * z1Sch; %todos MIN
mu_2(n)= (1-Q1Sch) * H2Sch * Q2Sch * z1Sch;
mu_3(n)= Q1Sch * (1-H2Sch) * Q2Sch * z1Sch;
mu_4(n)= Q1Sch * H2Sch * (1-Q2Sch) * z1Sch;

mu_5(n)= Q1Sch * H2Sch * Q2Sch * (1-z1Sch);
mu_6(n)= (1-Q1Sch) * (1-H2Sch) * (1-Q2Sch) * (1-z1Sch); %todos MAX
mu_7(n)= Q1Sch * (1-H2Sch) * (1-Q2Sch) * (1-z1Sch);
mu_8(n)= (1-Q1Sch) * H2Sch * (1-Q2Sch) * (1-z1Sch);

mu_9(n)= (1-Q1Sch) * (1-H2Sch) * Q2Sch * (1-z1Sch);
mu_10(n)= (1-Q1Sch) * (1-H2Sch) * (1-Q2Sch) * z1Sch;
mu_11(n)= (1-Q1Sch) * (1-H2Sch) * Q2Sch * z1Sch;
mu_12(n)= (1-Q1Sch) * H2Sch * (1-Q2Sch) * z1Sch;
 
mu_13(n)= (1-Q1Sch) * H2Sch * Q2Sch * (1-z1Sch);
mu_14(n)= Q1Sch * H2Sch * (1-Q2Sch) * (1-z1Sch);
mu_15(n)= Q1Sch * (1-H2Sch) * Q2Sch * (1-z1Sch);
mu_16(n)= Q1Sch * (1-H2Sch) * (1-Q2Sch) * z1Sch;           

            %% Interpolación y cálculo de L_interp (corregido)
%            L_interp = 0;
musL1(n)=mu_1(n)+mu_2(n)+mu_3(n)+mu_4(n)+mu_5(n)+mu_6(n)+mu_7(n)+mu_8(n)+mu_9(n)+mu_10(n)+mu_11(n)+mu_12(n)+mu_13(n)+mu_14(n)+mu_15(n)+mu_16(n);

 L_interp            = mu_1(n)*L1 + mu_2(n)*L2 + mu_3(n)*L3 + mu_4(n)*L4...
                       + mu_5(n)*L5 + mu_6(n)*L6 + mu_7(n)*L7 + mu_8(n)*L8...
                       + mu_9(n)*L9 + mu_10(n)*L10 + mu_11(n)*L11 + mu_12(n)*L12...
                       + mu_13(n)*L13 + mu_14(n)*L14 + mu_15(n)*L15 + mu_16(n)*L16;
                   
          Kk=L_interp;
            dX = dXm + Kk * ([Qin; Qout] - [dQ1; dQ2]);

            %% Actualización de las variables dinámicas
            x1 = dX(1);
            x2 = dX(2);
            x3 = dX(3);
            x4 = dX(4);
            x5 = dX(5);

    else
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
          errorz1(n) = abs(49.9-(z1Fuga1Observador(n)*(Lr/Leq)));
            if enclave == 1 
            errorcuadraticoz1 = errorcuadraticoz1 + abs(49.9-(z1Fuga1Observador(n)*(Lr/Leq)));
            end
        n=n+1;
end
RMSEz1 = sqrt(errorcuadraticoz1 / (70001-10300));

output_folder = 'E:\Tesis_Codigos_DatosReales\FUNCIONANDO_BIEN\ImagenesTesis\Dos Fugas Izquierdo';
tic
Ts = 0.01; % tiempo de muestreo
t = 0:Ts:(length(QinFiltrado)-1)*Ts;
mu = {mu_1, mu_2, mu_3, mu_4, mu_5, mu_6, mu_7, mu_8, ...
      mu_9, mu_10, mu_11, mu_12, mu_13, mu_14, mu_15, mu_16};


%% GRAFICO LA EVOLUCION DE MU'S
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','MUS EVOLUTION    ','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[75 700]); 
ylim(axes1,[-0.01 0.45]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Genera colores RGB
colors = lines(64);  % Genera 64 colores distintos
% Plotea cada mu_Leak2 en un color diferente usando un ciclo for
for l = 1:16
    % Plotea cada mu_Leak2 con su color correspondiente
    plot(t, mu{l}, 'Color', colors(l, :)); hold on;
    tic
end
grid on;
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
%print(gcf, fullfile(output_folder, 'musEvolutionL1.eps'), '-depsc');
tic

tic
tamano = length(lamda1Fuga1Observador);
for i =1 :1 : tamano
    Ql(i)=lamda1Fuga1Observador(i)*sqrt(H2Fuga1Observador(i));
end

z1prom = sum(z1Fuga1Observador(1,(10500:end)))/length(z1Fuga1Observador(1,(10500:end)));
tic
%z1prom = z1Real(1,53730:end)
%z1prom = sum(z1prom)/length(z1prom);
% lambda1prom = lamda1Fuga1Observador(1,53730:end);
% lambda1prom  = sum(lambda1prom)/length(lambda1prom)
lambda1prom = sum(lamda1Fuga1Observador(1,(11000:end)))/length(lamda1Fuga1Observador(1,(11000:end)));
z1prom = z1prom;
H2prom = sum(H2Fuga1Observador(1,(11000:end)))/length(H2Fuga1Observador(1,(11000:end)));


%%-------------------------------------------------------------------------

%% Z1 ESTÁ EN LONGITUD EQUIVALENTE


z1=z1prom;
lam1 = lambda1prom;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DATOS FILTRADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QoutFiltrado = Datos.data(:,1)*1e-03;
QinFiltrado = Datos.data(:,2)*1e-03;
HinFiltrado = Datos.data(:,5);
HoutFiltrado = Datos.data(:,6);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       DATOS FILTRADOS 
%%                        PARA UNA FUGA
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


TemperaturaFiltrada = Datos.data(:, 9);
tic
%Temp_av = sum(TemperaturaFiltrada(10000,1))/length(TemperaturaFiltrada(10000,1));
Temp_av = sum(TemperaturaFiltrada((1:end),1))/length(TemperaturaFiltrada((1:end),1));
Temp_av = TemperaturaFiltrada(1);
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



z0=Leq*0.57;
x1 = Q1Fuga1Observador(end);
x2 = 15;
x3 = Q2Fuga1Observador(end);
x4 = 17;
x5 = Q2Fuga1Observador(end);
x6 = z0;
x7 = 0;
tic
x1i = x1;
x2i = x2;
x3i = x3;
x4i = x4;
x5i = x5;
x6i = x6;
x7i = x7;

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

x2max=16.5;
x2min=10;

tic
%% ESTAS TAMBIÉN FUNCIONAN, PERO EN TEORÍA DESCONOZCO EL FLUJO EN Q2 (SACADO DE SIMULINK)
%%x3max = 0.020
%%x3min = 0.015

x3max=(x1max-(lam1*sqrt(x2max)))+0.0004;
x3min= (x1min-(lam1*sqrt(x2min)))-0.0001;

x4max=19;
x4min=13.5;

x5max = max(QoutFiltrado)+0.0001;
x5min = min(QoutFiltrado)-0.0001;

x6max = 60;
x6min = 30;

%%90,39 es el mínimo. Pero las estimaciones salen con mucho ruido
%%80,50 es lo mejor

% x3max = (x1max);
% x3min = (x5min)+0.0005;

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
Reinter = ((Qintermediate/A)*D)/(Viscosidad_av);
Reout=((Q0in_av/A)*D)/(Viscosidad_av);

% f1=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Rein)))^-1)^2;
f1=0.25 / ( log10( eps/(3.7*D) + 5.74/(Rein) ) ).^2;
% f2=((-1.8*log10(((eps/D)/3.7)^1.11+(6.9/Reout)))^-1)^2;
f2=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reout) ) ).^2;
f3=0.25 / ( log10( eps/(3.7*D) + 5.74/(Reinter) ) ).^2;
miu1=f1/(2*D*A);
miu2=f2/(2*D*A);
miu3=f3/(2*D*A);

E=E;
ro = Densidad_av;
D = D;
e=e;
K=MEA_av;

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


%% BUENA
Qk=[1*10^-5 0 0 0 0 0 0;
    0 1*10^-8 0 0 0 0 0;
    0 0 1*10^-5 0 0 0 0;
    0 0 0 1*10^-8 0 0 0; 
    0 0 0 0 1*10^-5 0 0;
    0 0 0 0 0 8000 0;
    0 0 0 0 0 0 1*10^-8];

Qk=[1*10^-5 0 0 0 0 0 0;
    0 1*10^-2 0 0 0 0 0;
    0 0 1*10^-5 0 0 0 0;
    0 0 0 1*10^-2 0 0 0;
    0 0 0 0 1*10^-5 0 0;
    0 0 0 0 0 1000 0;
    0 0 0 0 0 0 1*10^-6];

% Qk=[1*10^-5 0 0 0 0 0 0;
%     0 1*10^-2 0 0 0 0 0;
%     0 0 1*10^-5 0 0 0 0;
%     0 0 0 1*10^-2 0 0 0;
%     0 0 0 0 1*10^-5 0 0;
%     0 0 0 0 0 10000 0;
%     0 0 0 0 0 0 1*10^-7];

% %% QK ORIGINAL - NO ES MUY BUENA
% Qk=[1*10^-5 0 0 0 0 0 0;
%     0 1*10^-2 0 0 0 0 0;
%     0 0 1*10^-5 0 0 0 0;
%     0 0 0 1*10^-2 0 0 0;
%     0 0 0 0 1*10^-5 0 0;
%     0 0 0 0 0 2000 0;
%     0 0 0 0 0 0 1*10^-7];


% Qk=[1*10^-4 0 0 0 0 0 0;
%     0 0.001 0 0 0 0 0;
%     0 0 1*10^-6 0 0 0 0;
%     0 0 0 1*10^-3 0 0 0; 
%     0 0 0 0 1*10^-6 0 0;
%     0 0 0 0 0 2000 0;
%     0 0 0 0 0 0 1*10^-8];

% Qk=[1*10^-5 0 0 0 0 0 0;
%     0 0.01 0 0 0 0 0;
%     0 0 1*10^-5 0 0 0 0;
%     0 0 0 1*10^-2 0 0 0; 
%     0 0 0 0 1*10^-4 0 0;
%     0 0 0 0 0 8000 0;
%     0 0 0 0 0 0 1*10^-5];

%% AHORA INTENTO ELIMINAR EL PICO DE H3
% %% BUENA TAMBIÉN
% Qk=[1*10^-4 0 0 0 0 0 0;
%     0 1 0 0 0 0 0;
%     0 0 8*10^-6 0 0 0 0;
%     0 0 0 1 0 0 0;
%     0 0 0 0 1*10^-3 0 0;
%     0 0 0 0 0 4000 0;
%     0 0 0 0 0 0 1*10^-5];
% 
% %% VAMOS AVANZANDO, ESTA BUENA ESTA
% Qk=[1*10^-4 0 0 0 0 0 0;
%     0 1 0 0 0 0 0;
%     0 0 8*10^-6 0 0 0 0;
%     0 0 0 1 0 0 0;
%     0 0 0 0 2.5*10^-3 0 0;
%     0 0 0 0 0 6000 0;
%     0 0 0 0 0 0 1*10^-5];
% 
% Qk=[1*10^-4 0 0 0 0 0 0;
%     0 0.01 0 0 0 0 0;
%     0 0 8*10^-6 0 0 0 0;
%     0 0 0 0.01 0 0 0;
%     0 0 0 0 2*10^-3 0 0;
%     0 0 0 0 0 6000 0;
%     0 0 0 0 0 0 5*10^-9];


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
mu=0;
errorcuadraticoz2=0;
while(n<=iteraciones)
tic
    Qin = QinFiltrado(n,1);
    Qout = QoutFiltrado(n,1);
    Hin = HinFiltrado(n,1);
    Hout = HoutFiltrado(n,1);

    if (Qin-Qout>=6e-04)  && enclave ==0 && n>=68000
        enclave =1;
        tic
    end

    if enclave ==1
        Q1max=max(QinFiltrado)+0.0001;
        Q1min=min(QinFiltrado)-0.0001;

        H2max=16.5;
        H2min=10;

        % 
        Q2max=(Q1max-(lam1*sqrt(H2max)))+0.0004;
        Q2min= Q1min-(lam1*sqrt(H2min))-0.0001;

      %  Q2max = 0.025;
      %  Q2min = 0.005;

        H3max=19;
        H3min=13.5;

        Q3max = max(QoutFiltrado)+0.0001;
        Q3min = min(QoutFiltrado)-0.0001;


        z2max = 60;
        z2min = 30;

        Q1Sch = (Q1max - x1) / (Q1max - Q1min);
        H2Sch = (H2max - x2) / (H2max - H2min);
        Q2Sch = (Q2max - x3) / (Q2max - Q2min);
        H3Sch = (H3max - x4) / (H3max - H3min);
        Q3Sch = (Q3max - x5) / (Q3max - Q3min);
        z2Sch = (z2max - x6) / (z2max - z2min);

        if Q1Sch < 0 || H2Sch <0 || Q2Sch <0 || H3Sch<0 || Q3Sch <0 || z2Sch <0
        tic
        end
        if n >=76400
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
        L=Leq;
        dz1=(z2);
        dz2=(z1-z2);
        dz3=(L-z1);

        %% HAGO MIS Q ESTIMADAS (Qi gorro j+1)
        Q11=Q1+dt*((-((g*A)*(H3-Hin))/(dz1))-miu1*Q1^2);
        Q21=Q2+dt*((-((g*A)*(H2-H3))/(dz2))-miu2*Q2^2);
        Q31=Q3+dt*((-((g*A)*(Hout-H2))/(dz3))-miu3*Q3^2);
        
        %% HAGO MIS H ESTIMADAS (Hi gorro j+1)
        H21=H2+dt*(-(b^2)*(Q3-Q2+lam1*sqrt(abs(H2)))/(g*A*dz2));
        H31=H3+dt*(-(b^2)*(Q2-Q1+lam*sqrt(abs(H3)))/(g*A*dz1));

        phi1 = -miu1*Q11;
        phi2 = (g*Ar*(H31-1))/z2;
        phi3 = (-g*Ar*H31^2)/(z2^2);

        phi4 = (b^2/(g*Ar*(z1-z2)))*(1-((lam1*sqrt(abs(H21)))/(2*Q21)));
        phi5 = (-b^2/(g*Ar*(z1-z2)))*(1+((lam1*sqrt(abs(H21)))/(2*Q31)));

        phi6 = -(g*Ar)/(z1-z2);
        phi7 = -miu2*Q21;
        phi8 = (g*Ar)/(z1-z2);

        phi9 = (b^2)/(g*Ar*(z2));
        phi10 = (-b^2)/(g*Ar*(z2));
        phi11 = (-b^2*sqrt(abs(H31)))/(g*Ar*(z2));
        
        phi12 = (g*Ar)/(L-z1);
        phi13 = -miu3*Q31;
        
        Aact1 = [phi1 ,0, 0, phi2 ,0 ,phi3, 0;
        0 ,0, phi4, 0, phi5, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        phi9, 0, phi10 ,0 ,0, 0, phi11;
        0, phi12, 0, 0, phi13, 0, 0;
        zeros(2,7)];
        
        xact1 = [Q11, H21, Q21, H31, Q31, z2, lam]' * (dt / 2);
        Bact1=[(g*Ar)/x6, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z1);
        0, 0;
        0, 0];
        u = [Hin * (dt / 2); Hout * (dt / 2)];
        Ax1 = Aact1 * xact1;
        Bu1 = Bact1 * u;

        phi1 = -miu1*Q1;
        phi2 = (g*Ar*(H3-1))/z2;
        phi3 = (-g*Ar*H3^2)/(z2^2);

        phi4 = (b^2/(g*Ar*(z1-z2)))*(1-((lam1*sqrt(abs(H2)))/(2*Q2)));
        phi5 = (-b^2/(g*Ar*(z1-z2)))*(1+((lam1*sqrt(abs(H2)))/(2*Q3)));

        phi6 = -(g*Ar)/(z1-z2);
        phi7 = -miu2*Q2;
        phi8 = (g*Ar)/(z1-z2);

        phi9 = (b^2)/(g*Ar*(z2));
        phi10 = (-b^2)/(g*Ar*(z2));
        phi11 = (-b^2*sqrt(abs(H3)))/(g*Ar*(z2));

        phi12 = (g*Ar)/(L-z1);
        phi13 = -miu3*Q3;
        
        Aact2 = [phi1 ,0, 0, phi2 ,0 ,phi3, 0;
        0 ,0, phi4, 0, phi5, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        phi9, 0, phi10 ,0 ,0, 0, phi11;
        0, phi12, 0, 0, phi13, 0, 0;
        zeros(2,7)];
        
        xact2 = [Q1, H2, Q2, H3, Q3, z2, lam]' * (dt / 2);
        Bact2=[(g*Ar)/x6, 0;
        0, 0;
        0, 0;
        0, 0;
        0, -(g*Ar)/(L-z1);
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
            errorcuadraticoz2 = errorcuadraticoz2 + abs(33.47-(z2Fuga2Observador(n)*(Lr/Leq)));
            end
    n=n+1;
end
RMSEz2 = sqrt(errorcuadraticoz2 / (130000-70200));
tic

tic
Ts = 0.01; % tiempo de muestreo
t = 0:Ts:(length(QinFiltrado)-1)*Ts;

Qleak1 = lamda1Fuga1Observador.*sqrt(H2Fuga1Observador);
Qleak1Total = Qleak1;

Qleak2 = lamda2Fuga2Observador.*sqrt(H3Fuga2Observador);
Qleak2Total = Qleak2(70001:end)+lam1*sqrt(H2prom);

QLeakTotalObs = [Qleak1Total, Qleak2Total];

tic
%% GRAFICO LA EVOLUCION DE MU'S
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','MUS EVOLUTION    ','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[600 1300]); 
ylim(axes1,[-0.01 0.2]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
% Genera colores RGB
colors = lines(64);  % Genera 64 colores distintos
% Plotea cada mu_Leak2 en un color diferente usando un ciclo for
for l = 1:64
    % Plotea cada mu_Leak2 con su color correspondiente
    plot(t, mu_Leak2{l}, 'Color', colors(l, :)); hold on;
    tic
end
grid on;
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
%print(gcf, fullfile(output_folder, 'musEvolution.eps'), '-depsc');
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

musL1G = [musL1,ones(1,60001)];
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
%print(gcf,fullfile(output_folder, 'MUS.eps'), '-depsc');
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
ylim(axes1,[0.0085 0.0096]); 
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
plot(t,QLeakTotalObs,'--','color',[0.85, 0.33, 0.10],'linewidth',1), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m^{3}/s)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$|Q_{in}-Q_{out}|$';
%print(gcf,fullfile(output_folder, 'QRESTA.eps'), '-depsc');


tic

%% GRÁFICA DE LAM1,LAM2,LAM3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','magnitudleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
Ultvalor = mean(lamda1Fuga1Observador(11000:end));
Ultvalor = Ultvalor(end)
extendido = ones(1,60001)*Ultvalor;
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





% 
% 
% Delta = 1000;
% N = 2 * Delta + 1;
% 
% 
% %% SELECCIONO DESDE QUE YA HIZO CONVERGENCIA
% data_filterrorz2 = data_filtz2(75000:end);
% zl2 = 33.5 * ones(size(data_filterrorz2));  
% data_filtPromz2 = sum(data_filtz2(1,(75000:end)))/length(data_filtz2(1,(75000:end)));
% error_normz2 = sqrt( sum( (zl2 - data_filterrorz2).^2 )/55000 );

%% GRÁFICA DE z1,z2,z3
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 1300]); 
ylim(axes1,[20 60]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on
z1Real = ones(1,130001)*49.9;
z2Real = ones(1,130001)*33.5;
z1Obs = z1prom;
exten = ones(1,60001)*z1Obs;
%exten2 = ones(1,60001)*data_filtPromz1;
z1Obs = [z1Fuga1Observador,exten]*(Lr/Leq);
%data_filtObs = [data_filtz1,exten2];
z2Obs = [z2Fuga2Observador]*(Lr/Leq);
plot(t, z1Real, '--','color',[0.2 0.2 0.2], 'linewidth', 2), hold on
plot(t, z1Obs,'color', [0, 0.45, 0.74], 'linewidth', 1), hold on
%plot(t, data_filtObs, 'color', [1, 0.0, 0.0], 'linewidth', 1)
plot(t, z2Real, '--', 'color',[0.4 0.4 0.4],'linewidth', 2), hold on
plot(t, z2Obs,  'color', [0.85, 0.33, 0.10], 'linewidth', 1), hold on
%plot(t, data_filtz2, 'color', [0.0, 0.5, 0.00], 'linewidth', 1), hold on

% plot(t,z2Real,'k--','linewidth',0.5), hold on
grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 =  '$z_{1}(t)$';
leyenda2 =  '$\hat{z}_{1}(t)$';
leyenda3 =  '$z_{2}(t)$';
leyenda4 =  '$\hat{z}_{2}(t)$';
%print(gcf,fullfile(output_folder, 'leakspositions.eps'), '-depsc');


%% GRÁFICA DE ERRORES
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','Errores','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Configuración de los límites y ticks del eje X
xlim(axes1,[0 1300]); 
ylim(axes1,[0 25]);
%yticks(0:1:10);   % <-- añade esta línea
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

z1Obs = mean(z1Fuga1Observador(10500:end)) * (Lr/Leq);

errorz1mean = abs(z1Obs(end)-49.9);
exten = ones(1,60001)*errorz1mean;
errorz1G = [errorz1,exten];

% Colores en tonos de gris
plot(t, errorz1G, 'color', [0, 0.45, 0.74], 'linewidth', 1), hold on  % Gris oscuro
plot(t, errorz2, 'color', [0.85, 0.33, 0.10], 'linewidth', 1), hold on  % Gris intermedio

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
ylabel('(m)','FontName','Times New Roman','FontSize', 20);
leyenda1 = '$\sum \psi_{Leak_1}(t)$';
leyenda2 = '$\sum \psi_{Leak_2}(t)$';
%leyenda = legend(leyenda1,leyenda2,'Location','southeast');
%print(gcf,fullfile(output_folder, 'Errores.eps'), '-depsc');
tic





% figure(1)
% plot(Q1Fuga2Observador)
% hold on
% plot(Q1Fuga1Observador)
% 
% figure(2)
% plot(Q3Fuga2Observador)
% hold on
% plot(Q2Fuga1Observador)
% 
% figure(3)
% plot(H2Fuga2Observador)
% hold on
% plot(H2Fuga1Observador)
% % 
% figure(4)
% plot(H3Fuga2Observador)
% 
% figure(5)
% colors = lines(64);
% for l = 1:64
%     Plotea cada mu_Leak2 con su color correspondiente
%     plot( mu_Leak2{l}, 'Color', colors(l, :)); hold on;
%     tic
% end
%% --------------------------------------------------------------------------------
% sum(z2Fuga2Observador(73300:122000))/length(z2Fuga2Observador(73300:122000))*(Lr/Leq)
% 
% Valoresz2 = z2Fuga2Observador(90000:end);
% n = length(z2Fuga2Observador(90000:end));
% counter =0;
% for i = 1:n
%     suma = abs ( 33.5- (Valoresz2(i)*(Lr/Leq)) ) / 33.5 ;
%     counter = counter + suma;
%     tic
% end
% 
% Errorporcentaje = (counter*100)/n;
% 
% z2prom = sum(z2Fuga2Observador(1,(72000:end)))/length(z2Fuga2Observador(1,(72000:end)));
% lambda2prom = sum(lamda2Fuga2Observador(1,(72000:end)))/length(lamda2Fuga2Observador(1,(72000:end)));
% 
% z2 = z2prom;
% lam2= lambda2prom;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       DATOS FILTRADOS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QoutFiltrado = Datos.data(:,1)*1e-03;
% QinFiltrado = Datos.data(:,2)*1e-03;
% HinFiltrado = Datos.data(:,5);
% HoutFiltrado = Datos.data(:,6);
% tic
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       DATOS FILTRADOS 
% %%                        PARA UNA FUGA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QoutFiltrado = QoutFiltrado((1:end),1);
% QinFiltrado = QinFiltrado((1:end),1);
% HoutFiltrado = HoutFiltrado((1:end),1);
% HinFiltrado = HinFiltrado((1:end),1);
% 
% iteraciones = length(QinFiltrado);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       Saco Temperatura Promedio
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TemperaturaFiltrada = Datos.data(:, 9);
% Temp_av = sum(TemperaturaFiltrada((1:end),1))/length(TemperaturaFiltrada(1:end));
% % Parametro = ((Temp_av-55)/30.277);
% % 
% % YoungsModulus = (34.247*Parametro.^5 + 94.875*Parametro.^4 - 323.75*Parametro.^3 + 799.98*Parametro.^2 - 2549.7*Parametro + 3368.6)*g*1e4;
% tic
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                      Hago Tablas Para Interpolar
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% DATOS DEL LIBRO PAG 288 1 BAR (COMPRENSIBILIDAD)
% Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
% Comprensibilidad = [0.50881e-9 0.49167e-9 0.4779e-9 0.46691e-9 0.45825e-9 0.45157e-9 0.4466e-9 0.44314e-9 0.441e-9 0.44005e-9 0.4402e-9 0.44341e-9 0.45014e-9 0.4601e-9 0.47311e-9];
% %% m^2/N
% %% DATOS DEL LIBRO PAG 296 (\emph{v} Viscocidad cinemática)
% Temperatura = [0 5 10 15 20 25 30 35 40 45 50 60 70 80 90];
% v = [1.7920 1.5182 1.3063 1.1386 1.0034 0.89266 0.80070 0.72344 0.65785 0.60166 0.55313 0.474 0.41273 0.36433 0.32547];
% v = v*1e-06; %% 
% %% m^2/s
% 
% %% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (vapor–liquid phase)
% Temperatura1 = [273.160 278 284 288 294 298 304 308 314 318 324 336 346 356 366];
% Temperatura1 = Temperatura1-(ones*273.15);
% Densidad1 = [999.793 999.919 999.575 999.079 997.983 997.042 995.346 994.042  991.848  990.235 987.610  981.671  976.086 969.972 963.363];
% %% kg/m^3
% %% DATOS PAPER https://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf (one-single phase)
% Temperatura2 = [273.15 275 280 285 295 300 305 310 315 325 335 340 345 355 365];
% Temperatura2 = Temperatura2-(ones*273.15);
% Densidad2 = [999.843 999.937 999.910  999.517 997.807  996.556 995.075 993.383  991.496  987.187  982.233  979.535  976.699  970.628 964.057];
% %% kg/m^3
% 
% for i=1:1:length(Comprensibilidad)
% MEA(i)=1/Comprensibilidad(i);
% Visabs(i)=v(i)*Densidad1(i);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                      Interpolo los valores
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Viscosidad_av = interp1(Temperatura,  v,Temp_av, 'spline','extrap'); % Interpolación si es necesario
% Comprensibilidad_av = interp1(Temperatura, Comprensibilidad, Temp_av, 'spline','extrap'); % Interpolación si es necesario
% Densidad_av = interp1(Temperatura, Densidad1, Temp_av, 'spline','extrap'); % Interpolación si es necesario
% MEA_av = interp1(Temperatura,  MEA,Temp_av, 'spline','extrap'); % Interpolación si es necesario
% 
% 
% 
% z2prom = sum(z2Fuga2Observador(1,(74500:end)))/length(z2Fuga2Observador(1,(74500:end)));
% 
% 





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
g=9.81;

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
phi2 = (g*Ar*(x4-1))/x6;
phi3 = (-g*Ar*x4^2)/(x6^2);

phi4 = (b^2/(g*Ar*(z1-x6)))*(1-((lam1*sqrt(abs(x2)))/(2*x3)));
phi5 = (-b^2/(g*Ar*(z1-x6)))*(1+((lam1*sqrt(abs(x2)))/(2*x5)));

phi6 = -(g*Ar)/(z1-x6);
phi7 = -miu2*x3;
phi8 = (g*Ar)/(z1-x6);

phi9 = (b^2)/(g*Ar*(x6));
phi10 = (-b^2)/(g*Ar*(x6));
phi11 = (-b^2*sqrt(abs(x4)))/(g*Ar*(x6));

phi12 = (g*Ar)/(L-z1);
phi13 = -miu3*x5;

matricesA{i} = [phi1 ,0, 0, phi2 ,0 ,phi3, 0;
        0 ,0, phi4, 0, phi5, 0, 0;
        0 ,phi6 ,phi7 ,phi8, 0 ,0, 0;
        phi9, 0, phi10 ,0 ,0, 0, phi11;
        0, phi12, 0, 0, phi13, 0, 0;
    zeros(2,7)];


        matricesB{i}=[(g*Ar)/x6, 0;
    0, 0;
    0, 0;
    0, 0;
    0, -(g*Ar)/(L-z1);
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

