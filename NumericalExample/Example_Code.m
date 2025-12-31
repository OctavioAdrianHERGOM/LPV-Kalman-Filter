%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            EXAMPLE               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ADD YALMIP DIRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath(genpath('---'));

%EXAMPLE
addpath(genpath('C:\Users\adria\Desktop\DOCUMENTOS TESIS\solvers\YALMIP-master'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             ADD SOLVER/S DIRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath(genpath('---'));

%EXAMPLE
addpath(genpath('C:\Program Files\Mosek'));
addpath(genpath('C:\Users\adria\Desktop\DOCUMENTOS TESIS\solvers\SDPT3-4.0'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Scheduling Bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2max = 0.3;
x2min = 0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Continous time - LPV Model
%              System Vertex - Continous time polytopic Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A1 = [-1 x2min; 0 -1];
A2 = [-1 x2max; 0 -1];

B = [0 1].';
C = [1 0];
D = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             System Vertex - Discrete time polytopic Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;
SYS1=ss(A1,B,C,D);
SYSd1=c2d(SYS1,dt);
[adis1,bdis1,cdis1,ddis1]=ssdata(SYSd1);

SYS2=ss(A2,B,C,D);
SYSd2=c2d(SYS2,dt);
[adis2,bdis2,cdis2,ddis2]=ssdata(SYSd2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            STOP HERE AND RUN SIMULATION TO VISUALICE THE BEHAVIOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Ensure Observability each vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

O1 = rank(obsv(adis1,cdis1));
O2 = rank(obsv(adis2,cdis2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Declare LMI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P       = sdpvar(2,2);
Y       = sdpvar(2,2,'symmetric');
W1  = sdpvar(1,2);
W2  = sdpvar(1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Declare Arbitrary Qk and Rk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qk=[1 0; 0 1];
Rk=1000; 
% TRY TO RANGE RK 1 - 10000 
% (Rk = 10,000 ALL NOISE ELIMINATED)
% (RK = 1 NO ATTENUATED NOISE)


R_obs=Rk;
H = Qk.^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Solve LMIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdpvar gamma_LQR_obs;

Trace_cond  = [gamma_LQR_obs*eye(2) eye(2); eye(2) Y];


H1  = [ -Y , Y*adis1 - W1'*C , Y*H' , W1';
        adis1'*Y - C'*W1 , -Y , zeros(2,2) , zeros(2,1);
        H*Y , zeros(2,2) , -eye(2) , zeros(2,1);
        W1 , zeros(1,2) , zeros(1,2) , -inv(R_obs)];

H2  = [ -Y , Y*adis2 - W2'*C , Y*H' , W2';
        adis2'*Y - C'*W2 , -Y , zeros(2,2) , zeros(2,1);
        H*Y , zeros(2,2) , -eye(2) , zeros(2,1);
        W2 , zeros(1,2) , zeros(1,2) , -inv(R_obs)];

ops= sdpsettings('solver','MOSEK','verbose',0);

[opt_info]  = optimize([[H1<=0]+[H2<=0]+[Y>=0]+[Trace_cond>=0]],[gamma_LQR_obs],ops)

Y_sol  = value(Y);%  Y solution -> WARNING: LQR solves for -Y!!
W1_sol  = value(W1); % W1 solution
W2_sol  = value(W2); % W1 solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Gains for each vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L1=(inv(Y_sol))*W1_sol.'; 
L2=(inv(Y_sol))*W2_sol.'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Poles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poles1_obs  = eig(adis1 - L1*C);
poles2_obs  = eig(adis2 - L2*C);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Load Data Needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These data correspond to the original simulation performed using the proposed parameters. 
% Both the input and output signals are required. 
% Additionally, the real values of x1 ​ and x2 ​ are saved in order to plot 
% them and compare them with their corresponding estimated states. 
% 
%
% If it is desired to modify the parameters, the gray area must be uncommented 
% and the data saved again.

Datos = importdata('C:\Users\adria\Desktop\COREA\Input_Output.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   OBSERVER Implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = length(Datos.u.signals.values);
n=1;
%Observer Initial Conditions
x1 = 0;
x2= 0.1;


while (n <= steps)

            %Input and Output
            u = Datos.u.signals.values(n);
            y = Datos.y.signals.values(n);

%% THE FIRST STEP USES THE INITIAL CONDITION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Estimated Scheduling Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             %Scheduling Variable Bounds
            x2maxObserver = x2max;
            x2minObserver = x2min;
                
             SchedulingMin = (x2maxObserver - x2)/(x2maxObserver-x2minObserver);
             SchedulingMax = (1-SchedulingMin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Convex Property
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                SumMus(n)=SchedulingMin+SchedulingMax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Interpolated Gain | Observer | System update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%  STATE VECTOR
            X = [x1;x2];
            
            %%  PREDICTED STATE VECTOR (Heuns)
            x1Predict = (-x1+x2^2)*dt+x1;
            x2Predict = (-x2+u)*dt+x2;
            XPredict = [x1Predict;x2Predict];


            %% LPV
            Act1 = [-1 x2; 0 -1]*X;
            B1 = [0;1]*u;

            %% Predicted LPV (Heuns)
            Act2 = [-1 x2Predict; 0 -1]*XPredict;
            B2 = [0;1]*u;

            %% State Representation
            System_dx = [x1;x2]+(dt/2)*((Act1+Act2)+(B1+B2));

            %% OUTPUT
            dx1 = System_dx(1);
            dx2 = System_dx(2);

             L_interp            = SchedulingMin*L1 + SchedulingMax*L2;
                   
            Kk=L_interp;

            LPV_KALMAN_OBSERVER = System_dx + Kk * (y - dx1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Updated States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            x1 = LPV_KALMAN_OBSERVER(1);
            x2 = LPV_KALMAN_OBSERVER(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   IGNORE THIS, THIS IS FOR PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            EstimatedState_x1(n)=x1;
            EstimatedState_x2(n)=x2;

            RealState_x2(n) = Datos.Realx2.signals.values(:,:,n);
            RealState_x1(n) = y;
            psi_1(n) = SchedulingMin;
            psi_2(n) = SchedulingMax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   NEXT STEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n=n+1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts = 0.01;               
N = n-1;              
t = (0:N-1)*Ts;         

%% StatesPlot
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 30]); 
ylim(axes1,[-0 0.35]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on

h1 = plot(t, RealState_x1,'color',[0 0.45 0.74],'linewidth',1); hold on
h2 = plot(t, RealState_x2,'color',[0.85 0.33 0.10],'linewidth',1);
h3 = plot(t, EstimatedState_x1,'--','color',[0.2 0.2 0.2],'linewidth',2);
h4 = plot(t, EstimatedState_x2,'--','color',[0.4 0.4 0.4],'linewidth',2);

legend([h1 h2 h3 h4], ...
       {'$x_1$ Real','$x_2$ Real','$\hat{x}_1$','$\hat{x}_2$'}, ...
       'Location','best','FontSize',18);

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);


%% SchedulingFunctionsPlot
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 30]); 
ylim(axes1,[-0.1 1.1]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on

h1 = plot(t, psi_1,'color',[0 0.45 0.74],'linewidth',1); hold on
h2 = plot(t, psi_2,'color',[0.85 0.33 0.10],'linewidth',1);

legend([h1 h2], ...
       {'${\psi}_1(\hat{\varphi}(k))$','${\psi}_2(\hat{\varphi}(k))$'}, ...
       'Location','best','FontSize',18);

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);



%% ConvexProperty
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
figure('Name','posicionesleaks','NumberTitle','off','color', [1 1 1])
axes1 = axes('FontSize',20);
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 30]); 
ylim(axes1,[-0.1 1.1]); 
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');
hold on

h1 = plot(t, SumMus,'color',[0 0.45 0.74],'linewidth',1); hold on

legend([h1], ...
       {'$\sum_{i=1}^2{{\psi}_i(\hat{\varphi}(k))}$'}, ...
       'Location','best','FontSize',18);

grid on
xlabel('(s)','FontName','Times New Roman','FontSize', 20);
