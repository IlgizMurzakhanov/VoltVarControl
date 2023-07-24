% Code for Linear Constraint Quadratic Problem (LCQP)
clear; clear all; clc; close all;

tic

% Define case
structA = load('CasePv30_v7.mat'); 
caseA0 = structA.caseD;
PV_ind = [2:31];
WT_ind = [];
VmPFExcSla_List = [];

FlagSubsetActInv = 0;

%% 30 PVs, but only 5 have non-zero Q limits and placed at different depths
% structA = load('PV30_5Inv_Depth10_OtherNAS.mat'); 
% caseA0 = structA.caseB_Depth10;
% PV_ind = [2:31];
% WT_ind = [];

% structA = load('PV30_5Inv_Depth15_OtherNAS.mat'); 
% caseA0 = structA.caseB_Depth15;
% PV_ind = [2:31];
% WT_ind = [];

% structA = load('PV30_5Inv_Depth20_OtherNAS.mat'); 
% caseA0 = structA.caseB_Depth20;
% PV_ind = [2:31];
% WT_ind = [];

% This flag should be 1 for n out of m inverters to be active, otherwise 0
% FlagSubsetActInv = 1;

%% Set Pgen=Pload (to decrease total generation, so undervoltage problem becomes more evident)
caseA1 = caseA0;

NBusIncSla = size(caseA1.bus,1);
NGenIncSla = size(caseA1.gen,1);

NBusExcSla = NBusIncSla-1;
NGenExcSla = NBusIncSla-1;

for i = 2:NBusIncSla
    for j = 2:NGenIncSla
        BusId_BusStruc = caseA1.bus(i,1);
        BusId_GenStruc = caseA1.gen(j,1);
        if BusId_BusStruc == BusId_GenStruc
            caseA1.gen(j,[2,4,9,10]) = caseA1.bus(i,3);
            caseA1.gen(j,5) = -caseA1.bus(i,3);
        end
    end
end

%% Choose caseA1 w/ 1) modified Pgen as above 2) w/o changes
% %1)
% caseA = caseA1;
% %2)
caseA = caseA0;

%% All variable parameters
k_load = 2.5; % values k_load < 1 decrease the load

k_gen = 1.3; % values k_gen < 1 decrease the generation

% Maximum number of iterations for each outer problem
MaxIt = 1000;

% Gradient descent step
mu0 = 10^-1; 

% Diminishing step size
DimStep = 0; % if 1 - diminishing step size;

% Threshold for stopping criteria
dFz_threshold = 10^-6;

% Epsilon value
epsilon = 0.01;

% Incremental/Non-incremental (if IncremFlag=1, then incremental; if
% IncremFlag=0, then non-incremental)
IncremFlag = 0;

%% Parameters further do not change between scenarios
% We assume that topology does not change between scenarios

% Use 2 morning hours for deriving control rules
load('141Load_30Pv_24Sc.mat');
% Seen scenarios
PV_coef_TV = Solar30_MorTrain_Sc;
Load_coef_TV = Load141_MorTrain_Sc;
% In-sample scenarios
% PV_coef_TV = Solar30_Mor_InSample_Sc;
% Load_coef_TV = Load141_Mor_InSample_Sc;

% Use 2 day hours for deriving control rules 
% load('141Load_30Pv_24Sc_Day.mat');
% Seen scenarios
% PV_coef_TV = Solar30_DayTrain_Sc;
% Load_coef_TV = Load141_DayTrain_Sc;
% In-sample scenarios
% PV_coef_TV = Solar30_Day_InSample_Sc;
% Load_coef_TV = Load141_Day_InSample_Sc;

% Changing loading of the original system
caseA.bus(:,3) = k_load*caseA.bus(:,3);
caseA.bus(:,4) = k_load*caseA.bus(:,4);

% Changing generation of the original system
caseA.gen(:,[2,3,4,9,10]) = k_gen*caseA.gen(:,[2,3,4,9,10]); % Pg
% Qg is changed via qU

% Changing base power and R, X for comparison with LPF 
caseA.baseMVA = 1; 
caseA.branch(:,3:4) = caseA0.branch(:,3:4)/10; 

% Number of scenarios
T = size(Load_coef_TV,1);

% Runpf option
mpopt = mpoption('verbose',0,'out.all', 0,'pf.enforce_q_lims',1);

% Defining dimensions and indices
[NBusIncSla,NGenIncSla,NBusExcSla,NGenExcSla,IndGenExcSla,IndGenExcSlaOrd] = fIndDim(caseA);

%% For the case when only a subset of inverters are active and rest follow NAS
if FlagSubsetActInv == 1
    IndGenExcSla = [];

    for k = 2:NGenIncSla
        if caseA.gen(k,4) > 0
            IndGenExcSla = [IndGenExcSla; caseA.gen(k,1)];
        end
    end

    IndGenExcSlaOrd = sort(IndGenExcSla);
    NGenExcSla = size(IndGenExcSlaOrd,1);
    NGenIncSla = NGenExcSla+1;
end

%% Defining RX matrices
[RBusExcSla,XBusExcSla,IndGenExcSlaInt,R,X,XColRed] = fBigSmallRX(caseA,IndGenExcSlaOrd);

% Creating vectors
V0ExcSla = ones(NBusExcSla,1);
PgExcSla = zeros(NBusExcSla,1);
QgExcSla = zeros(NBusExcSla,1);

% Reactive power limit: NOW WE CAN SATURATE EARLIER
GenInd = (caseA.gen(:,4) > 0) & (caseA.gen(:,4) < 100);
qU = 0.45*caseA.gen(GenInd,4);

caseA.gen(2:end,4) = 0.45*caseA.gen(2:end,9);
caseA.gen(2:end,5) = -0.45*caseA.gen(2:end,9);

% Reference voltage
Vr = V0ExcSla(IndGenExcSlaInt);
VrExcSla = V0ExcSla;

% Assignment matrices
Eye = eye(NBusExcSla);
Agen = Eye(:,IndGenExcSlaInt);

%% Providing initial z0 by solving optimization problem:

% Keep variables' names from Vassilis for convention
G = NGenExcSla;
qmax = qU;
GEN = IndGenExcSlaInt;
S = T;

% Separate z and a
[z0,a0,sol0] = fz0_SepZA(NGenExcSla,qU,IndGenExcSlaInt,XBusExcSla,epsilon,IncremFlag);

% initialize matrix to maintain rule parameters across all iterations
Z = zeros(4*G,1); 
Z(:,1) = z0; % initial estimate (feasible by construction)

A = zeros(G,1);
A(:,1) = a0;

%% Generating tV (vtilde in different scenarios)
% Define empty matrices
VTildeExcSla = zeros(NBusExcSla,T);
VTilde = zeros(NGenExcSla,T);
TotPL_24Sc = zeros(1,T);
TotPG_24Sc = zeros(1,T);
caseB = caseA; % case for time varying load and generation

for t = 1:T
    caseB.bus(:,3) = Load_coef_TV(t,:)'.*caseA.bus(:,3); % PL
    caseB.bus(:,4) = Load_coef_TV(t,:)'.*caseA.bus(:,4); % QL
    caseB.gen(2:end,2) = PV_coef_TV(t,:)'.*caseA.gen(2:end,2); % Pg
    
    [PgExcSla,QgExcSla,Qg,PcExcSla,QcExcSla,...
    VmPFExcSla,pfCaseA] = fOrigPgQgPcQcVmPF(NBusExcSla,NBusIncSla,caseB,IndGenExcSlaInt,G,mpopt);
    
    % Defining some parameters
    VTildeExcSla(:,t) = RBusExcSla*(PgExcSla - PcExcSla) - XBusExcSla*QcExcSla + V0ExcSla; 
    VTilde(:,t) = VTildeExcSla(IndGenExcSlaInt,t);

    % Total active load and generation
    TotPL_24Sc(t) = sum(caseB.bus(:,3));
    TotPG_24Sc(t) = sum(caseB.gen(2:end,2));    
end

%% Loop for inner problem (Varying PG, PL, QL scenarios)
% Initiliaze lists
obj3_list = [];

% Define empty matrices
Qs_list = zeros(NGenExcSla,T*MaxIt);
Fz_list = zeros(MaxIt,1);
StartSOCP = zeros(MaxIt,1);
EndSOCP = zeros(MaxIt,1);

% disp('All script before outer loop is finished')
% toc
%% 
for i = 1:MaxIt

    % Moved inside loop: at each iteration we start from empty lists
    VExcSla_matrix = []; 
    dVExcSla_dtheta_matrix = [];
    disp(['### Outer problem iteration: ', num2str(i)])
    dFsdz = zeros(4*G,1);
    Fz = 0;

%   disp(['Inner problem (over scenarios): ', num2str(t),' / ',num2str(T)]

    % Assigning corresponding z
    if i == 1
        z = z0;
    else
        z = double(z);
    end
    
    % Inner problem      
    % STEP 1: Given voltage scenarios and rule parameters, find equilibrium voltages and injections
    bv = double(z(1:G));
    delta = double(z(G+1:2*G));
    sigma = double(z(2*G+1:3*G));
    c = double(z(3*G+1:4*G));

if i == 1
    Qs_prev = 0; % provide just zero which will not be used
    StartQP = 0;
    EndQP = 0;
end

%     % Original inner problem (for our approach)
    [Qs,Vs,sol1,Cost,StartQP,EndQP,qU_hat] = fInPr(G,T,VTildeExcSla,c,XBusExcSla,IndGenExcSlaInt,delta,qU,Agen,bv,Vr,i,Qs_prev,sigma,StartQP,EndQP);

    Qs_prev = Qs;

    Qs_list(:,T*(i-1)+1:T*i) = Qs;

    % STEP 2: Compute gradients by differentiating through the control rules

    % New extended derivatives after Aug 5
    dfdv = diag(1./c)*((Vs(GEN,:)-bv*ones(S,1)'-sigma*ones(S,1)') >= 0 - (Vs(GEN,:)-bv*ones(S,1)'-delta*ones(S,1)') >= 0 ...
    - ((-Vs(GEN,:)+bv*ones(S,1)'-delta*ones(S,1)') >= 0) + ((-Vs(GEN,:)+bv*ones(S,1)'-sigma*ones(S,1)') >= 0));
    
    dfdvbar = diag(1./c)*((-(Vs(GEN,:)-bv*ones(S,1)'-sigma*ones(S,1)') >= 0) + (Vs(GEN,:)-bv*ones(S,1)'-delta*ones(S,1)') >= 0 ...
    + ((-Vs(GEN,:)+bv*ones(S,1)'-delta*ones(S,1)') >= 0) - ((-Vs(GEN,:)+bv*ones(S,1)'-sigma*ones(S,1)') >= 0));
    
    dfddelta = diag(1./c)*((Vs(GEN,:)-bv*ones(S,1)'-delta*ones(S,1)') >= 0 - (-Vs(GEN,:)+bv*ones(S,1)'-delta*ones(S,1)') >= 0);
    
    dfdsigma = diag(1./c)*(((-Vs(GEN,:)+bv*ones(S,1)'-sigma*ones(S,1)') >= 0) - (Vs(GEN,:)-bv*ones(S,1)'-sigma*ones(S,1)') >= 0);

    dfdc = -diag(1./c)*Qs;

    for s = 1:S        
        % Modified gradient calculation since August 5
        dfdvs = diag(dfdv(:,s));
        I = eye(G,G);
        dfdz = [diag(dfdvbar(:,s)), diag(dfddelta(:,s)), diag(dfdsigma(:,s)), diag(dfdc(:,s))];
%         dQsdz = inv(I - dfdvs*XBusExcSla(GEN,GEN))*dfdz;
        dQsdz = (I - dfdvs*XBusExcSla(GEN,GEN))\dfdz;
        dvsdz = XBusExcSla*Agen*dQsdz;
        dFsdz = dFsdz + dvsdz'*(Vs(:,s) - ones(NBusExcSla,1))/S;
        Fz = Fz + norm(Vs(:,s) - ones(NBusExcSla,1),2)^2/(2*S);
    end

    disp(['Fz: ', num2str(Fz)])
    g = dFsdz;

    Fz_list(i) = Fz;

    % STEP 3: Gradient update and projection onto mcZ
    if DimStep == 1
        mu = mu0/sqrt(i);
    else
        mu = mu0;
    end
      
    x = Z(:,i) - mu*g; % this is Z(:,i+1) before projection. ILGIZ: eq.(21)
    
    clear('yalmip');
    z = sdpvar(4*G,1);
    a = sdpvar(G,1);
    obj3 = norm(x(1:4*G)-z(1:4*G),2)^2; % Only z variables
    
    % z(1:G) -> v_bar
    con = (0.95*ones(G,1) <= z(1:G));
    con = con + (z(1:G) <= 1.05*ones(G,1));
    
    % z(G+1:2*G) -> delta
    con = con + (zeros(G,1) <= z(G+1:2*G));
    con = con + (z(G+1:2*G) <= 0.03*ones(G,1));
    
    % z(2*G+1:3*G) -> sigma
    con = con + (z(G+1:2*G) + 0.02*ones(G,1) <= z(2*G+1:3*G));
    con = con + (z(2*G+1:3*G) <= 0.18*ones(G,1));
    
    con = con + (z(2*G+1:3*G)-z(G+1:2*G) <= (diag(qmax))*z(3*G+1:4*G)); % z(3*G+1:4*G) -> с
    con = con + (z(3*G+1:4*G) >= sum(XBusExcSla(GEN,GEN))'/(1-epsilon)); 
    if IncremFlag == 0
        con = con + ((XBusExcSla(GEN,GEN))*a <= (1-epsilon)*ones(G,1)); 
        con = con + (diag(a)*z(3*G+1:4*G) >= ones(G,1));
    end
    settings = sdpsettings('verbose',0,'solver','gurobi');

%     disp('Started solving SOCP')
    StartSOCP(i) = toc;

    sol3 = optimize(con,obj3,settings);

%     disp('Completed solving SOCP')
    EndSOCP(i) = toc;

    if ~sol3.problem % checking if solver had any issues
        Z(:,i+1) = double(z);
        A(:,i+1) = double(a);
        obj3_list = [obj3_list;double(obj3)];
    end

    if i > 2
        if abs(Fz_list(i) - Fz_list(i-1))/Fz_list(i-1) < dFz_threshold
            break
        end
    end
   
end

%% Computing times
TimeQP = -StartQP + EndQP;
TimeSOCP = -StartSOCP + EndSOCP;
AvTimeQP = mean(TimeQP)
AvTimeSOCP = mean(TimeSOCP)
AvTimeSQPsNSOCP = S*AvTimeQP+AvTimeSOCP

%% Plots
close all;

figure(1)
plot(z0)
title('starting z0')
xlabel('ID of z0 variable')
ylabel('z0 values')

figure(2)
plot(Qs)
title('Qs at the final iteration')
xlabel('Inverter ID')
ylabel('Qs value')

figure(3)
plot(Vs)
title('Vs at the final iteration')
xlabel('Bus ID')
ylabel('Vs value')

figure(4)
plot(Z(1:G,1:i)')
title('Vbar convergence')
xlabel('iteration')
ylabel('Vbar values, p.u.')

figure(5)
plot(Z(G+1:2*G,1:i)')
title('delta convergence')
xlabel('iteration')
ylabel('delta values')

figure(6)
plot(Z(2*G+1:3*G,1:i)')
title('sigma convergence')
xlabel('iteration')
ylabel('sigma values')

figure(7)
plot(Z(3*G+1:4*G,1:i)')
title('c convergence')
xlabel('iteration')
ylabel('c values')

figure(8)
plot(A(:,1:i)')
title('a convergence')
xlabel('iteration')
ylabel('a values')

figure(9)
plot(Fz_list)
title('convergence of cost')
xlabel('iteration')
ylabel('convergence of cost, p.u.')

% disp('Plotting is finished')
% toc

%% save the simulation results
