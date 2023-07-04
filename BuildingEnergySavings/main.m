% By HongyiLi.
% Thanks to upenn!
%% Create an mlepProcess instance and configure it

ep = mlepProcess;
%ep.arguments = {'SmOffPSZ', 'USA_IL_Chicago-OHare.Intl.AP.725300_TMY3'};
ep.arguments = {'exp', 'CHN_Guangdong.Shenzhen.594930_SWERA'};
%ep.arguments = {'SmOffPSZ', 'CHN_Guangdong.Shenzhen.594930_SWERA'};
ep.acceptTimeout = 6000;

VERNUMBER = 2;  % version number of communication protocol (2 for E+ 7.2.0)

%% Start EnergyPlus cosimulation
[status, msg] = ep.start;  

if status ~= 0
    error('Could not start EnergyPlus: %s.', msg);
end

[status, msg] = ep.acceptSocket;

if status ~= 0
    error('Could not connect to EnergyPlus: %s.', msg);
end

%% 初始化模型
load linear01;
%parameters of the Building Model
A=aAdsun1;
Bu =aBdsun1;
Bd =aDdsun1;
C=zeros(1,9);
%other parameters
nx = length(A); %number of states
ny = size(C,1); %number of outputs
nu = size(Bu,2); %number of inputs
nd = size(Bd,2); %number of disturbances

%define simulation parameters
N = 32;
%T = 504;

%define objective parameters
yrefl=20;
yref=26;
umax = 1500*ones(nu,1); umin = 0*ones(nu,1); %input constraints
ymax = 24*ones(ny,1); ymin = 22*ones(ny,1); %output constraints
%define constraints matrix
param_constraints.Mu = [eye(nu); -eye(nu)];
param_constraints.My = [eye(ny); -eye(ny)];
param_constraints.mu = [umax; -umin];
param_constraints.my = [ymax; -ymin];

Mu = param_constraints.Mu;
My = param_constraints.My;
mu = param_constraints.mu;
my = param_constraints.my;

%if strcmp(cmd,'normal')
d = [20;80;0]; %disturbances
d(2)=80;
d(3)=0;

%schedual
morning=9;
night=19;
%% The main simulation loop

deltaT = 15*60;  % time step = 15 minutes
kStep = 1;  % current simulation step
MAXSTEPS = 8*24*4+1;  % max simulation time = 2 days
yrefh =zeros(MAXSTEPS+N, 1);
my = zeros(MAXSTEPS+N, 2);
%MAXSTEPS=20;
TCRooLow = 22;  % Zone temperature is kept between TCRooLow & TCRooHi
TCRooHi = 24;
TOutLow = 25;  % Low level of outdoor temperature
TOutHi = 32;  % High level of outdoor temperature
ratio = (TCRooHi - TCRooLow)/(TOutHi - TOutLow);

% logdata stores set-points, outdoor temperature, and zone temperature at
% each time step.
logdata = zeros(MAXSTEPS, 24);
logdataa = zeros(MAXSTEPS, 4);
logdatau=zeros(MAXSTEPS, N);
logdataPMVvaule=zeros(MAXSTEPS, N);
logdataTclvaule=zeros(MAXSTEPS, N);
logdataTmrvaule=zeros(MAXSTEPS, N);
logdataTavaule=zeros(MAXSTEPS, N);
logdataobjective=zeros(MAXSTEPS, 1);
logdatav=zeros(MAXSTEPS, N);
logdatay=zeros(MAXSTEPS, N);
logdatamz=zeros(MAXSTEPS, 1);
logdatacool=zeros(MAXSTEPS, 1);
logdatamc=zeros(MAXSTEPS, 1);
logdataset=zeros(MAXSTEPS, 1);
logdatapreset=zeros(MAXSTEPS, 1);
comfort=zeros(MAXSTEPS, 1);
%
logdatado=zeros(MAXSTEPS, 1);
logdatadl=zeros(MAXSTEPS, 1);
logdatade=zeros(MAXSTEPS, 1);
logdatads=zeros(MAXSTEPS, 1);
logdatads2=zeros(MAXSTEPS, 1);
%}
logdataPMV=zeros(MAXSTEPS, 1);
logdatacount=zeros(MAXSTEPS, 1);
logdataRH=zeros(MAXSTEPS, 1);
logdatashiwenpre=zeros(MAXSTEPS, 1);
logdatanengliangpre=zeros(MAXSTEPS, 1);
schedual=zeros(MAXSTEPS+N, 1);
iv=zeros(MAXSTEPS, 1);
fee=zeros(MAXSTEPS+N, 1);
for j=1:MAXSTEPS+N
    i=mod(j,4*24);
    if (i>4*0)&&(i<4*8+2)
        fee(j)=33.58*0.25*0.001;
        yrefh(j)=23;
        my(j,:)=[0.1 0.1];
    end
    if (i>4*8+1)&&(i<4*14+2)
        fee(j)=66.29*0.25*0.001;
        yrefh(j)=24;
        my(j,:)=[0.6 0.6];
    end
    if (i>4*14+1)&&(i<4*17+2)
        fee(j)=108.81*0.25*0.001;
        yrefh(j)=25;
        my(j,:)=[0.7 0.7];
    end
    if (i>4*17+1)&&(i<4*19+2)
        fee(j)=66.29*0.25*0.001;
        yrefh(j)=24;
        my(j,:)=[0.6 0.6];
    end
    if (i>4*19+1)&&(i<4*22+2)
        fee(j)=108.81*0.25*0.001;
        yrefh(j)=25;
        my(j,:)=[0.7 0.7];
    end
    if (i>4*22+1)&&(i<4*24+2)
        fee(j)=66.29*0.25*0.001;
        yrefh(j)=24;
        my(j,:)=[0.6 0.6];
    end
    if (i>4*24+1)&&(i<MAXSTEPS+N+1)
        fee(j)=33.58*0.25*0.001;
        yrefh(j)=23;
        my(j,:)=[0.1 0.1];
    end   
    j=j+1;
end
for j=1:MAXSTEPS
    i=mod(j,4*24);
    if (i>=1+4*morning)&&(i<=1+4*night)
        schedual(j)=1;
    end
    j=j+1;
end
%PMV
M=60;
W=0;
%Ta=25;
%RH=37;
%Pa=RH*6.1094*exp(17.625*Ta/(Ta+243.04))*10e-4;
%Tmr=28.4;
Icl=0.155*0.5;
fcl=1+1.29*Icl;%
Var=0.1;
hc=12.1*sqrt(Var);%
a=35.7-0.0275*(M-W);
b=Icl*fcl*hc;
c=Icl*3.96*10e-9*fcl;
%Tcl=35.7-0.0275*(M-W)-Icl*(3.96*10e-9*fcl*((Tcl+273)^4-(Tmr+273)^4))-Icl*fcl*hc*(Tcl-Ta);%

%PMV=(0.303*exp(-0.036*M)+0.0275)*(M-W-3.05*(5.733-0.007*(M-W)-Pa)...
%-0.42*(M-W-58.2)-0.0173*M*(5.867-Pa)-0.0014*M*(34-Ta)-3.96*10e-9*fcl*((Tcl+273)^4-(Tmr+273)^4)-fcl*hc*(Tcl-Ta));



Tcl0=30;
Tmr0=27;
Ta0=25;
for i = 1:50
Tcl0=35.7-0.0275*(60-0)-0.0775*(3.96*10e-9*1.1*((Tcl0+273)^4-(Tmr0+273)^4))-0.0775*1.1*3.8264*(Tcl0-Ta0);
end
a=35.7-0.0275*(M-W);
b=Icl*fcl*hc;
c=Icl*3.96*10e-9*fcl;
dTa=b/(1+b+4*c*(Tcl0+273)^3);
dTcl=4*c*(Tmr0+273)^3/(1+b+4*c*(Tcl0+273)^3);

lattice = rPWLtoLattice(PWLS,Z);
tic;
while kStep <= MAXSTEPS    
    % Read a data packet from E+
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % Parse it to obtain building outputs
    [flag, eptime, outputs] = mlepDecodePacket(packet);
    if flag ~= 0, break; end
        
    % BEGIN Compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day
    
    %x0=[0 0 0 0 0];
    totalcooling=0;

    x0(1,1)=outputs(1);
    x0(2,1)=outputs(2);
    x0(3,1)=outputs(3);
    x0(4,1)=outputs(4);
    x0(5,1)=outputs(5);
    x0(6,1)=outputs(6);
    x0(7,1)=outputs(7);
    x0(8,1)=outputs(8);
    x0(9,1)=outputs(9);
    relativehum=outputs(22);
    d(1)=outputs(10);
    d(2)=outputs(11);
    d(3)=outputs(12);
    d(4)=outputs(19);
    d(5)=outputs(20);
    d(6)=0;
    totalcooling=outputs(13);
    mz=outputs(14);
    rain=outputs(17);
    windspeed=outputs(18);
    %d=eplus_out_prev.dist;
    logdatado(kStep, :)=outputs(10);
    logdatadl(kStep, :)=outputs(11);
    logdatade(kStep, :)=outputs(12);
    logdatads(kStep, :)=outputs(19);
    logdatads2(kStep, :)=outputs(20);
    logdataPMV(kStep, :)=outputs(21);
    logdataRH(kStep, :)=outputs(22);
    u0 = 1000* rand(1,N);
    uintal=u0;
    ugood =5000*ones(1,N);
    count=0;
    while count<5 && ((uintal(1:N-1)-ugood(1:N-1))*(uintal(1:N-1)-ugood(1:N-1))'>1*N )
        d(1)=logdatado(kStep);
        d(2)=logdatadl(kStep);
        d(3)=logdatade(kStep);
        d(4)=0.5*logdatads(kStep);
        d(5)=0.5*logdatads2(kStep);
        x1(:,1) = x0;
        x1(:,2) = A*x0 - Bu*u0(1) + Bd*d;
        TMR0(:,1)=(x1(2,1)+x1(3,1)+x1(4,1)+x1(5,1))/4;
        
        for i=2:N-1
            d(1)=logdatado(kStep+i-1);
            d(2)=logdatadl(kStep+i-1);
            d(3)=logdatade(kStep+i-1);
            d(4)=0.5*logdatads(kStep+i-1);
            d(5)=0.5*logdatads2(kStep+i-1);
            x1(:,i+1) = A*x1(:,i) - Bu*u0(i) + Bd*d;
            TMR0(:,i)=(x1(2,i)+x1(3,i)+x1(4,i)+x1(5,i))/4;
        end
        for i=1:N
            z0(:,i) = [x1(1,i);(x1(2,i)+x1(3,i)+x1(4,i)+x1(5,i))/4] ;
            [Aw(:,i),Bw(i)]=rLatticeValue(lattice,z0(:,i));
        end
        
        TCL0(:,1)=Aw(1,1)*x1(1,1)+Aw(2,1)*TMR0(:,1)+Aw(3,1);
        Pa0(:,1)=relativehum*(1.78333*x1(1,1)-12.75160)*10e-4;
        %PMV0(:,1)=0.02265*x1(1,1)+0.80516*TCL0(:,1)-25.41573;
        PMV0(:,1)=0.0624*(4.088*Pa0(:,1)+0.084*x1(1,1)+12.9032*TCL0(:,1)-405.2614);

        
        for i=2:N-1
            TCL0(:,i)=Aw(1,i)*x1(1,i)+Aw(2,i)*TMR0(:,i)+Aw(3,i);
            %PMV0(:,i)=0.02265*x1(1,i)+0.80516*TCL0(:,i)-25.41573;
            Pa0(:,i)=relativehum*(1.78333*x1(1,i)-12.75160)*10e-4;
            PMV0(:,i)=0.0624*(4.088*Pa0(:,i)+0.084*x1(1,i)+12.9032*TCL0(:,i)-405.2614);

        end
            
        ob=0;
        for i=1:N-1
            if schedual(kStep+i-1)==1
                ob = ob+ fee(kStep+i-1)*(u0(i))^2 + 25000*(PMV0(:,i))^2;
            else
                if i<N-1
                    ob = ob + fee(kStep+i-1)*(u0(i))^2+(u0(i+1)-u0(i))^2 ;
                end
            end          
        end
        %求解
        %define optimization variables
        yalmip('clear');
        x = sdpvar(nx, N, 'full'); %states
        x0 = sdpvar(nx, 1, 'full'); %initial states
        y = sdpvar(ny, N, 'full'); %outputs
        u = sdpvar(nu, N, 'full'); %inputs
        d = sdpvar(nd, 1, 'full');
        Ta = sdpvar(1, N, 'full');
        Tmr = sdpvar(1, N, 'full');
        Pa = sdpvar(1, N, 'full');
        Tcl = sdpvar(1, N, 'full');
        PMV = sdpvar(1, N, 'full');
        Tx = sdpvar(1, N, 'full');
        z0 = sdpvar(2, 1, 'full');
        A1 = sdpvar(1, 3, 'full');
        B1 = sdpvar(1, 1, 'full');
        x0(1)=outputs(1);
        x0(2)=outputs(2);
        x0(3)=outputs(3);
        x0(4)=outputs(4);
        x0(5)=outputs(5);
        x0(6)=outputs(6);
        x0(7)=outputs(7);
        x0(8)=outputs(8);
        x0(9)=outputs(9);
        d(1)=outputs(10);
        d(2)=outputs(11);
        d(3)=outputs(12);
        d(4)=outputs(19);
        d(5)=outputs(20);
        d(6)=0;
        relativehum=outputs(22);
        %define constraints and objective
        
        constraints = [];
        objective = 0;
        for i = 1:N-1
            if i == 1 
                d(1)=logdatado(kStep+i-1);
                d(2)=logdatadl(kStep+i-1);
                d(3)=logdatade(kStep+i-1);
                d(4)=0.5*logdatads(kStep+i-1);
                d(5)=0.5*logdatads2(kStep+i-1);
                constraints = [constraints, x(:,i) == x0];
                constraints = [constraints, x(:,i+1) == A*x0 - Bu*u(:,i) + Bd*d]; %system dynamics - first step
                constraints = [constraints,Ta(:,i)==x(1,i)];
                constraints = [constraints,Tmr(:,i)==(x(2,i)+x(3,i)+x(4,i)+x(5,i))/4];
                constraints = [constraints, Pa(:,i)==relativehum*(1.78333*Ta(:,i)-12.75160)*10e-4];
                %constraints = [constraints,Tcl(:,i)==Aw(1,i)*x(1,i)+Aw(2,i)*Tmr(:,i)+Aw(3,i)];
                %constraints = [constraints,PMV(:,i)==0.02265*Ta(:,i)+0.80516*Tcl(:,i)-25.41573];
            else
                if kStep<MAXSTEPS -N
                    d(1)=logdatado(kStep+i-1);
                    d(2)=logdatadl(kStep+i-1);
                    d(3)=logdatade(kStep+i-1);
                    d(4)=0.5*logdatads(kStep+i-1);
                    d(5)=0.5*logdatads2(kStep+i-1);
                end
                constraints = [constraints, x(:,i+1) == A*x(:,i) - Bu*u(:,i) + Bd*d]; %system dynamics
                constraints = [constraints,Ta(:,i)==x(1,i)];
                constraints = [constraints,Tmr(:,i)==(x(2,i)+x(3,i)+x(4,i)+x(5,i))/4];
                constraints = [constraints, Pa(:,i)==relativehum*(1.78333*Ta(:,i)-12.75160)*10e-4];
                %constraints = [constraints,Tcl(:,i)==Aw(1,i)*x(1,i)+Aw(2,i)*Tmr(:,i)+Aw(3,i)];
                %constraints = [constraints,PMV(:,i)==0.02265*Ta(:,i)+0.80516*Tcl(:,i)-25.41573];
            end
            constraints = [constraints, Mu*u(:,i) <= mu]; %input constraints        
            if schedual(kStep+i-1)==1
                %objective = objective + fee(kStep+i-1)*(u(:,i))^2 + 25000*(0.02265*Ta(:,i)+0.80516*(Aw(1,i)*x(1,i)+Aw(2,i)*Tmr(:,i)+Aw(3,i))-25.41573)^2;
                objective = objective + fee(kStep+i-1)*(u(:,i))^2 + 25000*(0.0624*(4.088*Pa(:,i)+0.084*Ta(:,i)+12.9032*(Aw(1,i)*x(1,i)+Aw(2,i)*Tmr(:,i)+Aw(3,i))-405.2614))^2;
            else
                if i<N-1
                    objective = objective + fee(kStep+i-1)*(u(:,i))^2+(u(:,i+1)-u(:,i))^2 ;
                end
            end          
        end
        ops = sdpsettings('verbose', 1);
        result = optimize(constraints, objective, ops);
        ugood =value(u);
        uintal=u0;
        disp(uintal);
        disp(ugood);
        if  ob>value(objective)
            u0 = ugood;
        else
            up =  (u0+value(u))/2;
            ob1=999999999999;
            count_go=0;
            while (ob1>ob+1000)&&(count_go<5)
                d(1)=logdatado(kStep);
                d(2)=logdatadl(kStep);
                d(3)=logdatade(kStep);
                d(4)=0.5*logdatads(kStep);
                d(5)=0.5*logdatads2(kStep);
                count_go=count_go+1;
                disp(count_go);
                x1(:,1) = x0;
                x1(:,2) = A*x0 - Bu*up(1) + Bd*d;
                TMR0(:,1)=(x1(2,1)+x1(3,1)+x1(4,1)+x1(5,1))/4;
        
                for i=2:N-1
                    d(1)=logdatado(kStep+i-1);
                    d(2)=logdatadl(kStep+i-1);
                    d(3)=logdatade(kStep+i-1);
                    d(4)=0.5*logdatads(kStep+i-1);
                    d(5)=0.5*logdatads2(kStep+i-1);
                    x1(:,i+1) = A*x1(:,i) - Bu*up(i) + Bd*d;
                    TMR0(:,i)=(x1(2,i)+x1(3,i)+x1(4,i)+x1(5,i))/4;
                end
   
                for i=1:N
                    z0(:,i) = [x1(1,i);(x1(2,i)+x1(3,i)+x1(4,i)+x1(5,i))/4] ;
                    [Aw(:,i),Bw(i)]=rLatticeValue(lattice,z0(:,i));
                end
                TCL0(:,1)=Aw(1,1)*x(1,1)+Aw(2,1)*TMR0(:,1)+Aw(3,1);
                %PMV0(:,1)=0.02265*x1(1,1)+0.80516*TCL0(:,1)-25.41573;
                Pa0(:,1)=relativehum*(1.78333*x1(1,1)-12.75160)*10e-4;
                PMV0(:,1)=0.0624*(4.088*Pa0(:,1)+0.084*x1(1,1)+12.9032*TCL0(:,1)-405.2614);

        
                for i=2:N-1
                    TCL0(:,i)=Aw(1,i)*x(1,i)+Aw(2,i)*TMR0(:,i)+Aw(3,i);
                    %PMV0(:,i)=0.02265*x1(1,i)+0.80516*TCL0(:,i)-25.41573;
                    Pa0(:,i)=relativehum*(1.78333*x1(1,i)-12.75160)*10e-4;
                    PMV0(:,i)=0.0624*(4.088*Pa0(:,i)+0.084*x1(1,i)+12.9032*TCL0(:,i)-405.2614);
                end

                ob1=0;
                for i=1:N-1
                    if schedual(kStep+i-1)==1
                        ob1 = ob1+ fee(kStep+i-1)*(up(i))^2 + 25000*(PMV0(:,i))^2;
                    else
                        if i<N-1
                            ob1 = ob1 + fee(kStep+i-1)*(up(i))^2+(up(i+1)-up(i))^2 ;
                        end
                    end          
                end
                count=count+0.001;
                up =  (u0+up)/2;
                disp(up);
                disp(u0);
                errx = ob1-ob;
                disp(ob1);
                disp(ob);
            end
            u0=up;
            uintal=u0;
            count=count+0.01;
        end
        
        count=count+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    u0=value(u);
        d(1)=logdatado(kStep);
        d(2)=logdatadl(kStep);
        d(3)=logdatade(kStep);
        d(4)=0.5*logdatads(kStep);
        d(5)=0.5*logdatads2(kStep);
        x1(:,1) = x0;
        x1(:,2) = A*x0 - Bu*u0(1) + Bd*d;
        TMR0(:,1)=(x1(2,1)+x1(3,1)+x1(4,1)+x1(5,1))/4;
        
        for i=2:N-1
            d(1)=logdatado(kStep+i-1);
            d(2)=logdatadl(kStep+i-1);
            d(3)=logdatade(kStep+i-1);
            d(4)=0.5*logdatads(kStep+i-1);
            d(5)=0.5*logdatads2(kStep+i-1);
            x1(:,i+1) = A*x1(:,i) - Bu*u0(i) + Bd*d;
            TMR0(:,i)=(x1(2,i)+x1(3,i)+x1(4,i)+x1(5,i))/4;
        end
        for i=1:N-1
            TCL0(:,i)=30;
            for j=1:30
                TCL0(:,i)=35.7-0.0275*(M-W)-Icl*(3.96*10e-9*fcl*((TCL0(:,i)+273)^4-(TMR0(:,i)+273)^4))-Icl*fcl*hc*(TCL0(:,i)-x1(1,i));
            end
        end
        %PMV0(:,1)=0.02265*x1(1,1)+0.80516*TCL0(:,1)-25.41573;
        %Pa0(:,1)=relativehum*(1.78333*x1(1,1)-12.75160)*10e-4;
        Pa0(:,1)=relativehum*6.1094*exp(17.625*x1(1,1)/(x1(1,1)+243.04))*10e-4;
        PMV0(:,1)=0.0624*(4.088*Pa0(:,1)+0.084*x1(1,1)+12.9032*TCL0(:,1)-405.2614);

        
        for i=2:N-1
            %PMV0(:,i)=0.02265*x1(1,i)+0.80516*TCL0(:,i)-25.41573;
            %Pa0(:,i)=relativehum*(1.78333*x1(1,i)-12.75160)*10e-4;
            Pa0(:,i)=relativehum*6.1094*exp(17.625*x1(1,i)/(x1(1,i)+243.04))*10e-4;
            PMV0(:,i)=0.0624*(4.088*Pa0(:,i)+0.084*x1(1,i)+12.9032*TCL0(:,i)-405.2614);
        end
            
        obb=0;
        for i=1:N-1
            if schedual(kStep+i-1)==1
                obb = obb+ fee(kStep+i-1)*(u0(i))^2 + 25000*(PMV0(:,i))^2;
            else
                if i<N-1
                    obb = obb + fee(kStep+i-1)*(u0(i))^2+(u0(i+1)-u0(i))^2 ;
                end
            end          
        end
        logdatacostoo(kStep)=obb;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(count);
    logdatacount(kStep,1)=count;
    d(1)=outputs(10);
    d(2)=outputs(11);
    d(3)=outputs(12);
    d(4)=0.5*outputs(19);
    d(5)=0.5*outputs(20);
    d(6)=0;
    uo=ugood(1);
    sp1=value(A*x0 - Bu*uo + Bd*d);
    setpoint(kStep)=value(sp1(1));
    if uo<1
        settt=[35 2 3 4 5 6 7 8 9 10 11 12 13];
    else
        settt = value(A*x0 - Bu*uo + Bd*d);
    end
    %controller = optimizer(constraints, objective, ops, parameters_in, solutions_out);
    %simulate the system
    %option = 1; %simulation with no night-setbacks and no variable cost
    %[xt, yt, ut, t, ~] = simBuild(controller, T, @shiftPred, N, option);
    logdatau(kStep, :)=value(u);
    logdataPMVvaule(kStep, :)=value(PMV);
    logdataTclvaule(kStep, :)=value(Tcl);
    logdataTmrvaule(kStep, :)=value(Tmr);
    logdataTavaule(kStep, :)=value(Ta);
    logdataobjective(kStep, :)=value(objective);
    logdatay(kStep, :)=value(y);
    logdatamz(kStep, :)=mz;
    logdatacool(kStep, :)=totalcooling;
    logdatamc(kStep, :)=totalcooling/(1007*mz);
    uone=value(u(:,1)); %取第一个值送到mpc
    if mz==0
        preset=x0(1);
        tttset=x0(1);
    else
        preset=(uone/(1007*mz))+13;
        tttset=(totalcooling/(1007*mz))+13;
    end
    logdataset(kStep, :)=tttset;
    logdatapreset(kStep, :)=preset;
    %tttset=25;
    %ttt = result(1);
    %SP = [10 settt(1)];
    if (dayTime >= 0*3600) && (dayTime <= 23*3600)
            % It is day time (6AM-6PM)
        
            % The Heating set-point: day -> 20, night -> 16
            % The cooling set-point is bounded by TCRooLow and TCRooHi
       SP = [10 settt(1)];
      % SP = [10 22+8*rand(1)];
     %   SP = [10, max(TCRooLow, ...
                  %  min(TCRooHi, TCRooLow + (outputs(10) - TOutLow)*ratio))];
        %SP = [10 22];
    else
            % The Heating set-point: day -> 20, night -> 16
            % The Cooling set-point: night -> 30
        %SP = [10 settt(1)];
        %Sp = [10 22+3*rand(1)]
       SP = [10 35];
    end
    %simple
    %{
    SP = [10 36];
    if (dayTime >= 10*3600) && (dayTime <= 14*3600)
        SP = [10 24];
    end
    if (dayTime >= 17*3600) && (dayTime <= 19*3600)
        SP = [10 24];
    end
    if (dayTime >= 14*3600) && (dayTime <= 17*3600) 
           SP = [10 25];
    end
    if  (dayTime >= 19*3600) && (dayTime <= 20*3600)
           SP = [10 25];
    end   
   %}
    if (dayTime >= morning*3600) && (dayTime <= night*3600)
        if outputs(1)>yrefh
            comfort(kStep,:)=outputs(1)-yrefh(kStep);
        end
    end
        % END Compute next set-points
    
        % Write to inputs of E+
        logdatacost(kStep)=value(objective);
    ep.write(mlepEncodeRealData(VERNUMBER, 0, (kStep-1)*deltaT, SP));    
        % Save to logdata
    logdata(kStep, :) = [SP outputs];
    logdataa(kStep, :) = [SP outputs(10) outputs(1)];
    kStep = kStep + 1;
    disp(['#################################################',num2str(kStep)]);
end
toc;
icomf=0;
money=0;
elemoney=0;
iiv=0;
iPMV=0;
for i=96:MAXSTEPS 
    icomf=icomf+comfort(i,:);
    iiv=iiv+iv(i,:);
    money=money+logdatacool(i,:);
    elemoney=elemoney+fee(i)*logdatacool(i,:);
end
for i=96:MAXSTEPS
    if schedual(i)==1
        iPMV=iPMV+abs(logdataPMV(i));
    end
end
iPMV=iPMV/(41*7);
% Stop EnergyPlus
ep.stop;

disp(['Stopped with flag ' num2str(flag)]);

% Remove unused entries in logdata
kStep = kStep - 1;
if kStep < MAXSTEPS
    logdata((kStep+1):end,:) = [];
end
%{
figure(1);
plot([0:(kStep-2)]'*deltaT/3600, logdatashiwenpre(2:97,1));
hold on;
plot([0:(kStep-2)]'*deltaT/3600, logdatashiwenpre(1:96,2));
legend('actual', 'predictive');
xlabel('Time (hour)');
ylabel('Temperature (\circC)');
figure(2);
plot([0:(kStep-2)]'*deltaT/3600, logdatanengliangpre(2:97,1));
hold on;
plot([0:(kStep-2)]'*deltaT/3600, logdatanengliangpre(1:96,2));
legend('actual', 'predictive');
xlabel('Time (hour)');
ylabel('Power (W)');
% Plot results
%}
figure(3);
plot([0:(kStep-1)]'*deltaT/3600, logdataa);
%plot(([96:(kStep-1)]'*deltaT/3600-24), logdataa(96:(kStep-1),:));
legend('Heat SP', 'Cool SP', 'Outdoor', 'Zone');
%title('Temperatures');
xlabel('Time (hour)');
ylabel('Temperature (\circC)');
grid minor;
disp(['time: ',num2str(toc)]);
disp(['uncomfort index: ',num2str(icomf)]);
disp(['energy index: ',num2str(money)]);
disp(['electricity index: ',num2str(elemoney)]);
disp(['PMV: ',num2str(iPMV)]);
logdatacostoo=logdatacostoo';
%plot([0:(kStep-1)]'*deltaT/3600, logdatacool);
%title('Energy');
%xlabel('Time (hour)');
%ylabel('power(W)');
% ==========FLAGS==============
% Flag	Description
% +1	Simulation reached end time.
% 0	    Normal operation.
% -1	Simulation terminated due to an unspecified error.
% -10	Simulation terminated due to an error during the initialization.
% -20	Simulation terminated due to an error during the time integration.