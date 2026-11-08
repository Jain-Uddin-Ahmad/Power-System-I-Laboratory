clc;
clear all;
close all;

tic
%% variables of System Modeling (Time-domain Simulation)
global Vref ws PC Yabs Yang m n PL QL a0 x0 V0 alpha beta H Xd Xdp Xdpp Xq Xqp Xqpp Td0p Td0pp Tq0p Tq0pp Rs Xls Dm KA TA KE TE KF TF Ax Bx TCH TSV RD MH bb1 bb2 aa1 aa2  x01 x02 x03 f
m = 3; n = 9;
baseMVA = 100; 
ws = 2*pi*50;
ws_vector = ws*ones(1,m);  
w0=ws_vector;

mpopt = mpoption('PF_ALG', 1,'ENFORCE_Q_LIMS',0);
[RESULTS1, SUCCESS] =runpf('IEEE9Bus',mpopt);
mpc=IEEE9Bus;
bus1=RESULTS1.bus;branch1=RESULTS1.branch;gen1=RESULTS1.gen;
Ybus1 = makeYbus(baseMVA, bus1, branch1) ;Ybus1=full(Ybus1);

Ybus=Ybus1; Y = Ybus; bus=bus1; branch=branch1; gen=gen1;
Yabs = abs(Y); Yang = angle(Y);
IC1 = bus(:,8);
IC2 = bus(:,9); gen0 = zeros(n,1); gen0(1:m) = gen(1:m,2); genP = gen0; IC3 = genP;
IC3 = IC3/baseMVA; gen0 = zeros(n,1); gen0(1:m) = gen(1:m,3); genQ = gen0; genQ = genQ+bus(:,6); IC4 = genQ;
IC4 = IC4/baseMVA;
IC5 = bus(:,3); IC5 = IC5/baseMVA;
IC6 = bus(:,4); IC6 = IC6/baseMVA;
IC  = [ IC1 IC2 IC3 IC4 IC5 IC6];
PL  = IC(:,5); QL = IC(:,6);
PG  = IC(:,3); QG = IC(:,4);
TH0 = IC(:,2)*pi/180; TH0 = TH0'; V0 = IC(:,1); V0 = V0';
VG0 = V0(1:m); THG0 = TH0(1:m);
%% Dynamic Paprameters
MD=[ % machine data
    23.640                 6.4000                3.0100;            % 1 -  H
    0.1460                 0.8958                1.3125;            % 2 -  Xd
    0.0608                 0.1198                0.1813;            % 3 -  Xdp
    0.0489                 0.0881                0.1133;            % 4 -  Xdpp
    0.0969                 0.8645                1.2578;            % 5 -  Xq
    0.0969                 0.1969                0.2500;            % 6 -  Xqp
    0.0396                 0.0887                0.0833;            % 7 -  Xqpp
    8.9600                 6.0000                5.8900;            % 8 -  Tdop
    0.1150                 0.0337                0.0420;            % 9 -  Td0pp
    0.3100                 0.5350                0.6000;            % 10-  Tqop
    0.0330                 0.0780                0.1875;            % 11-  Tq0pp
    0.0041                 0.0026                0.0035;            % 12-  RS
    0.1200                 0.1020                0.0750;            % 13-  Xls
    0.1*(2*23.64)/ws       0.2*(2*6.4)/ws      0.3*(2*3.01)/ws      % 14-  Dm
    ];
ED=[ % excitation data
    20.000 *ones(1,m);       % 1- KA
    0.2000 *ones(1,m);       % 2- TA
    1.0000 *ones(1,m);       % 3- KE
    0.3140 *ones(1,m);       % 4- TE
    0.0630 *ones(1,m);       % 5- KF
    0.3500 *ones(1,m);       % 6- TF
    0.0039 *ones(1,m);       % 7- Ax
    1.5550 *ones(1,m);       % 8- Bx
    ];
TD=[ %turbine data
    0.10 *ones(1,m);         % 1- TCH
    0.05 *ones(1,m);         % 2- TSV
    0.05 *ones(1,m);         % 3- RD
    ];
H     = MD(1,:);
Xd    = MD(2,:);
Xdp   = MD(3,:);
Xdpp  = MD(4,:);
Xq    = MD(5,:);
Xqp   = MD(6,:);
Xqpp  = MD(7,:);
Td0p  = MD(8,:);
Td0pp = MD(9,:);
Tq0p  = MD(10,:);
Tq0pp = MD(11,:);
Rs    = MD(12,:);
Xls   = MD(13,:);
Dm    = MD(14,:);
KA    = ED(1,:);
TA    = ED(2,:);
KE    = ED(3,:);
TE    = ED(4,:);
KF    = ED(5,:);
TF    = ED(6,:);
Ax    = ED(7,:);
Bx    = ED(8,:);
TCH   = TD(1,:);
TSV   = TD(2,:);
RD    = TD(3,:);
MH=2*H./ws;
%% Initial Values Calculation (State Variables)
Vphasor  =  VG0.*exp(1i*(THG0));
Iphasor  =  conj((PG(1:m)'+1i*QG(1:m)')./Vphasor);
E0       =  Vphasor+(Rs+1i*Xq).*Iphasor;
Em       =  abs(E0);
D0       =  angle(E0);
Id0      =  real(Iphasor.*exp(-1i*(D0-pi/2)));
Iq0      =  imag(Iphasor.*exp(-1i*(D0-pi/2)));
Edp0     = (Xq-Xqp).*Iq0;
Si2q0    = (Xls-Xq).*Iq0;
Eqp0     =  Rs.*Iq0+Xdp.*Id0+V0(1:m).*cos(D0-TH0(1:m));
Si1d0    =  Eqp0-(Xdp-Xls).*Id0;
Efd0     =  Eqp0+(Xd-Xdp).*Id0;
TM0      = ((Xdpp-Xls)./(Xdp-Xls)).*Eqp0.*Iq0+((Xdp-Xdpp)./(Xdp-Xls)).*Si1d0.*Iq0+((Xqpp-Xls)./(Xqp-Xls)).*Edp0.*Id0-((Xqp-Xqpp)./(Xqp-Xls)).*Si2q0.*Id0+(Xqpp-Xdpp).*Id0.*Iq0;
VR0      = (KE+Ax.*exp(Bx.*Efd0)).*Efd0;
RF0      =  (KF./TF).*Efd0;
Vref     =  V0(1:m)+VR0./KA;
PSV0     =  TM0;
PC       =  PSV0;
%
x0(1:m)          = Eqp0;
x0(m+1:2*m)      = Si1d0;
x0(2*m+1:3*m)    = Edp0;
x0(3*m+1:4*m)    = Si2q0;
x0(4*m+1:5*m)    = D0;
x0(5*m+1:6*m)    = ws_vector;
x0(6*m+1:7*m)    = Efd0;
x0(7*m+1:8*m)    = RF0;
x0(8*m+1:9*m)    = VR0;
x0(9*m+1:10*m)   = TM0;
x0(10*m+1:11*m)  = PSV0;
x011=x0;                      
a0 = [Id0 Iq0 V0 TH0];     
x01=[x011 a0];                 
%% Solving Equation
[bb1, aa1] = ndgrid(1:m, 1:n);                           % Create matrix grid to eliminate for-loops (generator buses)
[bb2, aa2] = ndgrid(m+1:n, 1:n);                         % Create matrix grid to eliminate for-loops (load buses)
M1 = ones(1,11*m); M2 = zeros(1,2*(m+n)); MT = [M1 M2];  % Mass Matrix
M0 = diag(MT);
%% fault or any disturbance 
f=8;
alpha=2; beta=1;
tf1=1; tf2=1.1;

st1 =1e-1;    tsim1 = 1;                                                              % st = step size, tf1 = fault start, tf2 = fault end
tspan1 = [0 , tsim1]; opt = odeset('Mass',M0,'MaxStep',st1);
[t1, x1]=ode15s(@DAEs1,tspan1,x01,opt);  x02=x1(end,:); x02(13*m+f)=1e-3;             %x02(13*m+n+f)=0.0001; 
st2 =1e-4;    tsim2 = tf2-tf1 ;                                                       % st = step size, tf1 = fault start, tf2 = fault end
tspan2 = [0 , tsim2]; opt = odeset('Mass',M0,'MaxStep',st2);
[t2, x2]=ode15s(@DAEs2,tspan2,x02,opt);  x03=[x2(end,1:11*m)  x1(end,11*m+1:end)];    %x02(13*m+f)=0.001;  x02(13*m+n+f)=0.001;
st3 =1e-2;    tsim3 = 10-tf2 ;                                                       % st = step size, tf1 = fault start, tf2 = fault end
tspan3 = [0 , tsim3]; opt = odeset('Mass',M0,'MaxStep',st3);
[t3, x3]=ode15s(@DAEs3,tspan3,x03,opt);  

x=[x1; x2; x3;];
t=[t1; t2 + tf1; t3 + tf2;];
toc
%% ploting
Eqp   =  x(:,1:m);
Si1d  =  x(:,m+1:2*m);
Edp   =  x(:,2*m+1:3*m);
Si2q  =  x(:,3*m+1:4*m);
Delta =  x(:,4*m+1:5*m);
w     =  x(:,5*m+1:6*m);
Efd   =  x(:,6*m+1:7*m);
RF    =  x(:,7*m+1:8*m);
VR    =  x(:,8*m+1:9*m);
TM    =  x(:,9*m+1:10*m);
PSV   =  x(:,10*m+1:11*m);
Id    =  x(:,11*m+1:12*m);
Iq    =  x(:,12*m+1:13*m);
V     =  x(:,13*m+1:13*m+n);
TH    = x(:,13*m+n+1:13*m+2*n);
figure(2); plot(t,Id, 'linewidth',2);grid on;set(gca,'fontsize',16); legend({'G1','G2','G3'}, 'FontSize',12);xlabel('Time (sec.)','FontSize',16,'FontWeight','bold');ylabel('Generator Id (pu)','FontSize',16,'FontWeight','bold');
figure(3); plot(t,Iq, 'linewidth',2);grid on;set(gca,'fontsize',16); legend({'G1','G2','G3'}, 'FontSize',12);xlabel('Time (sec.)','FontSize',16,'FontWeight','bold');ylabel('Generator Iq (pu)','FontSize',16,'FontWeight','bold');
figure(4); plot(t,TH, 'linewidth',2);grid on;set(gca,'fontsize',16); legend({'Bus 1','Bus 2','Bus 3','Bus 4','Bus 5','Bus 6','Bus 7','Bus 8','Bus 9'}, 'FontSize',12);xlabel('Time (sec.)','FontSize',16,'FontWeight','bold');ylabel(' Voltage angle (degree)','FontSize',16,'FontWeight','bold');
figure(5); plot(t,V, 'linewidth',2);grid on;set(gca,'fontsize',16); legend({'Bus 1','Bus 2','Bus 3','Bus 4','Bus 5','Bus 6','Bus 7','Bus 8','Bus 9'}, 'FontSize',12);xlabel('Time (sec.)','FontSize',16,'FontWeight','bold');ylabel(' Voltage magnitude (pu)','FontSize',16,'FontWeight','bold');
figure(6); plot(t,Delta*180/pi, 'linewidth',2);grid on;set(gca,'fontsize',16); legend({'G1','G2','G3'}, 'FontSize',12);xlabel('Time (sec.)','FontSize',16,'FontWeight','bold');ylabel(' Rotor angle (degree)','FontSize',16,'FontWeight','bold');
figure(7); plot(t,w/2/pi, 'linewidth',2);grid on;set(gca,'fontsize',16); legend({'G1','G2','G3'}, 'FontSize',12);xlabel('Time (sec.)','FontSize',16,'FontWeight','bold');ylabel(' Frequency (Hz.)','FontSize',16,'FontWeight','bold');

%END