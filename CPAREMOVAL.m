clear;

%%
%%%%%%%%%%% Structure parameters %%%%%%%%%%%%%%%%%
%%%% subscript b: blood/cell suspension
%%%% subscript d: dialysate/replacement fluid
Hc_b=500*10^(-6);               % Hc: channel height, unit: m
Wc_b=600*10^(-6);               % Wc: channel width, unit: m
Hc_d=500*10^(-6);                
Wc_d=600*10^(-6); 

Ac_b=Hc_b*Wc_b;                 % Ac: channel area, unit: m
Ac_d=Hc_d*Wc_d;

L=0.1;                          % L: channel length, unit: m
%%
%%%%%%%%%%% Mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CVNum=2000;
nodeN=CVNum+2;        
edgeN=CVNum+2;
Dx=zeros(1,nodeN);             % x distance between two edges
dx=zeros(1,edgeN);             % x distance between two nodes

for i=1:nodeN    
    Dx(i)=L/CVNum;
    Dx(1)=0;
    Dx(nodeN)=0;
end

for i=1:edgeN    
    if i>1
        dx(i)=(Dx(i)+Dx(i-1))/2;
    else
    	dx(i)=0;
    end
end
%%
%%%%%%%%%%% Physical parameters %%%%%%%%%%%%%%%%%
R=8.314472;                     % R
T=295;                          % Temperature 

den_w=1.0*10^3;                 % density of water
Lp=1.74*10^(-12);               % Memberane permeability, unit: m/(Pa*s)*************************************************
Ks=6.61*10^(-8);                % Transmembrane mass transfer coefficient, unit: m/s *****************************************************
deta=0.841;                     % Solute reflection coefficient ******************************************************************
nDg=0.5*10^(-9);                % Diffusion coefficient of glycerol %%%%%%%

global Vs;
Vs=73.03237*10^(-6);            % partial volume of glycerol;Vs=Ms/(1000*den_g);
Ms=92.09382;                    % Molar mass of glycerol, g/mol
den_g=1.261*10^3;               % density of glycerol
n=1*10^-6;                      % viscous coefficient     ***************
% ks=den_w*Ms/den_g;            % 
%%
%%%%%%%%%%% Boudary/Initial condition %%%%%%%%%%%%%%%%%
cs_b0=400;                        % Initial glyceroal concentration, unit: g/L  
den_b0=cs_b0+(1-cs_b0/den_g)*den_w;
den_b=zeros(1,nodeN)+den_b0;      % 

cs_d0=0;
den_d0=cs_d0+(1-cs_d0/den_g)*den_w;
den_d=zeros(1,nodeN)+den_d0;

Nsall=1000*cs_b0/Ms;             % mole number, unit: mole per m^3
ms_b0=Nsall;                     % Initial osmolarity of solution, unit: mole per m^3.
ms_b(1,:)=zeros(1,nodeN)+ms_b0;
ms_d0=0;
ms_d(1,:)=zeros(1,nodeN)+ms_d0;   % Assumed a priming with glycerol solution

Qb=(20/36)*10^(-8);                  % Inlet flow rate, unit: m^3/s
Qd=(25/36)*10^(-8);                  % 

v_b0=Qb/Ac_b;                     % Inlet velocity, unit: m/s 
v_d0=Qd/Ac_d;

v_b=zeros(1,nodeN)+v_b0;          % 
v_d=zeros(1,nodeN)+v_d0;
dv_b=zeros(1,nodeN);              % unit variation of velocity along x
dv_d=zeros(1,nodeN);

Jw_b(1,:)=zeros(1,nodeN);         % water flux, unit: m^3/(m^2 s), m/s
Js_b(1,:)=zeros(1,nodeN);         % solute flux, unit: kg/(m^2 s)
Jw_d(1,:)=zeros(1,nodeN);
Js_d(1,:)=zeros(1,nodeN); 

P_b=zeros(1,nodeN);
P_d=zeros(1,nodeN);

dms=zeros(1,nodeN);
avg_ms=zeros(1,nodeN);
dP=zeros(1,nodeN);


%%
%%%%%%%%%%%% Numerical Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;
kMax=10000;
t(k)=0;
t_V=0.01;                       % relative volume variation per iteration
% dt=0.001;

while k<=kMax
    %% Update transmembrane mass flux %%
    for i=1:nodeN
        dms(i)=ms_b(k,i)-ms_d(k,nodeN-i+1);
        avg_ms(i)=(ms_b(k,i)+ms_d(k,nodeN-i+1))/2;
        dP(i)=P_b(i)-P_d(nodeN-i+1);
        Jw_b(i)=-Lp*(dP(i)-deta*R*T*dms(i));             % K-K
        Js_b(i)=-((1-deta)*Jw_b(i)*avg_ms(i)+Ks*dms(i));
        
        Jw_d(nodeN-i+1)=-Jw_b(i);
        Js_d(nodeN-i+1)=-Js_b(i);
    end    
    
    %% Update velocity of cell suspension %%%%%%%
    dv_b=(Js_b*Vs+Jw_b)/Hc_b;    
    for j=2:nodeN-1        
        v_b(j)=v_b(j-1)+dv_b(j)*dx(j);
    end
    v_b(nodeN)=v_b(nodeN-1);

    dt_b=t_V/max(abs(dv_b));
    
    %% Update velocity of dialysate %%%%%%%
    dv_d=(Js_d*Vs+Jw_d)/Hc_d;    
    for j=2:nodeN-1        
        v_d(j)=v_d(j-1)+dv_d(j)*dx(j);
    end
    v_d(nodeN)=v_d(nodeN-1);

    dt_d=t_V/max(abs(dv_d));
    
    %% Move to the next time level %%%%%%%
    k=k+1;
    dt=min(0.001,min(dt_b,dt_d));   
    
    %% Update the solute concentration: explict method %%
    ms_b(k,1)=ms_b0;
    ms_d(k,1)=ms_d0;
    
    for i=2:nodeN
        if i<nodeN
            dms_b_1=nDg*dt*((ms_b(k-1,i+1)-ms_b(k-1,i))/dx(i+1)-(ms_b(k-1,i)-ms_b(k-1,i-1))/dx(i))/Dx(i);           %% variation due to diffusion
            dms_d_1=nDg*dt*((ms_d(k-1,i+1)-ms_d(k-1,i))/dx(i+1)-(ms_d(k-1,i)-ms_d(k-1,i-1))/dx(i))/Dx(i);
            
            dms_b_2=(v_b(i-1)*ms_b(k-1,i-1)-v_b(i)*ms_b(k-1,i))*dt/Dx(i);                                           %% variation due to convection
            dms_d_2=(v_d(i-1)*ms_d(k-1,i-1)-v_d(i)*ms_d(k-1,i))*dt/Dx(i);                                           %% first order upwind **
        else
            dms_b_1=0;
            dms_d_1=0;
            
            dms_b_2=0;
            dms_d_2=0;
        end
        
        dms_b_3=dt*(Js_b(i)-(Jw_b(i)+Js_b(i)*Vs)*ms_b(k-1,i))/(Hc_b+(Jw_b(i)+Js_b(i)*Vs)*dt);                      %% variation due to transmembrane flux
        dms_d_3=dt*(Js_d(i)-(Jw_d(i)+Js_d(i)*Vs)*ms_d(k-1,i))/(Hc_d+(Jw_d(i)+Js_d(i)*Vs)*dt);
        
        ms_b(k,i)=ms_b(k-1,i)+dms_b_1+dms_b_2+dms_b_3;
        ms_d(k,i)=ms_d(k-1,i)+dms_d_1+dms_d_2+dms_d_3;
    end

    
    %% Update fluid density %%%%%%%
    for i=1:nodeN        
         den_b(i)=ms_b(k,i)*Ms/1000+(1-ms_b(k,i)*Vs)*den_w;
         den_d(i)=ms_d(k,i)*Ms/1000+(1-ms_d(k,i)*Vs)*den_w;
    end
    
    %% Update pressure feild %%%%%%%
    P_b(nodeN)=0;
    P_d(nodeN)=0;
    for i=nodeN-1:-1:2
        P_b(i)=P_b(i+1)+den_b(i+1)*v_b(i+1)*v_b(i+1)-den_b(i)*v_b(i)*v_b(i)-n*(dv_b(i+1)-dv_b(i));
        P_d(i)=P_d(i+1)+den_d(i+1)*v_d(i+1)*v_d(i+1)-den_d(i)*v_d(i)*v_d(i)-n*(dv_d(i+1)-dv_d(i));
    end    
end
