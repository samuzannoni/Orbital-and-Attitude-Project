clc
clear all
close all 

addpath('C:\Users\sferl\OneDrive - Politecnico di Milano\POLI\SPACE ENGINEERING\Orbital 22_23\MATLAB functions\');

%%%%%%%%%%%%%%%%%%%%%% Departure time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_dep_in = date2mjd2000([2003 04 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_dep_fin = date2mjd2000([2003 08 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_dep = linspace(t_dep_in,t_dep_fin, t_dep_fin - t_dep_in); % time vector for departure window (be careful it's in days!!!)




%%%%%%%%%%%%%%%%%%%%%%% Arrival time window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_arr_in = date2mjd2000([2003 09 01 12 00 00]); %initial departure time given in number of days since 01 01 2000 at 12:00
t_arr_fin = date2mjd2000([2004 03 01 12 00 00]); %final departure time given in number of days since 01 01 2000 at 12:00

t_span_arr = linspace(t_arr_in,t_arr_fin, t_arr_fin - t_arr_in); % time vector for departure window (be careful it's in days!!!)




%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Earth's ephemerides %%%%%%%%%%%%%%%%%%%%%%

DEP_eph_E = [];
for k1 = 1:length(t_span_dep)

[kep_dep_E,ksun] = uplanet(t_span_dep(k1), 3); % ephemerides of Earth at departure

a = kep_dep_E(1);
e = kep_dep_E(2);
i = kep_dep_E(3);
Om =  kep_dep_E(4);
om =  kep_dep_E(5);
f0 =  kep_dep_E(6);

mu = ksun; 
[r_dep_E,v_dep_E] = kep2car_mod(a,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_dep_E = [r_dep_E' v_dep_E']; % Matrix containing Earth's state for departure window

DEP_eph_E = [DEP_eph_E; s_dep_E]; % matrix containing all Earth's states for departure time window

end

% ARR_eph_E = []; 
% for k2 = 1:length(t_span_arr)
% 
% [kep_arr_E,ksun] = uplanet(t_span_arr(k2), 3); % ephemerides of Earth at arrival
% 
% a = kep_arr_E(1);
% e = kep_arr_E(2);
% i = kep_arr_E(3);
% Om =  kep_arr_E(4);
% om =  kep_arr_E(5);
% f0 =  kep_arr_E(6);
% 
% [r_arr_E,v_arr_E] = kep2car_mod(a,e,i,Om,om,f0); % converted in cartesian coordinates
% s_arr_E = [r_arr_E' v_arr_E']; % Matrix containing Earth's state for arrival window
% 
% 
% ARR_eph_E = [ARR_eph_E; s_arr_E]; % matrix containing all Earth's states for arrival time window
% end







%%%%%%%%%%%%%%%%%%%%%%%% Loops to compute Mars's ephemerides %%%%%%%%%%%%%%%%%%%%%%

% DEP_eph_M = [];
% for k1 = 1:length(t_span_dep)
% 
% [kep_dep_M,ksun] = uplanet(t_span_dep(k1), 3); % ephemerides of Earth at departure
% 
% a = kep_dep_M(1);
% e = kep_dep_M(2);
% i = kep_dep_M(3);
% Om =  kep_dep_M(4);
% om =  kep_dep_M(5);
% f0 =  kep_dep_M(6);
% 
% [r_dep_M,v_dep_M] = kep2car(a,e,i,Om,om,f0); % converted in cartesian coordinates
% s_dep_M = [r_dep_M' v_dep_M']; % Matrix containing Mars's state for departure window
% 
% 
% DEP_eph_M = [DEP_eph_M; s_dep_M]; % matrix containing all Mars's states for departure time window
% 
% end

ARR_eph_M = [];
mu = ksun;
for k2 = 1:length(t_span_arr)

[kep_arr_M,ksun] = uplanet(t_span_arr(k2), 3); % ephemerides of Mars at arrival window

a = kep_arr_M(1);
e = kep_arr_M(2);
i = kep_arr_M(3);
Om =  kep_arr_M(4);
om =  kep_arr_M(5);
f0 =  kep_arr_M(6);

[r_arr_M,v_arr_M] = kep2car_mod(a,e,i,Om,om,f0,mu); % converted in cartesian coordinates
s_arr_M = [r_arr_M' v_arr_M']; % Matrix containing Mars's state for arrival window


ARR_eph_M = [ARR_eph_M; s_arr_M]; % matrix containing all Mars's states for arrival time window
end

%%%%%%%%%%%%%%%%%%%%% COMPUTING OF LAMBERT ARCS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Lambert arcs between departure windows %%%%%%%%%%%%%

% DV1 = [];
% MU = ksun;
% 
% for k = 1:length(t_span_dep)
% 
%     TOF = (t_span_arr(k) - t_span_dep(k)) *24*3600; % Time of flight converted in seconds
% 
%     RI = DEP_eph_E(k,1:3); % position vector of Earth
%     RF = DEP_eph_M(k,1:3); % position vector of Mars
% 
%     orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
%     Nrev = 0; %Zero-revolution case
%     Ncase = 0; %Not used for the zero-revolution case
%     optionsLMR = 0; %No display
% 
% 
%     [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);
% 
%     vi_dep = DEP_eph_E(k,4:6); % velocity vector of Earth
%     vf_dep = DEP_eph_M(k,4:6); % velocity vector of Mars 
% 
%    
%     Delta_V1 = VI - vi_dep;
%     Delta_V2 = vf_dep - VF;
% 
%     Delta_V = norm(Delta_V1) + norm(Delta_V2);
% 
%     DV1 = [DV1; Delta_V]; % Matrix containing DeltaV
% end
% 
% 
% %%%%%%%%%% Lambert arcs between departure and arrival windows %%%%%%%%%%%%
    

MU = ksun;
DV2 = [];

for k = 1:length(t_span_dep)

    TOF = (t_span_arr(k) - t_span_dep(k)) *24*3600; % Time of flight converted in seconds

    RI = DEP_eph_E(k,1:3); % position vector of Earth
    RF = ARR_eph_M(k,1:3); % position vector of Mars
    RI = RI';
    RF = RF';

    orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
    Nrev = 0; %Zero-revolution case
    Ncase = 0; %Not used for the zero-revolution case
    optionsLMR = 0; %No display


    [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);

    vi_dep = DEP_eph_E(k,4:6); % velocity vector of Earth
    vf_arr = ARR_eph_M(k,4:6); % velocity vector of Mars 


    Delta_V1 = VI - vi_dep;
    Delta_V2 = vf_arr - VF;

    Delta_V = norm(Delta_V1) + norm(Delta_V2);

    DV2 = [DV2; Delta_V]; % Matrix containing DeltaV
end



%%%%%%%%%%%%%%% Lambert arcs between both arrival windows %%%%%%%%%%%%%%%%
             
% DV3 = [];
% 
% for k = 1:length(t_span_dep)
% 
%     TOF = (t_span_arr(k) - t_span_dep(k)) *24*3600; % Time of flight converted in seconds
% 
%     RI = ARR_eph_E(k,1:3); % position vector of Earth
%     RF = ARR_eph_M(k,1:3); % position vector of Mars
% 
%     orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
%     Nrev = 0; %Zero-revolution case
%     Ncase = 0; %Not used for the zero-revolution case
%     optionsLMR = 0; %No display
% 
% 
%     [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);
% 
%     vi_arr = DEP_eph_E(k,4:6); % velocity vector of Earth
%     vf_arr = ARR_eph_M(k,4:6); % velocity vector of Mars 
% 
% 
%     Delta_V1 = VI - vi_arr;
%     Delta_V2 = vf_arr - VF;
% 
%     Delta_V = norm(Delta_V1) + norm(Delta_V2);
% 
%     DV3 = [DV3; Delta_V]; % Matrix containing DeltaV
% 
% end


%%%%%%% plot the porkchop plot %%%%%%%
