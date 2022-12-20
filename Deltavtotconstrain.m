function [DV_tot,DV1,min_DV_tot,dv1_verificata,date1,date2,ROW,COL] = Deltavtotconstrain(iplanet1,iplanet2,t_span_dep,t_span_arr,mu,n1,n2,V_inf)

%%% INPUT:
% t_span_dep: departure window (vector)
% t_span_arr: arrival window ( vector)
% iplanet: planet markers of uplanet.m
% mu is related to the celestial body around which the planets rotate
% n:step time span

%%% OUPUT:
% DV_tot: Matrix containing DeltaV
% min_DV_tot: minimum DeltaV, minimum cost of the manouvre
% date1: date of departure
% date2: date of arrival
DV1=[]; 
DV_tot=[];
DEP_eph = [];
ARR_eph = [];
A=[]; 
row=[];
col=[];
for k1 = 1:length(t_span_dep)
    [kep_dep_E,~] = uplanet(t_span_dep(k1), iplanet1); % ephemerides of Earth at departure
    [r_dep_E,v_dep_E] = kep2car(kep_dep_E(1),kep_dep_E(2),kep_dep_E(3),kep_dep_E(4),kep_dep_E(5),kep_dep_E(6),mu); % converted in cartesian coordinates
    s_dep_E = [r_dep_E' v_dep_E']; % Matrix containing Earth's state for departure window
    DEP_eph = [DEP_eph; s_dep_E]; % matrix containing all Earth's states for departure time window
end
    
for k2 = 1:length(t_span_arr)
    [kep_arr_M,~] = uplanet(t_span_arr(k2), iplanet2); % ephemerides of Mars at arrival window
    [r_arr_M,v_arr_M] = kep2car(kep_arr_M(1),kep_arr_M(2),kep_arr_M(3),kep_arr_M(4),kep_arr_M(5),kep_arr_M(6),mu); % converted in cartesian coordinates
    s_arr_M = [r_arr_M' v_arr_M']; % Matrix containing Mars's state for arrival window
    ARR_eph = [ARR_eph; s_arr_M]; % matrix containing all Mars's states for arrival time window
end
    
for i = 1:n1
    for j = 1:n2
        RI = DEP_eph(i,1:3); % position vector of Earth
        RF = ARR_eph(j,1:3); % position vector of Mars
        vi_dep = DEP_eph(i,4:6); % velocity vector of Earth
        vf_arr = ARR_eph(j,4:6); % velocity vector of Mars 
             
        TOF = (t_span_arr(j)-t_span_dep(i))*24*3600; % Time of flight converted in seconds
        
        orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
        Nrev = 0; %Zero-revolution case
        Ncase = 0; %Not used for the zero-revolution case
        optionsLMR = 0; %No display
        [~,~,~,~,VI,VF,~,~] = lambertMR(RI,RF,TOF,mu,orbitType,Nrev,Ncase,optionsLMR);

        Delta_V1 = VI - vi_dep;
        Delta_V2 = vf_arr - VF;
        
        Delta_Vtot = norm(Delta_V1) + norm(Delta_V2);
        DV1(i,j) = [norm(Delta_V1)];
      
        DV_tot(i,j) = [Delta_Vtot]; % Matrix containing DeltaV
        
    end
end

[row,col]=find(DV1<=V_inf);

DV_totbest=[];
DV1_constrain= [];
for i = 1:length(row)
    for j = 1:length(col)
        A=DV_tot(row(i),col(j));
        B=DV1(row(i),col(j));
  
        DV_totbest(i,j)= [A];
        DV1_constrain(i,j)=[B];
       end
end 

dv_totbest=diag(DV_totbest);
dv1_constrain=diag(DV1_constrain);

min_DV_tot = min(dv_totbest);

[k]=find(dv_totbest==min_DV_tot);

dv1_verificata= dv1_constrain(k); % verifica che sia minore di V_inf


if dv1_verificata<=V_inf
    
    ROW = row(k);
    COL= col(k);
 
   date1 = mjd20002date(t_span_dep(ROW));
   date2 = mjd20002date(t_span_arr(COL));
else
    disp('error: vincolo non verificato')
end

end
