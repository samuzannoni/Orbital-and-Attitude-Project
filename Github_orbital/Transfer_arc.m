function [DV_constrained,DV_unconstrained,min_V,date_dep,date_arr,t_dep_index,t_arr_index,DV1,DV2,VI_T,RI_T,RF_T] = Transfer_arc(t_span_dep,t_span_arr,MU,DEP_eph,ARR_eph,V_inf)

%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_span_dep and t_span_arr has to be in julian days (or fraction of it)
% Mu: gravitational constant of the main attractor
% DEP_eph: Matrix containing ephemerides of the first planet at departure
% ARR_eph: Matrix containing ephemerides of the second planet at arrival 
% V_inf = Launcher constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DV_constrained: Matrix containing all the DeltaV satisfying the launcher
% constraint
% min_V: Minimum Delta_V 
% date_dep: departure date corresponding to the minimum DeltaV
% date_arr: arrival date corresponding to the minimum DeltaV
% DV1: Matrix 1000x1000 containing initial Delta V
% DV2: Matrix 1000x1000 containing final Delta V
% VI_T: Initial velocity of the Transfer arc
% RI_T: Initial position vector of Transfer arc
% RF_T: Final position vector of Transfer arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for k = 1:length(t_span_dep)


    RI = DEP_eph(k,1:3); % position vector of Earth
    RI = RI'; % transpose because becomes input of LambertMR

        for j = 1:length(t_span_arr)
    
            RF = ARR_eph(j,1:3); % position vector of Mars
            
            RF = RF'; % transpose because becomes input of LambertMR
        
            TOF = (t_span_arr(j) - t_span_dep(k)) *24*3600; % Time of flight converted in seconds
            orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
            Nrev = 0; %Zero-revolution case
            Ncase = 0; %Not used for the zero-revolution case
            optionsLMR = 0; %No display
    
    
            [~,~,~,~,VI_T,VF_T,~,~] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR); % for every RI I get N transfer arc (N length of time vector)
        
            vi_dep = DEP_eph(k,4:6); % velocity vector of Earth at every instant of time
            vf_arr = ARR_eph(j,4:6); % velocity vector of Mars at every instant of time
        
        
            Delta_V1 = VI_T - vi_dep;
            Delta_V2 = vf_arr - VF_T;
        
            Delta_V = norm(Delta_V1) + norm(Delta_V2);

            DV1(k,j) = norm(Delta_V1);
            DV2(k,j) = norm(Delta_V2);

            DV_unconstrained(k,j) = Delta_V; % Matrix of Delta Vs without satisfying the Launcher constraint
            

        end

           
end


min_DV1 = min(min(DV1)); % Minimum Delta V1




% Here is evaluated the constraint

if min_DV1 <= V_inf

     [t_dep_index, ~] = find(DV1 == min_DV1);   % index where the minimum Delta V1 fulfill the constraint


    
     Delta_V1 = min_DV1;
     VI_T = Delta_V1 - DEP_eph(t_dep_index,4:6); % Initial velocity of the transfer

  % computing all the transfer arcs with position departure fixed


     RI_T = DEP_eph(t_dep_index,1:3); % Initial position vector of transfer arc 
     RI_T = RI_T'; % transpose because becomes input of LambertMR

     for j = 1:length(t_span_arr)
    
            RF_T = ARR_eph(j,1:3); % position vector of Mars
            
            RF_T = RF_T'; % transpose because becomes input of LambertMR
        
            TOF = (t_span_arr(j) - t_span_dep(t_dep_index)) *24*3600; % Time of flight converted in seconds
            orbitType = 0; %Direct orbit (0: direct; 1: retrograde)
            Nrev = 0; %Zero-revolution case
            Ncase = 0; %Not used for the zero-revolution case
            optionsLMR = 0; %No display

   
    
           [a_T,~,~,~,VI_T,VF_T,~,~] = lambertMR(RI_T,RF_T,TOF,MU,orbitType,Nrev,Ncase,optionsLMR);



            vf_arr = ARR_eph(j,4:6); % velocity vector of the planet at every instant of time
            Delta_V2 = vf_arr - VF_T;
        
            Delta_V = Delta_V1 + norm(Delta_V2); % now Delta_V1 is always the same
        
            DV_constrained(:,j) = Delta_V;  % Matrix of Delta Vs satisfying the Launcher constraint (should be a matrix)


     end



     %min_V = min(min(DV_constrained));
     min_V = min(DV_constrained); % minimum Delta V satisfying the launcher constraint

     [~,t_arr_index] = find(DV_constrained==min_V);
     RF_T = ARR_eph(t_arr_index,1:3); % Final position vector of transfer arc


     date_dep = mjd20002date(t_span_dep(t_dep_index));
     date_arr = mjd20002date (t_span_arr(t_arr_index));


     disp('The cheapest transfer (satisfying the launcher constraint) has been computed:');
     fprintf('\n');
     fprintf('The departure date is [yyyy, mm, dd, UTC] %d %d %d %d:%d:%.2g \n', date_dep);
     fprintf('The arrival date is [yyyy, mm, dd, UTC] %d %d %d %d:%d:%.2g \n', date_arr);
     fprintf('The associated initial minimum cost is %.5g [Km/s] \n' ,min_DV1);
     fprintf('The associated total minimum cost is %.5g [Km/s] \n' ,min_V);




else

    disp('The transfer can''t be computed because of Launcher constraint')

end





