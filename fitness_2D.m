function error=fitness_2D(laser_real,laser_estimate,version_fitness,NUM_MEASUREMENTS,alpha)
%--------------------------------------------------------------------------
%   Function: fitness_3d
%   Author: Fernando Martin Monar.
%   Date: November, 2015
%--------------------------------------------------------------------------
% -> Description: fitness function that is optimized by the DE-based Global
% localization filter. Laser measurements from the true location are
% compared to laser measurements from an estimate to compute a cost value
% for an estimate. This cost value will be minimized to obtain the solution
% of the GL problem.
%   Polar matching between laser orientations applied. Closest points could
%   be also computed (implemented at the end of the script).
%--------------------------------------------------------------------------
% -> Inputs:
%       -laser_real: Vector with NUM_MEASUREMENTS elements containing the
%       distances of the laser measurements from the true location.
%       -laser_estimate: Vector with NUM_MEASUREMENTS elements containing
%       the estimated distances of the laser measurements from a candidate
%       solution.
%       -version_fitness: Type of fitness function, chosen via keyboard.
%       -err_dis: Sensor noise, standard deviation, in percentage over the
%       distance weighted.
%       -NUM_MEASUREMENTS: Number of horizontal measurements in a scan
% -> Output: 
%       -error: fitness value.
%--------------------------------------------------------------------------
%  - Fitness Function Options (given by version_fitness):
%    1) L2 Norm. Sum of the squared errors (Default).
%    2) L1 Norm. Sum of the absolute values of the error (Mahalanobis).
%    3) Kullback-Leibler Divergence based.
%    4) Density Power Divergence based.
%    5) Hellinger Distance based.
%    6) L2 Norm from Probability Distributions
%    7) L(Variable Exponent) Norm from Probability Distributions
%    8) Generalized Kullback-Leibler Divergence based.
%    9) Itakura-Saito Divergence based.
%    10) Jensen-Shannon Divergence
%    11) Jeffreys Divergence
%--------------------------------------------------------------------------
% -> See also: Global_Local_2D alg_genet_2D
%--------------------------------------------------------------------------

% L2 Norm. (Sum of the squared errors).
if (version_fitness==1)
    % First option to compute it, Sum of the Squared errors.
    error=sum((laser_real - laser_estimate).^2); 
    
    % Second option to compute it, each term is weighted by its variance.
    % error=sum(sum(((laser_real - laser_estimate).^2)./(1+2*(err_dis*laser_real).^2)));
    
    % Third option to compute it, divided by twice the variance to check
    % convergence according to the proposed criteria.
    % difference=laser_real - laser_estimate;
    % std_diff=std(difference);  
    % error=sum((laser_real - laser_estimate).^2)./(2*std_diff^2);
%     L2_T=0;
%     OCL_T=0;
%     for i=1:NUM_MEASUREMENTS
%       
%         L2=0;
%         if floor(laser_real(i)+8) < floor(laser_estimate(i)),
%             %      posible oclusion
%             %      medida       ______^_____
%             %      estimacion   ___________^
%             L2=L2+ floor(laser_real(i))*(0.1^2-0.05^2)  ; % antes 0.1
%             L2=L2+ (0.9^2-0.05^2)  ;
%             L2=L2+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(0.15^2-0.05^2)  ; %antes log(0.15/0.05)*0.15
%             L2=L2+ (0.15^2-0.95^2)  ; 
%         elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
%             %      medida       _________^--
%             %      estimacion   ___________^
%             L2=L2+ floor(laser_real(i))*(0.1^2-0.05^2) ; % antes 0.1
%             L2=L2+ (0.9^2-0.05^2) ;
%             L2=L2+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(0.5^2-0.05^2)  ;
%             L2=L2+ (0.5^2-0.95^2) ;
%         elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
%             % imposible
%             L2=L2+floor(laser_real(i))*(0.95^2-0.05^2);        
%         elseif floor(laser_real(i))>floor(laser_estimate(i)),    
%             %      medida       ___________^
%             %      estimacion   ________^---
%             L2=L2+  floor(laser_estimate(i))*(0.1^2-0.05^2)  ;      
%             L2=L2+  (0.9^2-0.05^2)   ;
%             L2=L2+  (floor(laser_real(i))-ceil(laser_estimate(i)))*(0.99^2-0.05^2) ; 
%             L2=L2+  (0.5^2-0.95^2);
%         elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
%             %      medida       ___________^
%             %      estimacion   ___________^
%             L2=L2+  floor(laser_real(i))*(0.1^2-0.05^2) ;
%             L2=L2+  (0.9^2-0.95^2) ;   
%         end
%         
%         if (laser_real(i)+8 <= laser_estimate(i)),
%             % oclusion
%             OCL_T=OCL_T+1;
%         end
%         L2_T=L2_T+L2;
%     end
%     
%     error=L2_T;
%     error=L2_T*(exp(OCL_T/NUM_MEASUREMENTS));
end

% L1 Norm. (Sum of the absolute values of the error). 
if (version_fitness==2) 
    % First option to compute it, Sum of the absolute values of the error.
    error=sum(abs(laser_real - laser_estimate));
    
    % Second option to compute it, each term is weighted by its std dev.
    % error=sum(sum(abs(laser_real - laser_estimate)./(1+err_dis*laser_real)));
    
    % Third option to compute it, divided by standard deviation to check
    % convergence according to the proposed criteria.
    % difference=laser_real - laser_estimate;
    % std_diff=std(difference);  
    % error=sum(abs(laser_real - laser_estimate))/std_diff;

end

% KL Divergence
if (version_fitness==3)
    
    KLD_T=0;
    OCL_T=0;
    for i=1:NUM_MEASUREMENTS
      
        KLD=0;
        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            KLD=KLD+ floor(laser_real(i))*log(0.1/0.05)*0.1  ; % antes 0.1
            KLD=KLD+ log(0.9/0.05)*0.9  ;
            KLD=KLD+ (floor(laser_estimate(i))-ceil(laser_real(i)))*log(0.15/0.05)*0.15  ; %antes log(0.15/0.05)*0.15
            KLD=KLD+ log(0.15/0.95)*0.15  ; 
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            KLD=KLD+ floor(laser_real(i))*log(0.1/0.05)*0.1 ; % antes 0.1
            KLD=KLD+ log(0.9/0.05)*0.9 ;
            KLD=KLD+ (floor(laser_estimate(i))-ceil(laser_real(i)))*log(0.5/0.05)*0.5  ;
            KLD=KLD+ log(0.5/0.95)*0.5 ;
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            KLD=KLD+floor(laser_real(i))*log(0.95/0.05)*0.95 ;        
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            KLD=KLD+  floor(laser_estimate(i))*log(0.1/0.05)*0.1  ;      
            KLD=KLD+  log(0.9/0.05)*0.9   ;
            KLD=KLD+  (floor(laser_real(i))-ceil(laser_estimate(i)))*log(0.99/0.05)*0.99 ;  
            KLD=KLD+  log(0.5/0.95)*0.5 ;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^
            KLD=KLD+  floor(laser_real(i))*log(0.1/0.05)*0.1 ;
            KLD=KLD+  log(0.9/0.95)*0.9 ;   
        end
        
        if (laser_real(i)+8 <= laser_estimate(i)),
            % oclusion
            OCL_T=OCL_T+1;
        end
        KLD_T=KLD_T+KLD;
    end
    
    error=KLD_T;
    error=KLD_T*(exp(OCL_T/NUM_MEASUREMENTS));
end

% Density Power Divergence
if (version_fitness==4)
    
    DPD=0;
    alpha=1;
    OCL_T=0;
    for i=1:NUM_MEASUREMENTS

        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            DPD=DPD+floor(laser_real(i))*((0.05^alpha)-(1+1/alpha)*0.1*0.05^alpha+(1/alpha)*0.1^(1+alpha));
            DPD=DPD+((0.05^alpha)-(1+1/alpha)*0.9*0.05^alpha+(1/alpha)*0.9^(1+alpha));
            DPD=DPD+(floor(laser_estimate(i))-ceil(laser_real(i)))*((0.05^alpha)-(1+1/alpha)*0.15*0.05^alpha+(1/alpha)*0.15^(1+alpha));
            DPD=DPD+((0.95^alpha)-(1+1/alpha)*0.15*0.95^alpha+(1/alpha)*0.15^(1+alpha));
            OCL_T=OCL_T+1;
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            DPD=DPD+floor(laser_real(i))*((0.05^alpha)-(1+1/alpha)*0.1*0.05^alpha+(1/alpha)*0.1^(1+alpha));
            DPD=DPD+((0.05^alpha)-(1+1/alpha)*0.9*0.05^alpha+(1/alpha)*0.9^(1+alpha));
            DPD=DPD+(floor(laser_estimate(i))-ceil(laser_real(i)))*((0.05^alpha)-(1+1/alpha)*0.5*0.05^alpha+(1/alpha)*0.5^(1+alpha));
            DPD=DPD+((0.95^alpha)-(1+1/alpha)*0.5*0.95^alpha+(1/alpha)*0.5^(1+alpha));
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible  
            DPD=DPD+floor(laser_real(i))*((0.95^alpha)-(1+1/alpha)*0.05*0.95^alpha+(1/alpha)*0.05^(1+alpha));
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            DPD=DPD+floor(laser_estimate(i))*((0.05^alpha)-(1+1/alpha)*0.1*0.05^alpha+(1/alpha)*0.1^(1+alpha));
            DPD=DPD+((0.95^alpha)-(1+1/alpha)*0.1*0.95^alpha+(1/alpha)*0.1^(1+alpha));
            DPD=DPD+(floor(laser_real(i))-ceil(laser_estimate(i)))*((0.5^alpha)-(1+1/alpha)*0.1*0.5^alpha+(1/alpha)*0.1^(1+alpha));
            DPD=DPD+((0.5^alpha)-(1+1/alpha)*0.9*0.5^alpha+(1/alpha)*0.9^(1+alpha));
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^
            DPD=DPD+floor(laser_real(i))*((0.05^alpha)-(1+1/alpha)*0.1*0.05^alpha+(1/alpha)*0.1^(1+alpha));
            DPD=DPD+((0.95^alpha)-(1+1/alpha)*0.9*0.95^alpha+(1/alpha)*0.9^(1+alpha)); 
        end
    end
    error=DPD;
    error=DPD*(exp(OCL_T/NUM_MEASUREMENTS));
end

% Hellinger Distance
if (version_fitness==5)
    
    HD_Total=0;
    %OCL_T=0;
    for i=1:NUM_MEASUREMENTS
        
        HD=0;
        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            HD=HD+floor(laser_real(i))*(sqrt(0.1)-sqrt(0.05))^2;
            HD=HD+(sqrt(0.9)-sqrt(0.05))^2;
            HD=HD+(floor(laser_estimate(i))-ceil(laser_real(i)))*(sqrt(0.15)-sqrt(0.05))^2;
            HD=HD+(sqrt(0.15)-sqrt(0.95))^2;
            %OCL_T=OCL_T+1;
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            HD=HD+floor(laser_real(i))*(sqrt(0.1)-sqrt(0.05))^2;
            HD=HD+(sqrt(0.9)-sqrt(0.05))^2;
            HD=HD+(floor(laser_estimate(i))-ceil(laser_real(i)))*(sqrt(0.5)-sqrt(0.05))^2;
            HD=HD+(sqrt(0.5)-sqrt(0.95))^2;

        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            HD=HD+floor(laser_real(i))*(sqrt(0.95)-sqrt(0.05))^2;   
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            HD=HD+floor(laser_estimate(i))*(sqrt(0.1)-sqrt(0.05))^2;
            HD=HD+(sqrt(0.9)-sqrt(0.05))^2;
            HD=HD+(floor(laser_real(i))-ceil(laser_estimate(i)))*(sqrt(0.99)-sqrt(0.05))^2;
            HD=HD+(sqrt(0.5)-sqrt(0.95))^2;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^            
            HD=HD+floor(laser_real(i))*(sqrt(0.1)-sqrt(0.05))^2;
            HD=HD+(sqrt(0.9)-sqrt(0.95))^2;  
        end
        
        HD_Total=HD_Total+sqrt(HD);
    end
    error=1/sqrt(2)*HD_Total;
    %error=HD_Total*(exp(OCL_T/NUM_MEASUREMENTS));
end

% L2 Norm from Probability Distributions
if (version_fitness==6)
    
    L2_Total=0;
    %OCL_T=0;
    for i=1:NUM_MEASUREMENTS
        
        L2=0;
        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            L2=L2+floor(laser_real(i))*(0.1-0.05)^2;
            L2=L2+(0.9-0.05)^2;
            L2=L2+(floor(laser_estimate(i))-ceil(laser_real(i)))*(0.15-0.05)^2;
            L2=L2+(0.15-0.95)^2;
            %OCL_T=OCL_T+1;
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            L2=L2+floor(laser_real(i))*(0.1-0.05)^2;
            L2=L2+(0.9-0.05)^2;
            L2=L2+(floor(laser_estimate(i))-ceil(laser_real(i)))*(0.5-0.05)^2;
            L2=L2+(0.5-0.95)^2;
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            L2=L2+floor(laser_real(i))*(0.95-0.05)^2;   
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            L2=L2+floor(laser_estimate(i))*(0.1-0.05)^2;
            L2=L2+(0.9-0.05)^2;
            L2=L2+(floor(laser_real(i))-ceil(laser_estimate(i)))*(0.99-0.05)^2;
            L2=L2+(0.5-0.95)^2;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^            
            L2=L2+floor(laser_real(i))*(0.1-0.05)^2;
            L2=L2+(0.9-0.95)^2;  
        end
        
        L2_Total=L2_Total+sqrt(L2);
    end
    error=L2_Total;
    %error=L2_Total*(exp(OCL_T/NUM_MEASUREMENTS));
end

% L(Variable exponent) Norm from Probability Distributions
if (version_fitness==7)
    
    LV_Total=0;
    %OCL_T=0;
    for i=1:NUM_MEASUREMENTS
        
        LV=0;
        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            LV=LV+floor(laser_real(i))*(0.1-0.05)^alpha;
            LV=LV+(0.9-0.05)^alpha;
            LV=LV+(floor(laser_estimate(i))-ceil(laser_real(i)))*(0.15-0.05)^alpha;
            LV=LV+(0.15-0.95)^alpha;
            %OCL_T=OCL_T+1;
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            LV=LV+floor(laser_real(i))*(0.1-0.05)^alpha;
            LV=LV+(0.9-0.05)^alpha;
            LV=LV+(floor(laser_estimate(i))-ceil(laser_real(i)))*(0.5-0.05)^alpha;
            LV=LV+(0.5-0.95)^alpha;
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            LV=LV+floor(laser_real(i))*(0.95-0.05)^alpha;   
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            LV=LV+floor(laser_estimate(i))*(0.1-0.05)^alpha;
            LV=LV+(0.9-0.05)^alpha;
            LV=LV+(floor(laser_real(i))-ceil(laser_estimate(i)))*(0.99-0.05)^alpha;
            LV=LV+(0.5-0.95)^alpha;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^            
            LV=LV+floor(laser_real(i))*(0.1-0.05)^alpha;
            LV=LV+(0.9-0.95)^alpha;  
        end
        
        LV_Total=LV_Total+LV^(1/alpha);
    end
    error=LV_Total;
    %error=L2_Total*(exp(OCL_T/NUM_MEASUREMENTS));
end

% Generalized KL Divergence
if (version_fitness==8)
    
    KLD_T=0;
    OCL_T=0;
    for i=1:NUM_MEASUREMENTS
      
        KLD=0;
        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            KLD=KLD+ floor(laser_real(i))*(log(0.1/0.05)*0.1-0.1+0.05); % antes 0.1
            KLD=KLD+ log(0.9/0.05)*0.9-0.9+0.05  ;
            KLD=KLD+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(log(0.15/0.05)*0.15-0.15+0.05)  ; %antes log(0.15/0.05)*0.15
            KLD=KLD+ log(0.15/0.95)*0.15-0.15+0.95  ; 
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            KLD=KLD+ floor(laser_real(i))*(log(0.1/0.05)*0.1-0.1+0.05) ; % antes 0.1
            KLD=KLD+ log(0.9/0.05)*0.9-0.9+0.05 ;
            KLD=KLD+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(log(0.5/0.05)*0.5-0.5+0.05)  ;
            KLD=KLD+ log(0.5/0.95)*0.5-0.5+0.95 ;
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            KLD=KLD+floor(laser_real(i))*(log(0.95/0.05)*0.95-0.95+0.05);        
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            KLD=KLD+  floor(laser_estimate(i))*(log(0.1/0.05)*0.1-0.1+0.05)  ;      
            KLD=KLD+  log(0.9/0.05)*0.9-0.9+0-05   ;
            KLD=KLD+  (floor(laser_real(i))-ceil(laser_estimate(i)))*(log(0.99/0.05)*0.99-0.99+0.05) ;  
            KLD=KLD+  log(0.5/0.95)*0.5-0.5+0.95 ;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^
            KLD=KLD+  floor(laser_real(i))*(log(0.1/0.05)*0.1-0.1+0.05);
            KLD=KLD+  log(0.9/0.95)*0.9-0.9+0.95 ;   
        end
        
        if (laser_real(i)+8 <= laser_estimate(i)),
            % oclusion
            OCL_T=OCL_T+1;
        end
        KLD_T=KLD_T+KLD;
    end
    
    error=KLD_T*(exp(OCL_T/NUM_MEASUREMENTS));

end

% Itakura-Saito
if (version_fitness==9)
    
    IS=0;
    OCL_T=0;
    for i=1:NUM_MEASUREMENTS

        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^          
            IS=IS+ floor(laser_real(i))*(0.1/0.05-log(0.1/0.05)-1); % antes 0.1
            IS=IS+ 0.9/0.05-log(0.9/0.05)-1;
            IS=IS+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(0.15/0.05-log(0.15/0.05)-1)  ; %antes log(0.15/0.05)*0.15
            IS=IS+ 0.15/0.95-log(0.15/0.95)-1  ; 
            OCL_T=OCL_T+1;
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            IS=IS+ floor(laser_real(i))*(0.1/0.05-log(0.1/0.05)-1) ; % antes 0.1
            IS=IS+ 0.9/0.05-log(0.9/0.05)-1 ;
            IS=IS+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(0.5/0.05-log(0.5/0.05)-1)  ;
            IS=IS+ 0.5/0.95-log(0.5/0.95)-1 ;
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            IS=IS+floor(laser_real(i))*(0.95/0.05-log(0.95/0.05)-1);        
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            IS=IS+  floor(laser_estimate(i))*(0.1/0.05-log(0.1/0.05)-1)  ;      
            IS=IS+  0.9/0.05-log(0.9/0.05)-1  ;
            IS=IS+  (floor(laser_real(i))-ceil(laser_estimate(i)))*(0.99/0.05-log(0.99/0.05)-1) ;  
            IS=IS+  0.5/0.95-log(0.5/0.95)-1 ;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^
            IS=IS+  floor(laser_real(i))*(0.1/0.05-log(0.1/0.05)-1);
            IS=IS+  0.9/0.95-log(0.9/0.95)-1;   
        end
    end
    error=IS;
    error=IS*(exp(OCL_T/NUM_MEASUREMENTS));
end

% Jensen-Shannon Divergence
if (version_fitness==10)
    
    JS=0;
    OCL_T=0;
    for i=1:NUM_MEASUREMENTS
      
        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            JS=JS+ 0.5*floor(laser_real(i))*log(0.1/((0.05+0.1)/2))*0.1 + 0.5*floor(laser_real(i))*log(0.05/((0.05+0.1)/2))*0.05  ; % antes 0.1
            JS=JS+ 0.5*log(0.9/((0.05+0.9)/2))*0.9+0.5*log(0.05/((0.05+0.9)/2))*0.05   ;
            JS=JS+ 0.5*(floor(laser_estimate(i))-ceil(laser_real(i)))*log(0.15/((0.05+0.15)/2))*0.15 +0.5*(floor(laser_estimate(i))-ceil(laser_real(i)))*log(0.05/((0.05+0.15)/2))*0.05  ; %antes log(0.15/0.05)*0.15
            JS=JS+ 0.5*log(0.15/((0.95+0.15)/2))*0.15+ 0.5*log(0.95/((0.95+0.15)/2))*0.95  ; 
            OCL_T=OCL_T+1;
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            JS=JS+ 0.5* floor(laser_real(i))*log(0.1/((0.05+0.1)/2))*0.1 +0.5* floor(laser_real(i))*log(0.05/((0.05+0.1)/2))*0.05 ; % antes 0.1
            JS=JS+ 0.5*log(0.9/((0.05+0.9)/2))*0.9 +0.5*log(0.05/((0.05+0.9)/2))*0.05 ;
            JS=JS+ 0.5*(floor(laser_estimate(i))-ceil(laser_real(i)))*log(0.5/((0.05+0.5)/2))*0.5+ 0.5*(floor(laser_estimate(i))-ceil(laser_real(i)))*log(0.05/((0.05+0.5)/2))*0.05  ;
            JS=JS+ 0.5*log(0.5/((0.95+0.5)/2))*0.5 +0.5*log(0.95/((0.95+0.5)/2))*0.95;
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            JS=JS+0.5*floor(laser_real(i))*log(0.95/((0.05+0.95)/2))*0.95 + +0.5*floor(laser_real(i))*log(0.05/((0.05+0.95)/2))*0.05  ;        
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            JS=JS+  0.5*floor(laser_estimate(i))*log(0.1/((0.05+0.1)/2))*0.1 + 0.5*floor(laser_estimate(i))*log(0.05/((0.05+0.1)/2))*0.05  ;      
            JS=JS+  0.5* log(0.9/((0.05+0.9)/2))*0.9 + 0.5* log(0.05/((0.05+0.9)/2))*0.05   ;
            JS=JS+  0.5*(floor(laser_real(i))-ceil(laser_estimate(i)))*log(0.99/((0.05+0.99)/2))*0.99+ 0.5*(floor(laser_real(i))-ceil(laser_estimate(i)))*log(0.05/((0.05+0.99)/2))*0.05 ;  
            JS=JS+  0.5*log(0.5/((0.5+0.95)/2))*0.5 +0.5*log(0.95/((0.95+0.5)/2))*0.95 ;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^
            JS=JS+  0.5*floor(laser_real(i))*log(0.1/((0.05+0.1)/2))*0.1+ 0.5*floor(laser_real(i))*log(0.05/((0.05+0.1)/2))*0.05  ;
            JS=JS+  0.5*log(0.9/((0.95+0.9)/2))*0.9+ 0.5*log(0.95/((0.95+0.9)/2))*0.95 ;   
        end
       
    end
    error=JS;
    error=JS*(exp(OCL_T/NUM_MEASUREMENTS));

end


% Jeffreys
if (version_fitness==11)
    
    JF_T=0;
    OCL_T=0;
    for i=1:NUM_MEASUREMENTS
      
        JF=0;
        if floor(laser_real(i)+8) < floor(laser_estimate(i)),
            %      posible oclusion
            %      medida       ______^_____
            %      estimacion   ___________^
            JF=JF+ floor(laser_real(i))*(0.1-0.05)*(log(0.1)-log(0.05)); % antes 0.1
            JF=JF+ (0.9-0.05)*(log(0.9)-log(0.05))  ;
            JF=JF+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(0.15-0.05)*(log(0.15)-log(0.05))  ; %antes log(0.15/0.05)*0.15
            JF=JF+ (0.15-0.95)*(log(0.15)-log(0.95))  ; 
        elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
            %      medida       _________^--
            %      estimacion   ___________^
            JF=JF+ floor(laser_real(i))*(0.1-0.05)*(log(0.1)-log(0.05)); % antes 0.1
            JF=JF+ (0.9-0.05)*(log(0.9)-log(0.05)) ;
            JF=JF+ (floor(laser_estimate(i))-ceil(laser_real(i)))*(0.5-0.05)*(log(0.5)-log(0.05))  ;
            JF=JF+ (0.5-0.95)*(log(0.5)-log(0.95)) ;
        elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
            % imposible
            JF=JF+floor(laser_real(i))*(0.95-0.05)*(log(0.95)-log(0.05)) ;        
        elseif floor(laser_real(i))>floor(laser_estimate(i)),    
            %      medida       ___________^
            %      estimacion   ________^---
            JF=JF+  floor(laser_estimate(i))*(0.1-0.05)*(log(0.1)-log(0.05))  ;      
            JF=JF+  (0.9-0.05)*(log(0.9)-log(0.05))   ;
            JF=JF+  (floor(laser_real(i))-ceil(laser_estimate(i)))*(0.99-0.05)*(log(0.99)-log(0.05)) ;  
            JF=JF+  (0.5-0.95)*(log(0.5)-log(0.95)) ;
        elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
            %      medida       ___________^
            %      estimacion   ___________^
            JF=JF+  floor(laser_real(i))*(0.1-0.05)*(log(0.1)-log(0.05)) ;
            JF=JF+  (0.9-0.95)*(log(0.9)-log(0.95)) ;   
        end
        
        if (laser_real(i)+8 <= laser_estimate(i)),
            % oclusion
            OCL_T=OCL_T+1;
        end
        JF_T=JF_T+JF;
    end
    
    error=JF_T;
    error=JF_T*(exp(OCL_T/NUM_MEASUREMENTS));
end


% Havrda-Charvat
% if (version_fitness==5)
% 
%     HC=0;    
%     for i=1:NUM_MEASUREMENTS
%         if floor(laser_real(i)+8) < floor(laser_estimate(i)),
%             %      posible oclusion
%             %      medida       ______^_____
%             %      estimacion   ___________^
%            HC=HC+floor(laser_real(i))*(1/alpha/(1-alpha))*((0.1^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(1/alpha/(1-alpha))*((0.9^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(floor(laser_estimate(i))-ceil(laser_real(i)))*(1/alpha/(1-alpha))*((0.15^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(1/alpha/(1-alpha))*((0.15^alpha)*(0.95^(1-alpha))-1);
%         elseif floor(laser_real(i)+2)<floor(laser_estimate(i)),
%             %      medida       _________^--
%             %      estimacion   ___________^ 
%            HC=HC+floor(laser_real(i))*(1/alpha/(1-alpha))*((0.1^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(1/alpha/(1-alpha))*((0.9^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(floor(laser_estimate(i))-ceil(laser_real(i)))*(1/alpha/(1-alpha))*((0.5^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(1/alpha/(1-alpha))*((0.5^alpha)*(0.95^(1-alpha))-1);
%         elseif floor(laser_real(i))>floor(laser_estimate(i)+6),  
%             % imposible  
%             HC=HC+floor(laser_real(i))*(1/alpha/(1-alpha))*((0.95^alpha)*(0.05^(1-alpha))-1);
%         elseif floor(laser_real(i))>floor(laser_estimate(i)),    
%             %      medida       ___________^
%             %      estimacion   ________^--- 
%            HC=HC+floor(laser_estimate(i))*(1/alpha/(1-alpha))*((0.1^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(1/alpha/(1-alpha))*((0.1^alpha)*(0.95^(1-alpha))-1);
%            HC=HC+(floor(laser_real(i))-ceil(laser_estimate(i)))*(1/alpha/(1-alpha))*((0.1^alpha)*(0.5^(1-alpha))-1);
%            HC=HC+(1/alpha/(1-alpha))*((0.9^alpha)*(0.5^(1-alpha))-1);
%         elseif floor(laser_real(i))<=floor(laser_estimate(i)),   
%             %      medida       ___________^
%             %      estimacion   ___________^
%            HC=HC+floor(laser_real(i))*(1/alpha/(1-alpha))*((0.1^alpha)*(0.05^(1-alpha))-1);
%            HC=HC+(1/alpha/(1-alpha))*((0.9^alpha)*(0.95^(1-alpha))-1);
%         end
%     end
%     error=abs(HC);
% end



% To use kd trees searching closest points.
%ptrtree=BuildGLTree3D(cart_real_3d');
%[kNNG,dist]=KNNSearch3D(cart_real_3d',cart_est_3d',ptrtree,1);
%error=sqrt(dist'*dist);
%error=sqrt((dist'*dist)/size(cart_real_3d,1));

end


