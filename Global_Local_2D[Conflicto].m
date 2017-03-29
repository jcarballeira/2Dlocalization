function [Solution]=Global_Local_2D
%--------------------------------------------------------------------------
%   Main Program: Global_Local_2D
%   Author: Fernando Martin Monar.
%   Date: Nobember, 2015
%--------------------------------------------------------------------------
% -> Description: Global Localization module based on the Differential
%   Evolution filter. The objective is to estimate the robot's location
%   given a laser scan and the known map. Robot's motion is possible after
%   convergence.
%   The robot's true pose is given via keyboard (initialization.m) and the
%   main programs calls alg_genet_3D, which contains the DE-based GL
%   filter. After convergence, the robot's pose estimate is returned in
%   bestmem.
%--------------------------------------------------------------------------
% -> Inputs: no inputs are required.
%
% -> Outputs: 
%       -Solution:
%           -pose_estimate: Vector of D elements containing the solution of
%           the global localization filter (robot's location). The first
%           two coordinates are the cartesian coordinates, in cells, and
%           the third one is the orientation.
%           -error: Cost value of the solution.
%--------------------------------------------------------------------------
% -> Usage: Solution=Global_Local_2D
%--------------------------------------------------------------------------
% Different options for the GL algorithm can be selected via keyboard:
%  - DE Core Options:
%    1) Random Mutation, with Thresholding and Discarding (Default). 
%    2) Basic version, Random Mutation, without Thresholding, Discarding.
%    3) Mutation from Best candidate, with Thresholding and Discarding.
%    4) Random Mutation, with Thresholding and Discarding, NP is
%    drastically reduced (tracking) after convergence.
%  - Fitness Function Options:
%    1) L2 Norm. Sum of the squared errors (Default).
%    2) L1 Norm. Sum of the absolute values of the error (Mahalanobis). 
%    3) Kullback-Leibler Divergence based.
%    4) Density Power Divergence based.
%    5) Hellinger Distance based.
%    6) L2 Norm from Probability Distributions
%    7) L(Variable Exponent) Norm from Probability Distributions
%    8) Generalized Kullback-Leibler Divergence based.
%    9) Itakura-Saito Divergence based.
%    10) Jensen-Shannon Divergence based.
%--------------------------------------------------------------------------
% -> Data requirements: Maps, placed in the map folder 
%--------------------------------------------------------------------------
% -> See also: alg_genet_2D dist_est_2D dist_est_rob_2D fitness_2D 
%              initialization initiate_pop init_NP map_loading
%--------------------------------------------------------------------------
clc 

% The GL module works with different maps (placed in the maps folder). The
% user can choose one of them.
[map_real, map_known,num_map]=map_loading;
% The transposed map is used to estimate the sensor measurements.
map_known_tr=int8(map_known');
map_real_tr=int8(map_real'); 
% The map limits are computed
[m,n]=size(map_known_tr);
mapmax=[m,n,360];
mapmin=[1,1,0];

% The map is displayed in two different figures. 
figure(1); 
imagesc(map_real,'CDataMapping','scaled')
set(gcf,'Color','white');
set(gcf,'BackingStore','off');
set(gca,'DataAspectRatio',[1 1 1]);
colormap('gray')
hold on

figure(2);
imagesc(map_known,'CDataMapping','scaled')
set(gcf,'Color','white');
set(gcf,'BackingStore','off');
set(gca,'DataAspectRatio',[1 1 1]);
colormap('gray')
hold on

%--------------------------------------------------------------------------
%Initialization parameters:
NUM_MEASUREMENTS=61;        % Number of horizontal measurements in a scan.
SENSOR_RES=3*pi/180;        % Laser sensor angular resolution (radians)

if num_map==5.0
    CELL_SIZE=0.056;
elseif num_map==6.0
    CELL_SIZE=0.05;
else    
    CELL_SIZE=0.121;            % Cell size, in m.
end

T=1;                        % Translation constant. To estimate laser beams 
                            % in order to generate a laser scan,
                            % dist_est_2D considers increments of T units
                            % in the map.
SENSOR_RANGE=15/CELL_SIZE;  % Laser sensor range
D=3;                        % DE Number of Chromosomes. 
F=0.8;                     % Differential variations factor (mutation)
CR=0.5;                     % Crossover constant
% Variables introduced via keyboard
[position, err_dis, iter_max]=initialization(mapmin,mapmax,map_real');

%--------------------------------------------------------------------------
% The laser scan is generated according to the robot's true pose in the
% known map.
laser_real=dist_est_rob_2D(position,map_real_tr,mapmax,mapmin,err_dis,NUM_MEASUREMENTS,SENSOR_RES,SENSOR_RANGE,T);

%--------------------------------------------------------------------------
%   Initialization of the population size (NP). Two options: 
%                 1: Initialized by function init_NP
%                 Else: Fixed size given by code.
NP_opt=2;
if NP_opt==1
    NP=init_NP(laser_real,err_dis,mapmax,SENSOR_RES,NUM_MEASUREMENTS,CELL_SIZE);
else
    NP=100;
end
NP=round(NP);
fprintf(1,'\n Population size: %i \n',NP);

%--------------------------------------------------------------------------
% Different options for the GL algorithm can be selected via keyboard:
%  - DE Core Options:
%    1) Random Mutation, with Thresholding and Discarding (Default). 
%    2) Basic version, Random Mutation, without Thresholding, Discarding.
%    3) Mutation from Best candidate, with Thresholding and Discarding.
%    4) Random Mutation, with Thresholding and Discarding, NP is
%    drastically reduced (tracking) after convergence.
version_de=input('\ \n Introduce the DE version that you want to apply: \n 1) Random Mutation, with Thresholding and Discarding. \n 2) Basic version, Random Mutation, without Thresholding, Discarding. \n 3) Mutation from Best candidate, with Thresholding and Discarding. \n 4) Random Mutation, with Thresholding and Discarding, NP reduced (tracking) after convergence. \n');
if isempty(version_de),
    version_de=1;   
    fprintf(1,'\n \t Option 1 by default. \n');
end
%  - Fitness Function Options:
%    1) L2 Norm. Sum of the squared errors (Default).
%    2) L1 Norm. Sum of the absolute values of the error (Mahalanobis).
%    3) Kullback-Leibler Divergence based.
%    4) Density Power Divergence based.
%    5) Hellinger Distance based.
%    6) L2 Norm from Probability Distributions
%    7) L(Variable Exponent) Norm from Probability Distributions
%    8) Generalized Kullback-Leibler Divergence based.
%    9) Itakura-Saito Divergence based.
%    10) Jensen-Shannon Divergence based.
version_fitness=input('\ \n Introduce the Fitness Function that you want to apply: \n 1) L2 Norm. Sum of the squared errors (Default). \n 2) L1 Norm. Sum of the absolute values of the error. \n 3) Kullback-Leibler Divergence based. \n 4) Density Power Divergence based. \n 5) Hellinger Distance based. \n 6) L2 Norm from Probability Distributions \n 7) L(Variable Exponent) Norm from Probability Distributions. \n 8) Generalized Kullback-Leibler Divergence based. \n 9) Itakura-Saito Divergence based. \n 10) Jensen-Shannon Divergence based. \n');
if isempty(version_fitness),
    version_fitness=1;   
    fprintf(1,'\n \t Option 1 by default. \n');
end

N_SIMULATIONS=20;
Solution.best_estimate=zeros(N_SIMULATIONS,3);
Solution.error=zeros(N_SIMULATIONS,1);

for simul=1:N_SIMULATIONS
    
%--------------------------------------------------------------------------
% The initial population is randomly generated to cover the whole map.
D=3;
mapmax=[m,n,360];
mapmin=[1,1,0];
population=initiate_pop(mapmin,mapmax,NP,D);

%--------------------------------------------------------------------------
% Some parame change for the DPD and HC function, extra D is needed, with
% limits.
if (version_fitness==7)
    population=[population rand(NP,1)]; % Extra column with alpha
    D=4;
    mapmin=[mapmin 0.001]; % Lower limit for alpha is 0
    mapmax=[mapmax 1]; % Upper limit for alpha is 1
end

%--------------------------------------------------------------------------

% The robot motion simulation starts. In a single step, the robot tries to
% locate itself. After convergence, robot motion is allowed until an 'f'
% is introduced in dir_disp. In this case, the GL module ends its
% execution.
steps=0;
dir_disp=' ';
while (dir_disp~='f')
        
        fprintf(1,'\n Simulation: %d/%d ',simul,N_SIMULATIONS);
        fprintf(1,'\n Robot real pose  (x y theta) %f %f %f \n',position(1),position(2),position(3));
        tic
        % The DE-based GL filter is called
        [bestmem,error,population,F,NP]=alg_genet_2D(version_de,version_fitness,laser_real,map_known_tr,population,iter_max,mapmax,mapmin,err_dis,NUM_MEASUREMENTS,SENSOR_RES,NP,D,F,CR,SENSOR_RANGE,T);
        toc
        fprintf(1,'\n Robot real pose  (x y theta) %f %f %f %f:\n',position(1),position(2),position(3));
        fprintf(1,'\n Estimated pose given by the GL filter (x y theta) %f %f %f %f\n',bestmem(2),bestmem(3),bestmem(4));
        
        % The error between real pose and estimate is computed
        poserror=CELL_SIZE*100*sqrt((position(1)-bestmem(2))^2+(position(2)-bestmem(3))^2);
        orierror=abs(position(3)-bestmem(4));
        fprintf(1,'\n The position error is %f cm and the orientation error is %f grados\n',poserror,orierror);
        
        %----------------------------------------------------------------
        % Robot motion is allowed (after convergence) via keyboard 
        MOTION_TRANSL_RES=3;  % translation motion resolution in a step (cells)
        MOTION_ORIENT_RES=5;  % angular motion resolution in a step (degrees)
        disp=zeros(3);     % Vector that contains the displacements
        dir_disp='f';
        %dir_disp=input('\n Introduzce the movement direction: \n','s');
        if isempty(dir_disp)
            fprintf(1,'\n The default movement is zero \n');
        end
        if (dir_disp=='k')
            disp(1)=+MOTION_TRANSL_RES;
        elseif (dir_disp=='j')
            disp(1)=-MOTION_TRANSL_RES;
        elseif (dir_disp=='m')
            disp(2)=-MOTION_TRANSL_RES;
        elseif (dir_disp=='i')
            disp(2)=+MOTION_TRANSL_RES;
        elseif (dir_disp=='h')
            disp(3)=MOTION_ORIENT_RES;  %el angulo esta en radianes 
        elseif (dir_disp=='l')
            disp(3)=MOTION_ORIENT_RES;  
        end
        
        % The best solution is displaced according to the robot's motion.
        bestmem(2)=bestmem(2)+disp(1);
        bestmem(3)=bestmem(3)+disp(2);
        bestmem(4)=bestmem(4)+disp(3);

        % The real pose is displaced according to the robot's mootion, but
        % including a Normally distributed error with mean given by the
        % coordinate of disp and standard deviation equal to
        % disp*MOTION_ERROR. 
        MOTION_ERROR=0.03;      % Motion error standard dev, over disp.
        position(1)=position(1)+disp(1)*(1+MOTION_ERROR*randn(1));
        position(2)=position(2)+disp(2)*(1+MOTION_ERROR*randn(1));
        position(3)=position(3)+disp(3)*(1+MOTION_ERROR*randn(1));

         
        % The whole population is displaced, including the same type of
        % error.
        for i=1:NP
            population(i,2)= population(i,2)+disp(1)*(1+MOTION_ERROR*randn(1)); %+err_pos*randn(1);
            population(i,3)= population(i,3)+disp(2)*(1+MOTION_ERROR*randn(1)); %+err_pos*randn(1);
            population(i,4)= population(i,4)+disp(3)*(1+MOTION_ERROR*randn(1)); %+err_pos*randn(1);
            
            % It is checked that the population is inside the map limits. 
            for h=2:3
                if population(i,h)<mapmin(h-1) 
                    population(i,h)=mapmin(h-1);
                end         
                if population(i,h)>mapmax(h-1), 
                    population(i,h)=mapmax(h-1);
                end
            end
            if population(i,4)<mapmin(3)
                population(i,4)=population(i,4)+360.0;
            end
            if population(i,4)>mapmax(3)
                population(i,4)=population(i,4)-360.0;
            end
  
        end 
        
        %  The next option is used for tracking. NP is drastically reduced. 
        if version==6
            if error<50
                NP=10;
            end
        end
        
        steps=steps+1;
        
        % A new laser scan is read if the robot is in a new location.
        if dir_disp~='f'
            laser_real=dist_est_rob_2D(position,map_real_tr,mapmax,mapmin,err_dis,NUM_MEASUREMENTS,SENSOR_RES,SENSOR_RANGE,T);
        end
        
end

% The best solution and the error are returned.
Solution.best_estimate(simul,:)=bestmem(2:(D+1));
Solution.error(simul)=error;

end

%--------------------------------------------------------------------------
% Representation of results

% Obtaining coordinates of redings from true location
observations=zeros(1,NUM_MEASUREMENTS);
robot=position(1)+1i*position(2);
ang_robot=(position(3)-90)*pi/180;% se pasa a radianes
for index=1:NUM_MEASUREMENTS
    observations(index)=robot+laser_real(index)*exp(1i*(ang_robot+(index-1)*SENSOR_RES));
end
% Obtaining coordinates of redings from best estimate
est_observations=zeros(1,NUM_MEASUREMENTS);
est_robot=bestmem(2)+1i*bestmem(3);
est_ang_robot=(bestmem(4)-90)*pi/180;% se pasa a radianes
for index=1:NUM_MEASUREMENTS
    est_observations(index)=est_robot+laser_real(index)*exp(1i*(est_ang_robot+(index-1)*SENSOR_RES));
end

% Display robot's position
figure(1)
plot(position(1),position(2),'b+','MarkerSize',5);
% Display observations from true location
for i=1:NUM_MEASUREMENTS  
    plot(real(observations(i)),imag(observations(i)),'m.','MarkerSize',5);
end 
% Display population set
figure(2)
plot(population(:,2),population(:,3),'r.','MarkerSize',5);        
% Display observations from best estimate
for i=1:NUM_MEASUREMENTS  
    plot(real(est_observations(i)),imag(est_observations(i)),'m.','MarkerSize',5);
end 

end