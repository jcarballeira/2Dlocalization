function[position, err_dis, iter_max]=initialization(mapmin,mapmax,map)
%--------------------------------------------------------------------------
%   Function: initialization
%   Author: Fernando Martin Monar.
%   Date: November, 2015
%--------------------------------------------------------------------------
% -> Description: The initial configuration is introduced via keyboard
%--------------------------------------------------------------------------
% -> Inputs:
%       -mapmin: Minimum index in the map. Typically =[1,1,0].
%       -mapmax: Vector of 3 elements that corresponds to the map size. The
%       first two coordinates are the map dimensions, in cartesian
%       coordinates, and the third one is the orientation, typically 360
%       degrees.
%       -map: matrix that contains the 2D map of the environment.
%--------------------------------------------------------------------------
% -> Outputs: 
%       -position: robot's true pose. Vector coordinates are given in
%       Cartesian coordinates and plane orientation (3 dof).
%       -err_dis: sensor noise, standard deviation, in percentage over the
%       distance weighted.
%       -iter_max: maximum number of iterations.
%--------------------------------------------------------------------------
% -> File requirements: this function is called by Global_Local_2D.m
%--------------------------------------------------------------------------
% -> See also: Global_Local_2D
%--------------------------------------------------------------------------

% The robot initial orientation is introduced.
th0=input('\n Introduce the robot orientation (in degrees):');
if th0>=360,th0=th0-360;end
if th0<0,th0=th0+360;end 
if isempty(th0),    
    th0=0;  
    fprintf(1,'\n The default orientation is %.2fº \n',th0);
end
%--------------------------------------------------------------------------
position=[1,1,th0]; % th0 in degrees
%--------------------------------------------------------------------------
% The robot Cartesian coordinates, limited by the map borders, are read.
fprintf(1,'\n Introduce the robot coordinates:')
while(map(round(position(1)),round(position(2)))==0),
     position(1)=input('\n Introduzce the x coordinate:');
     if(position(1)<mapmin(1)||position(1)>mapmax(1)),
         fprintf(1,'\n Map limits in x are [%.2f,%.2f]:',mapmin(1),mapmax(1));
         position(1)=input('\n Please, introduce a different x coordinate:');
     end
     position(2)=input('\n Introduzce the y coordinate: ');
     if(position(2)<mapmin(2)||position(2)>mapmax(2)),
        fprintf(1,'\n Map limits in y are [%.2f,%.2f] ',mapmin(2),mapmax(2));
         position(2)=input('\n \ Please, introduce a different y coordinate:');
     end
     if(map(round(position(1)),round(position(2)))==0),
         fprintf(1,'\n This position corresponds to an obstacle in the map. \n You have to introduce a different location.\n');
     end
end
%--------------------------------------------------------------------------
%The Laser sensor error is chosen
err_dis=input('\n Introduce the laser sensor error, in % (standard deviation over the distance weighted): \ \');
err_dis=err_dis/100;
if isempty(err_dis),
    err_dis=0.1;   
    fprintf(1,'\n \t The default error is %.2f%%  \n',err_dis*100);
end
%--------------------------------------------------------------------------
% Genetic algorithm upper iterations limit
iter_max=input('\ \n Introduzce the maximum number of iterations: \ \');
iter_max=round(iter_max);
if isempty(iter_max),
    iter_max= 40;   
    fprintf(1,'\n \t The default iterations are %d \n',iter_max);
end
end