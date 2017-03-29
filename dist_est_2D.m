function [laser_scan]=dist_est_2D(position,map,mapmax,mapmin,err_dis,NUM_MEASUREMENTS,SENSOR_RES,SENSOR_RANGE,T)
%--------------------------------------------------------------------------
%   Function: dist_est_2D
%   Author: Fernando Martin Monar, Luis Moreno
%   Date: November, 2015
%--------------------------------------------------------------------------
% -> Description: a 2D laser scan is estimated for a given position in a
% known map.
%--------------------------------------------------------------------------
% -> Inputs:
%       -position: Pose from which the measurements are estimated. Vector
%       coordinates are given in Cartesian coordinates and orientation (3
%       dof).
%       -map: 2D Map of the environment. Known map to estimate readings. In
%       map, 0 is an obstacle, 1 represents free space and 0.5 is unknown.
%       -mapmax: Vector of 3 elements that corresponds to the map size. The
%       first two coordinates are the map dimensions, in Cartesian
%       coordinates, and the third one is the orientation, typically 360
%       degrees.
%       -mapmin: Minimum index in the map. Typically =[1,1,0].
%       -err_dis: Sensor noise, standard deviation, in percentage over the
%       distance weighted.
%       -NUM_MEASUREMENTS: Number of horizontal measurements in a scan.
%       -SENSOR_RES: Laser sensor angular resolution (radians).
%       -T: Translation constant. To estimate laser beams in order to
%       generate a laser scan, this function considers increments of T
%       units in the map.
%--------------------------------------------------------------------------
% -> Outputs: 
%       -laser_scan: Vector with NUM_MEASUREMENTS elements containing the
%       distances of the laser measurements.
%--------------------------------------------------------------------------
% -> File requirements: this function is called by Global_Local_2D.m
%--------------------------------------------------------------------------
% -> See also: Global_Local_2D
%--------------------------------------------------------------------------
laser_scan=zeros(1,NUM_MEASUREMENTS);

ths=+90*pi/180; % Angle of the first measurement.

for j=1:NUM_MEASUREMENTS
    
    x=position(1);
    y=position(2);
    thr=(-position(3)+90)*pi/180;   % from deg to radians   
    
    sin_sensor=sin(thr+ths);
    cos_sensor=cos(thr+ths);
    incr_x=T*sin_sensor;
    incr_y=T*cos_sensor;
    
    final=0;
    dis=0;

    while(final < 2)
        
        x_round=round(x);
        y_round=round(y);

        if map(x_round,y_round)==0
            if dis==0,
              dis=sqrt((x-position(1))^2+(y-position(2))^2);
              final=final+1;
            end           
        end    
        
        x=x+incr_x;
        y=y+incr_y;  
         
        if x>mapmax(1), final=2;end
        if x<mapmin(1), final=2;end
        if y>mapmax(2), final=2;end
        if y<mapmin(2), final=2;end    
    end
    
    if dis<0,dis=0;end
    if dis>SENSOR_RANGE, dis=SENSOR_RANGE;end

    laser_scan(j)=dis;
    ths=ths-SENSOR_RES; % Change of orientation between measurements.
end 

% Normally distributed noise is added to the measurements. 
%err_m=randn(size(laser_scan));    % A N(0,1) distribution is generated.
%err_level=laser_scan.*err_dis;    % The standard deviation is computed for 
%                                   each measurement.
%laser_scan=laser_scan+(err_m.*err_level);   % The final measurement 
%                   contains noise with a standard devation proportional to
%                   the distance weighted.

end