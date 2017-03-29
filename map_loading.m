function [map_real, map_known,num_map]=map_loading
%--------------------------------------------------------------------------
%   Function: map_loading
%   Author: Fernando Martin Monar.
%   Date: November, 2015
%--------------------------------------------------------------------------
% -> Description: The GL module works with different maps (placed in the
% maps folder). The user can choose between them via keyboard.
%--------------------------------------------------------------------------
% -> Inputs: no inputs are required.
% -> Output: 
%       -map_real: This map will be used to obtain the true measurements
%       (laser readings from the robot's true location)
%       If the user want to work with dynamic obstacles or occlusions,
%       these features can be included in this map, the GL module will not
%       consider these features because the known map will be used to
%       optmize       
%       -map_known: Learned map used by the GL module. The estimates will
%       use this map to simulate the laser readings.
%--------------------------------------------------------------------------
% -> See also: Global_Local_2D
%--------------------------------------------------------------------------

type_map=input('\n Introduce the map type: \n \ \ \ 1 - Total \n \ \ \ 2 - Partial \n \ \ \ 3 - Test \n \ \ \ 4 - Huge \n \ \ \ 5 - Real \n \ \ \ 6 - Real Intel \n');
if isempty(type_map)
    type_map=1.0;   
    fprintf(1,'\n \t The default map is Total \n');
end
if type_map==1.0
  map_known=imread('maps/uc3m_p3_open','bmp'); 
  map_real=imread('maps/uc3m_p3_obstacles43','bmp');  
  %map_real=imread('maps/uc3m_p3_obstacles','bmp');
elseif type_map==2.0
    map_known=imread('maps/uc3m_p3_parcial_open','bmp');
    map_real=imread('maps/uc3m_p3_parcial_open','bmp');
elseif type_map==3.0
    map_known=imread('maps/E_test_1_open','bmp');
    map_real=imread('maps/E_test_1_open','bmp');
elseif type_map==4.0
    map_known=imread('maps/e_test_3_open','bmp');
    map_real=imread('maps/e_test_3_open','bmp');
elseif type_map==5.0
    map_known=imread('maps/uc3m_real_map','bmp');
    map_real=imread('maps/uc3m_real_map','bmp');
elseif type_map==6.0
    map_known=imread('maps/intel3DmapV3','png');
    
    [nf,nc,plane]=size(map_known);
    for i1=1:nf
        for j1=1:nc
            if (map_known(i1,j1,1) > 230)
                mapa(i1,j1)=1.0;
            else
                mapa(i1,j1)=0.0;
            end;
        end
    end
    map_known=mapa;
    map_real=mapa;
else
  map_known=imread('maps/uc3m_p3_open','bmp');  
  map_real=imread('maps/uc3m_p3_open','bmp');
end  
save('map_known');
save('map_real');
load('map_known');
load('map_real');
num_map=type_map;

end