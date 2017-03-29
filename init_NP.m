function[NP]=init_NP(laser_real,err_dis,mapmax,SENSOR_RES,NUM_MEDIDAS,CELL_SIZE)
%--------------------------------------------------------------------------
%   Main Function: init_NP
%   Author: Fernando Martin Monar.
%   Date: November, 2015
%--------------------------------------------------------------------------
% -> Description: The population size is computed according to the
% information contained in a laser scan. The algorithm considers size,
% overlapping, symmetries... to estimate an optimum population size.
%--------------------------------------------------------------------------
% -> Inputs:
%       -laser real: vector with NUM_MEDIDAS components containing the
%       distances of the laser measurements.
%       -err_dis: sensor noise, standard deviation, in percentage over the
%       distance weighted.
%       -mapmax: Vector of 3 elements that corresponds to the map size. The
%       first two corrdinates are the map dimensions, in Cartesian
%       coordinates, and the thirs one is the orientation, typically 360
%       degrees.
%       -SENSOR_RES: Laser sensor angular resolution (radians).
%       -NUM_MEASUREMENTS: Number of horizontal measurements in a scan.
%       -CELL_SIZE: Cell size.
% -> Output: 
%       -NP: population size.
%--------------------------------------------------------------------------
% -> See also: Global_Local_2D
%--------------------------------------------------------------------------

%Overlapping
l=1/tan(SENSOR_RES*pi/180);
Aef=zeros(NUM_MEDIDAS); 
Aef(1)=laser_real(1);
nd=0; % Number of discontinuities
%--------------------------------------------------------------------------

for j=2:NUM_MEDIDAS

    if (laser_real(j)>l) && (laser_real(j-1)>l)
        Aef(j)=CELL_SIZE*(laser_real(j)-(l/2));
    else
        Aef(j)=CELL_SIZE*(laser_real(j)-(l/2)+(l-min(laser_real(j),laser_real(j-1)))^2/(2*l));
    end
    
    %Discontinuities
    dif=laser_real(j)-laser_real(j-1);
    dif_percent=(laser_real(j)-laser_real(j-1))/laser_real(j)*100;
    if (dif>50)||(dif_percent>60)
        nd=nd+1;
    end
end

NP=mapmax(1)*mapmax(2)/(sum(sum(Aef)));
NP=(1+sqrt(err_dis))*NP;

end

% Previouse version. Check

% %------------------------------------------------------------------------
% % % Determinacion automatica de NP
%  Total=sum(v0);
%  d_m=Total/num_sensores;
%  d_max=15/0.121;               % d_max= max_range/cell_size
%  i_sp=exp(-(d_m/d_max))  ;        % indice de solapamiento
%  area_B=Total*(1-i_sp);
%  area_T=dim_x*dim_y;
% % %NP=round((0.75*area_T)/area_B)
% 
% % log2 logaritmo en base 2
% % log  logaritmo natural en base e
% % log10 logaritmo en base 10
% % viejos
% %NP=round(area_T/1000)+round((area_T/area_B)*log10(area_T/area_B))
% %NP=round(pi*log(area_T))+round(0.5*(area_T/area_B)*log10(area_T/area_B))
% 
% %------------------------------------------------------------------------
% % Determinacion automatica de NP
% %
% % Nmax=2^((n/2)log(1+St/So))
% % n=dimensiones
% % S_t=superficie total
% % S_o= superficie observada
% %
% Dim=3;
% NP=round(0.5*2^((Dim/2)*log(1+area_T/area_B)))
