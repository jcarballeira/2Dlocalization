function [pop]=initiate_pop(mapmin,mapmax,NP,D)
%--------------------------------------------------------------------------
%   Function: initiate_pop
%   Author: Fernando Martin Monar.
%   Date: November, 2015
%--------------------------------------------------------------------------
% -> Description: Initial population generation. A random population of NP
% candiates, each one with D chromosomes, in randomly generated to cover
% the whole map according to the map limits (mapmin and mapmax).
%--------------------------------------------------------------------------
% -> Inputs:
%       -mapmin: Minimum index in the map. Typically =[1,1,0].
%       -mapmax: Vector of 3 elements that corresponds to the map size. The
%       first two coordinates are the map dimensions, in Cartesian
%       coordinates, and the third one is the orientation, typically 360
%       degrees.
%       -NP: population size.
%       -D:Number of chromosomes
%--------------------------------------------------------------------------
% -> See also: Global_Local_2D
%--------------------------------------------------------------------------

poppre=zeros(NP,D);
cost=zeros(NP,1);
for i=1:NP
   poppre(i,:) = mapmin + rand(1,D).*(mapmax- mapmin);
end
for i=1:NP
    for j=1:D
        if poppre(i,j)<mapmin(j), poppre(i,j)=mapmin(j);end
        if poppre(i,j)>mapmax(j), poppre(i,j)=mapmax(j);end
    end
end
pop=[cost poppre];
end