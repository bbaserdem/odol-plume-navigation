function [C, A] = get_adjacent(shp)
% GETADJACENT Get linear index of grid points that are adjacent in D space
%   Given a shape, (of L=prod(shp) points in D=size(shp,1) space) :
%       A=(L,3^D) where the row A(i,:) contains linear indices which are
%       adjacent to the i'th linear index in D-space

% Get a row which gives the linear index displacement
S = [1, cumprod(shp(1:(end-1)))];
% Generate list of positive change vectors
V = aux.get_positiveVec(length(shp));
% Calculate the row of linear index displacement for the ones in V
D = S * V';
% Generate the subscrip index for all the voxels
O = aux.i2s(shp,(1:prod(shp))');
% Calculate the linear index of all adjacents
A = ((1:prod(shp))') + D;

% Null the invalid indices;
for d = 1:length(shp)
    % When dim index;O(<i>,d) was max;shp(d) and tried to inc;V(<j>)== 1
    A(O(:,d)==shp(d),V(:,d)== 1) = 0;
    % When dim index;O(<i>,d) was min;1      and tried to dec;V(<j>)==-1
    A(O(:,d)==1,     V(:,d)==-1) = 0;
end

% Create a cell giving the index in rows, and trimming the access
C = cellfun(@(x) x(x>0), num2cell(A,2), 'UniformOutput', false);

end
