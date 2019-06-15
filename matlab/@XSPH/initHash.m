function initHash(o)
%INITHASH Initialize info for hashing

% Get geometry, and set xy limits
cor = [o.fluid.r; o.bound.r];
o.locMax = max(cor, [], 1);
o.locMin = min(cor, [], 1);
o.hashOffset = -1 + floor(o.locMin/o.sup);
o.hashSubDim = floor(o.locMax/o.sup) - o.hashOffset;

% Get the boundaries of hash states (for future plotting)
o.hashTicks = cell(o.dim, 1);
for d = 1:o.dim
    o.hashTicks{d} = ((1:o.hashSubDim(d))+o.hashOffset(d))*o.sup;
end

% Get adjacent partitions
o.hashAdj = aux.get_adjacent(o.hashSubDim);

% Calculate boundary particle hashes
o.bound.getHash(o.sup, o.hashOffset, o.hashSubDim, o.hashAdj);


end

