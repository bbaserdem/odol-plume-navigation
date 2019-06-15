function getHash(o, l, off, dim, adj)
%GETHASH Calculate the hash of the particles
%   l:   the hashing distance
%   off: offset to be used with the hash
%   dim: maximum dimensionality expected

% Get the subscript hash
if ~exist('off', 'var')
    off = -1 + min(floor(o.r/l), [], 1);
end
o.subHash = floor(o.r/l) - off;

% Get the linear hash
if ~exist('dim', 'var')
    dim = max(o.subHash, [], 1);
end
o.linHash = aux.s2i(dim, o.subHash);

% Seperate to partitions
o.parHash = aux.get_partitions(o.linHash, prod(dim));
o.parLen = cellfun(@length, o.parHash);

% Do adjaceny if info is given
if exist('adj', 'var')
    o.adjHash = cellfun( @(x) [o.parHash{x}], adj, ...
        'UniformOutput', false);
    o.adjLen = cellfun(@length, o.adjHash);
end

end

