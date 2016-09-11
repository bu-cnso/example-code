function x = trimEnds(x, N, dim)
%% trimEnds
% Author: Erik Roberts
%
% Purpose: trims start and end of data x by N along dim (default: 1)
%
% Usage: x = trimEnds(x, N)
%        x = trimEnds(x, N, dim)
%
% Inputs (required):
%   x: data to trim
%   N: amount to trim
%
% Inputs (optional):
%   dim: dim of indices to trim (default: 1)

if isvector(x)
  x = x(N+1:end-N);
else
  dim = setDefault('dim', 1);
  
  Index(1:ndims(x)) = {':'};
  xDimLen = size(x,dim);
  Index{dim} = N+1:xDimLen-N;
  x = x(Index{:});
end

end
