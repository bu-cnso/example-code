function y = reflectEnds(x, N, dim)
%% reflectEnds
% Author: Erik Roberts
%
% Purpose: reflects data at start and end of data x by N along dim (default: 1)
%
% Usage: x = reflectEnds(x, N)
%        x = reflectEnds(x, N, dim)
%
% Inputs (required):
%   x: vector to augment
%   N: amount to reflect on each end
%
% Inputs (optional):
%   dim: dim to work along (default: 1)

if isrow(x)
  xStartReflect = fliplr(x(1, 2:N+1));
  xEndReflect = fliplr(x(1, end-N:end-1));
  y = [xStartReflect, x, xEndReflect];
elseif iscolumn(x)
  xStartReflect = flipud(x(2:N+1, 1));
  xEndReflect = flipud(x(end-N:end-1, 1));
  y = [xStartReflect; x; xEndReflect];
else
  dim = setDefault('dim', 1);
  
  xDimLen = size(x,dim);
  
  Index(1:ndims(x)) = {':'};
  
  xStartReflectIndex = Index;
  xStartReflectIndex{dim} = 2 : N+1;
  xStartReflect = x(xStartReflectIndex{:});
  y = cat(dim, xStartReflect, x);
  clear xStartReflectIndex xStartReflect
  
  xEndReflectIndex = Index;
  xEndReflectIndex{dim} = xDimLen-N  : xDimLen-1;
  xEndReflect = x(xEndReflectIndex{:});
  y = cat(dim, y, xEndReflect);
  clear xEndReflectIndex xEndReflect
end

end
