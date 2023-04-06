% Copyright (C) 2017-present, Facebook, Inc.
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

function v = randrange(min_v, max_v, siz)
  % function v = randrange(min_v, max_v)
  %
  % Return a random number in the interval min_v max_v
  %
  if nargin < 3, siz = [1 1]; end

  v = (max_v - min_v).*rand(siz) + min_v;
end
