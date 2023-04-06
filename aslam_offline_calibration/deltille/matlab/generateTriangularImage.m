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

function img=generateTriangularImage(row, col, uv_grid, line1, line2, line3, uv_grid_sample)
    
    n_sample=size(uv_grid_sample,2);
    img=min(max(line1*uv_grid,-sqrt(0.5)),sqrt(0.5))/sqrt(0.5).*min(max(line2*uv_grid,-sqrt(0.5)),sqrt(0.5))/sqrt(0.5).*min(max(line3*uv_grid,-sqrt(0.5)),sqrt(0.5))/sqrt(0.5);
    img=reshape(img,[row,col]);
    
    % boundary
    id=find(img>-1&img<1);
    n_pixel=numel(id);
    u=uv_grid(1,id)';
    v=uv_grid(2,id)';
    u=repmat(u,1,n_sample)+repmat(uv_grid_sample(1,:),n_pixel,1);
    v=repmat(v,1,n_sample)+repmat(uv_grid_sample(2,:),n_pixel,1);
    intensity=min(max(line1(1)*u+line1(2)*v+line1(3),-sqrt(0.5)),sqrt(0.5))/sqrt(0.5).*min(max(line2(1)*u+line2(2)*v+line2(3),-sqrt(0.5)),sqrt(0.5))/sqrt(0.5).*min(max(line3(1)*u+line3(2)*v+line3(3),-sqrt(0.5)),sqrt(0.5))/sqrt(0.5);
    intensity=sum(intensity,2)/n_sample;
    img(id)=intensity;
end
