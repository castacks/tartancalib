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

classdef SyntheticCornerGenerator < handle

  properties (SetAccess = private)
    nrows_;
    ncols_;

    uv_grid_;
    uv_grid_sample_;
    cx_;
    cy_;
  end

  properties(SetAccess = private, Hidden)
    S_;
  end

  methods
    function this = SyntheticCornerGenerator(nrows, ncols)
      if nargin < 2, ncols = 200; end
      if nargin < 1, nrows = 200; end

      this.nrows_ = nrows;
      this.ncols_ = ncols;

      this.cx_ = (1 + ncols) / 2.0;
      this.cy_ = (1 + nrows) / 2.0;

      [u, v] = meshgrid((1:ncols) - this.cx_, (1:nrows) - this.cy_);
      this.uv_grid_ = [u(:), v(:), ones(nrows*ncols, 1)]';

      delta = 0.05;
      [u, v] = meshgrid(-0.5:delta:0.5, -0.5:delta:0.5);
      this.uv_grid_sample_ = [u(:) v(:)]';

      this.S_{1} = [1,0; 0,1; 0,0; 1,1];
      this.S_{2} = [1,cos(60/180*pi),cos(120/180*pi);
             0,sin(60/180*pi),sin(120/180*pi);
             0,0,0;
             1,1,1];
    end

    function x = getInitialPoints(this)
      [u,v] = meshgrid((this.cx_ - 0.5 : this.cx_ + 0.5), ...
        (this.cy_ - 0.5) : (this.cy_ + 0.5));
      x = [u(:), v(:)]';
    end

    function [I, xt] = makeImage(this, pattern_type, t, p, b, n)
      %function [I, xt] = makeImage(this, pattern_type, t, p, b, n)
      %

      dx = rand - 0.5;
      dy = rand - 0.5;

      xt = [this.cx_ + dx; this.cy_ + dy];

      rz=[-sin(t)*cos(p); -sin(t)*sin(p); -cos(t)*ones(size(p))];
      tmp=rz+[zeros(2,size(rz,2));ones(1,size(rz,2))];
      ry=cross(rz, tmp);
      ry=ry./repmat(sqrt(sum(ry.^2)),3,1);
      rx=cross(ry,rz);

      R=[rx,ry,rz]';
      T=[0;0;100];
      RT=[R,T];

      make_square = (isnumeric(pattern_type) && (pattern_type == 1)) || ...
        (ischar(pattern_type) && strcmp('square', pattern_type));
      if make_square
        I = this.makeSquare(RT, b, n, dx, dy);
      else
        I = this.makeTriangle(RT, b, n, dx, dy);
      end


      white_level = 0.9;
      black_level = 0.1;

      I = I/2 + 0.5;
      I = I*(white_level - black_level) + black_level;

      a=rand(1)*2*pi;
      d=reshape([sin(a),-cos(a),0]*this.uv_grid_, size(I));
      L = ((d-min(d(:)))/(max(d(:))-min(d(:)))-0.5)*rand(1)*0.2+1.0;

      I = I .* L;

      if b > 0
        bsize = ceil(b*6)+1-mod(ceil(b*6),2);
        k = fspecial('gaussian', [bsize bsize], b);
        I = imfilter(I, k);
      end

      if n > 0
        I = imnoise(I, 'gaussian', 0, n^2);
      end
    end
  end

  methods (Access = private)
    function [I, xt] = makeSquare(this, Rt, b, n, dx, dy)

      xsaddle=Rt*this.S_{1};
      vsaddle=xsaddle./repmat(xsaddle(3,:),3,1);
      asaddle=[
      atan2(vsaddle(2,1),vsaddle(1,1)); ...
        atan2(vsaddle(2,2),vsaddle(1,2))];

      line1=[sin(asaddle(1)),-cos(asaddle(1)),-dx*sin(asaddle(1))+dy*cos(asaddle(1))];
      line2=[sin(asaddle(2)),-cos(asaddle(2)),-dx*sin(asaddle(2))+dy*cos(asaddle(2))];

      I = generateRectangularImage(this.nrows_, this.ncols_, this.uv_grid_, ...
        line1, line2, this.uv_grid_sample_);

    end

    function [I, xt] = makeTriangle(this, Rt, b, n, dx, dy)

      xmonkey=Rt*this.S_{2};
      vmonkey=xmonkey./repmat(xmonkey(3,:),3,1);
      amonkey=[atan2(vmonkey(2,1),vmonkey(1,1));
            atan2(vmonkey(2,2),vmonkey(1,2));
            atan2(vmonkey(2,3),vmonkey(1,3))];

      line1=[sin(amonkey(1)),-cos(amonkey(1)),-dx*sin(amonkey(1))+dy*cos(amonkey(1))];
      line2=[sin(amonkey(2)),-cos(amonkey(2)),-dx*sin(amonkey(2))+dy*cos(amonkey(2))];
      line3=[sin(amonkey(3)),-cos(amonkey(3)),-dx*sin(amonkey(3))+dy*cos(amonkey(3))];

      I = generateTriangularImage(this.nrows_, this.ncols_, this.uv_grid_, ...
        line1, line2, line3, this.uv_grid_sample_);
    end
  end 

end
