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

classdef TargetDetector < handle

  properties (SetAccess = private, Hidden = true)
    handle_;
  end

  methods
    function this = TargetDetector(dsc_file, w, h)
      assert(isempty(dsc_file) || exist(dsc_file, 'file') > 0);
      if nargin < 3, h = 22; end
      if nargin < 2, w = 22; end

      this.handle_ = target_detector_mex('new', dsc_file, w, h);
    end

    function delete(this)
      target_detector_mex('delete', this.handle_);
    end

    function [corners, board_ids, is_ordered] = detectCorners(this, I)
      %function [corners, board_ids, is_orderd] = detectCorners(I)
      %
      % Runs the corner detector

      if size(I,3) > 1
        I = rgb2gray(I);
      end

      [corners, board_ids, is_ordered] = target_detector_mex(...
        'run', this.handle_, I);
    end

  end

end % TargetDetector

