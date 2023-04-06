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

function [x, b_id, indexed] = detectCalibrationBoards(I, board_width, board_height, is_triangular)
  % function [x, b_id, indexed] = detectCalibrationBoards(I, board_width, board_height, is_triangular)
  %
  % INPUT
  %   I   the input image
  %   board_width, board_height the size of the calibration board
  %   is_triangular  we are looking for triangular boards
  %
  % OUTPUT
  %   x     [u,v] corner locations (2xN)
  %   b_id  unique board id (1xN)
  %   indexed true if the corner was indexed properly (1xN)
  %
  
  if size(I,3) > 1
    % we'll assume rgb
    I = rgb2gray(I);
  end


  [x, b_id, indexed] = detect_calib_boards_mex(I, int32(board_width), int32(board_height), int32(is_triangular));

end % detectCalibrationBoards
