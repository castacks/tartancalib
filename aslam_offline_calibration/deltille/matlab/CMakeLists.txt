## Copyright (C) 2017-present, Facebook, Inc.
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
##
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY)
include(${CMAKE_SOURCE_DIR}/cmake/mexmat.cmake)

matlab_add_mex(
  NAME deltille_detector_mex
  SRC deltille_detector_mex.cc
  LINK_TO deltille mexmat)

matlab_add_mex(
  NAME saddle_subpix_refinement_mex
  SRC saddle_subpix_refinement_mex.cc
  LINK_TO deltille mexmat)

matlab_add_mex(
  NAME detect_calib_boards_mex.cc
  SRC detect_calib_boards_mex.cc
  LINK_TO deltille mexmat)




