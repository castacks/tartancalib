/**
* Copyright (C) 2017-present, Facebook, Inc.
*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#pragma once

#ifndef INCLUDE_CHECKERBOARD_DETECTOR_TRIANGLE_H
#define INCLUDE_CHECKERBOARD_DETECTOR_TRIANGLE_H

#include <utility>
#include <vector>

namespace orp {
namespace calibration {

template <typename PointType>
class Triangle {
public:
  typedef PointType point_type;
  //! Constructor
  Triangle(const std::vector<PointType> &p) {
    p0 = p[0];
    p01 = p[1] - p[0];
    p02 = p[2] - p[0];
  }

  //! Same as interpolate, except that the coordinates are interpreted between 0
  //! and 1, instead of -1 and 1.
  PointType interpolate01(float x, float y) const {
    return p0 + p01 * x + p02 * y;
  }

private:
  PointType p0, p01, p02;
};

} // namespace calibration
} // namespace orp

#endif
