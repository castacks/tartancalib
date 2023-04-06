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
#ifndef INCLUDE_CHECKERBOARD_DETECTOR_QUADANGLE_H
#define INCLUDE_CHECKERBOARD_DETECTOR_QUADANGLE_H

#include <utility>
#include <vector>

namespace orp {
namespace calibration {

template <typename PointType> class Quadangle {
public:
  typedef PointType point_type;

  //! Constructor
  Quadangle(const std::vector<PointType> &p) {
    p0 = p[0];
    p3 = p[3];
    p01 = p[1] - p[0];
    p32 = p[2] - p[3];
  }

  //! Interpolate given that the lower left corner of the lower left cell is at
  //! (-1,-1) and the upper right corner of the upper right cell is at (1,1).
  PointType interpolate(float x, float y) const {
    PointType r1 = p0 + p01 * ((x + 1.0f) / 2.0f);
    PointType r2 = p3 + p32 * ((x + 1.0f) / 2.0f);
    return r1 + (r2 - r1) * ((y + 1.0f) / 2.0f);
  }

  //! Same as interpolate, except that the coordinates are interpreted between 0
  //! and 1, instead of -1 and 1.
  PointType interpolate01(float x, float y) {
    return interpolate(2 * x - 1, 2 * y - 1);
  }

private:
  PointType p0, p3, p01, p32;
};

} // namespace calibration
} // namespace orp

#endif
