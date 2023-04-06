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

#include "mexmat.h"
#include <deltille/DetectorTools.h>
#include <deltille/GridDetectorContext.h>

#include <opencv2/core/core.hpp>

#include <iostream>
#include <vector>

static constexpr bool kFindBestOnly = false;

typedef uint8_t DefaultInputImageType;
typedef double DefaultFloatImageType;

template <class SaddlePointType,
          typename InputImageType = DefaultInputImageType,
          typename FloatImageType = DefaultFloatImageType>
void findBoards(const cv::Mat &I, const cv::Size &board_size,
                std::vector<orp::calibration::BoardObservation> &boards,
                bool best_only = kFindBestOnly) {
  using namespace orp::calibration;

  GridDetectorContext<SaddlePointType, InputImageType, FloatImageType> gd(I);
  gd.findBoards(board_size, boards, best_only);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  static const char *USAGE = "[x,b,i] = fn(I, w, h, t)";
  if (nrhs < 4) {
    mexError("USAGE: %s\n", USAGE);
  }

  const mex::Mat<uint8_t> I(prhs[0]);
  int board_width = mex::getNumber<int>(prhs[1]),
      board_height = mex::getNumber<int>(prhs[2]),
      is_triangular = mex::getNumber<int>(prhs[3]) == 1;

  std::vector<orp::calibration::BoardObservation> boards;
  const cv::Mat I_cv = mex::mex2cv<DefaultInputImageType>(I);

  cv::Size board_size(board_width, board_height);

  if (is_triangular) {
    findBoards<orp::calibration::MonkeySaddlePoint>(I_cv, board_size, boards);
  } else {
    findBoards<orp::calibration::SaddlePoint>(I_cv, board_size, boards);
  }

  mex::Cell x(1, boards.size());
  mex::Cell b_ids(1, boards.size());
  mex::Cell indexed(1, boards.size());
  for (size_t i = 0; i < boards.size(); ++i) {
    auto n = int(boards[i].corner_locations.size());
    mex::Mat<double> x_i(2, n);
    mex::Mat<int> b_i(1, n);
    mex::Mat<bool> indexed_i(1, n);
    for (int j = 0; j < n; ++j) {
      x_i(0, j) = boards[i].corner_locations[j].x;
      x_i(1, j) = boards[i].corner_locations[j].y;
      b_i[j] = boards[i].board_id;
      indexed_i[j] = boards[i].indexed;
    }

    x.set(i, x_i.release());
    b_ids.set(i, b_i.release());
    indexed.set(i, indexed_i.release());
  }

  if (nlhs > 0)
    plhs[0] = x.release();
  if (nlhs > 1)
    plhs[1] = b_ids.release();
  if (nlhs > 2)
    plhs[2] = indexed.release();
}
