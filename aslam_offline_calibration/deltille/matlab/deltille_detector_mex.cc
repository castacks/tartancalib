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
#include <deltille/target_detector.h>

#include <string>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  if (nrhs < 1)
    mex::error("must have at the least one argument");

  std::string command = mex::getString(prhs[0]);

  if ("new" == command) {
    mex::nargchk(4, 4, nrhs, "h = cmd('new', dsc_file, w, h)");
    auto dsc_file = mex::getString(prhs[1]);
    auto w = mex::getNumber<int>(prhs[2]);
    auto h = mex::getNumber<int>(prhs[3]);
    plhs[0] = mex::PtrToMex<TargetDetector>(new TargetDetector(dsc_file, w, h));
  } else if ("delete" == command) {
    mex::DeleteClass<TargetDetector>(prhs[1]);
  } else if ("run" == command) {
    auto ptr = mex::MexToPtr<TargetDetector>(prhs[1]);
    mex::Mat<uint8_t> I(prhs[2]);
    auto I_cv = mex::mex2cv(I);

    std::vector<CalibrationCorner> corners;
    ptr->run(I_cv, corners);

    if (nlhs > 0) {
      mex::Mat<double> x(2, mwSize(corners.size()));
      for (int i = 0; i < int(corners.size()); ++i) {
        x(0, i) = corners[i].x;
        x(1, i) = corners[i].y;
      }

      plhs[0] = x.release();
    }

    if (nlhs > 1) {
      mex::Mat<int> board_id(1, corners.size());
      for (int i = 0; i < int(corners.size()); ++i) {
        board_id[i] = corners[i].boardId;
      }

      plhs[1] = board_id.release();
    }

    if (nlhs > 2) {
      mex::Mat<uint8_t> is_ordered(1, corners.size());
      for (int i = 0; i < int(corners.size()); ++i) {
        is_ordered[i] = uint8_t(corners[i].isOrdered);
      }

      plhs[2] = is_ordered.release();
    }

  } else {
    mex::error("unknown command '" + command + "'");
  }
}
