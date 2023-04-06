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
#include <deltille/PolynomialFit.h>
#include <deltille/PolynomialSaddleDetectorContext.h>

#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>

using namespace orp::calibration;

template <class SaddlePointType, typename InputImageType = uint8_t,
          typename FloatImageType = double>
inline void refineCorners(const cv::Mat &I,
                          const std::vector<cv::Point2f> &x_init,
                          std::vector<SaddlePointType> &x_refined,
                          int half_kernel_size, int max_iterations) {
  PolynomialFit<SaddlePointType> pfit;
  pfit.initSaddleFitting(half_kernel_size);

  cv::Mat I_s;
  cv::filter2D(I, I_s, cv::DataType<FloatImageType>::depth,
               pfit.getSmoothingKernel());

  pfit.saddleSubpixelRefinement(I_s, x_init, x_refined, max_iterations, false);
}

template <class SaddlePointType>
mex::Mat<double> toMatlab(const std::vector<SaddlePointType> &x) {
  mex::Mat<double> ret(2, x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret(0, int(i)) = x[i].x;
    ret(1, int(i)) = x[i].y;
  }

  return ret;
}

template <class SaddlePointType, typename InputImageType = uint8_t,
          typename FloatImageType = double>
inline std::vector<SaddlePointType>
run(const cv::Mat &I, const std::vector<cv::Point2f> &x_init,
    int half_kernel_size, int max_iterations) {
  std::vector<SaddlePointType> ret(x_init.size());
  refineCorners<SaddlePointType, InputImageType, FloatImageType>(
      I, x_init, ret, half_kernel_size, max_iterations);

  return ret;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  static const char *USAGE =
      "x = fn('type', I, x_init, kernel_size, num_iters = 5)";
  mex::nargchk(4, 5, nrhs, USAGE);

  const auto pattern_type = mex::getString(prhs[0]);
  const mex::Mat<double> image(prhs[1]);
  const mex::Mat<double> x(prhs[2]);
  int half_kernel_size = mex::getNumber<int>(prhs[3]);
  int num_iters = nrhs > 4 ? mex::getNumber<int>(prhs[4]) : 20;

  const auto I = mex::mex2cv(image);
  std::vector<cv::Point2f> x_init(x.cols());
  for (size_t i = 0; i < x_init.size(); ++i) {
    x_init[i].x = float(x(0, int(i)));
    x_init[i].y = float(x(1, int(i)));
  }

  if (pattern_type == "square") {
    auto x_refined =
        run<SaddlePoint, double>(I, x_init, half_kernel_size, num_iters);
    plhs[0] = toMatlab(x_refined).release();
  } else {
    auto x_refined = run<MonkeySaddlePointSpherical, double>(
        I, x_init, half_kernel_size, num_iters);
    plhs[0] = toMatlab(x_refined).release();
  }
}
