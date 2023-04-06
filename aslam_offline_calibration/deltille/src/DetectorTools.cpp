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

/*
 * DetectorTools.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: mperdoch
 */

#include <deltille/DetectorParams.h>
#include <deltille/DetectorTools.h>
#include <opencv2/opencv.hpp>

namespace orp {
namespace calibration {

DetectorParams detector_params;

#ifdef DEBUG_INDEXING
cv::Mat DEBUG;
#endif
/*
void test()
{
    cv::Mat input;
    cv::Size board_size;
    std::vector<orp::calibration::BoardObservation> boards;

    orp::calibration::GridDetectorContext<orp::calibration::MonkeySaddlePoint>
grid_detector(input);
    grid_detector.findBoards(board_size, boards, false);

    orp::calibration::GridDetectorContext<orp::calibration::SaddlePoint>
grid_detector2(input);
    grid_detector2.findBoards(board_size, boards, false);

    orp::calibration::GridDetectorContext<orp::calibration::MonkeySaddlePoint,
float> grid_detectorf(input);
    grid_detectorf.findBoards(board_size, boards, false);

    orp::calibration::GridDetectorContext<orp::calibration::SaddlePoint, float>
grid_detector2f(input);
    grid_detector2f.findBoards(board_size, boards, false);
}
*/

template <typename ImageType>
static inline void hessianResponseImpl(const cv::Mat &inputImage,
                                       cv::Mat &outputImage) {
  const int rows = inputImage.rows;
  const int cols = inputImage.cols;
  const int stride = cols;

  // allocate output
  outputImage = cv::Mat::zeros(rows, cols, CV_64FC1);

  // setup input and output pointer to be centered at 1,0 and 1,1 resp.
  auto *in = inputImage.ptr<const ImageType>(1);
  auto *out = outputImage.ptr<double>(1) + 1;
  // double             norm = 1.0 / (255.0*255.0);
  /* move 3x3 window and convolve */
  for (int r = 1; r < rows - 1; ++r) {
    double v11, v12, v21, v22, v31, v32;
    /* fill in shift registers at the beginning of the row */
    v11 = in[-stride];
    v12 = in[1 - stride];
    v21 = in[0];
    v22 = in[1];
    v31 = in[+stride];
    v32 = in[1 + stride];
    /* move input pointer to (1,2) of the 3x3 square */
    in += 2;
    for (int c = 1; c < cols - 1; ++c) {
      /* fetch remaining values (last column) */
      const double v13 = in[-stride];
      const double v23 = *in;
      const double v33 = in[+stride];

      // compute 3x3 Hessian values from symmetric differences.
      double Lxx = (v21 - 2 * v22 + v23);
      double Lyy = (v12 - 2 * v22 + v32);
      double Lxy = (v13 - v11 + v31 - v33) / 4.0f;

      /* normalize and write out */
      *out = Lxx * Lyy - Lxy * Lxy;

      /* move window */
      v11 = v12;
      v12 = v13;
      v21 = v22;
      v22 = v23;
      v31 = v32;
      v32 = v33;

      /* move input/output pointers */
      in++;
      out++;
    }
    out += 2;
  }
}

void hessianResponse(const cv::Mat &inputImage, cv::Mat &outputImage) {
  switch (inputImage.depth()) {
  case cv::DataDepth<uint8_t>::value:
    hessianResponseImpl<uint8_t>(inputImage, outputImage);
    break;
  case cv::DataDepth<uint16_t>::value:
    hessianResponseImpl<uint16_t>(inputImage, outputImage);
    break;
  case cv::DataDepth<float>::value:
    hessianResponseImpl<float>(inputImage, outputImage);
    break;
  case cv::DataDepth<double>::value:
    hessianResponseImpl<double>(inputImage, outputImage);
    break;
  default:
    throw std::runtime_error("unsupported image depth");
  }
}

int findRoot(std::vector<int> &cluster_ids, int id) {
  if (id == -1)
    return id;
  // check root
  int root = cluster_ids[id];
  // if it is not me, and tree is at least 2 hops deep
  if (root != id && cluster_ids[root] != root) {
    // find root first...
    do {
      root = cluster_ids[root];
    } while (cluster_ids[root] != root);
    // full flatten tree
    while (cluster_ids[id] != root) {
      int temp = cluster_ids[id];
      cluster_ids[id] = root;
      id = temp;
    }
  }
  return root;
}

double sqr(double a) { return a * a; }

double mod(double a, double b) { return a - std::floor(a / b) * b; }

double mod2pi(double vin) {
  // MOD2PI returns a result in [-Pi, Pi]
  double twopi = 2.0 * M_PI;
  double twopi_inv = 1.0 / twopi;
  double absv = fabs(vin);
  double q = absv * twopi_inv + 0.5;
  double r = absv - floor(q) * twopi;
  if (vin < 0)
    return -r;
  else
    return r;
}

bool BoardObservation::orderById(const BoardObservation &a,
                                 const BoardObservation &b) {
  return a.board_id < b.board_id;
}

//! Returns a result in [-Pi, Pi]
float mod2pi(float vin) {
  const float twopi = 2 * (float)M_PI;
  const float twopi_inv = 1.f / (2.f * (float)M_PI);
  float absv = std::abs(vin);
  float q = absv * twopi_inv + 0.5f;
  int qi = (int)q;
  float r = absv - qi * twopi;
  return (vin < 0) ? -r : r;
}

void adjustGamma(cv::InputArray input, cv::OutputArray output, double gamma) {
  // build a lookup table mapping the pixel values [0, 255] to
  // their adjusted gamma values
  double invGamma = 1.0 / gamma;
  cv::Mat table(1, 256, CV_8U);
  for (int i = 0; i < 256; i++)
    table.at<uint8_t>(0, i) =
        cv::saturate_cast<uint8_t>(pow(double(i) / 255.0, invGamma) * 255);

  // apply gamma correction using the lookup table
  cv::LUT(input, table, output);
}

static void rescaleGrayLevelMat(const cv::Mat &inputMat, cv::Mat &outputMat,
                                const float histogramClippingLimit) {
  cv::Mat small;
  cv::resize(inputMat, small, cv::Size(), 0.125, 0.125, cv::INTER_LINEAR);
  int type = CV_MAT_TYPE(inputMat.type());
  if (type == CV_32F) {
    float *begin = small.ptr<float>();
    float *end = small.ptr<float>() + small.rows * small.cols;
    std::sort(begin, end);
    cv::normalize(inputMat, outputMat,
                  *(begin + int((end - begin) * histogramClippingLimit)), 255.f,
                  cv::NORM_MINMAX);
  }
}

void stretchIntensities(cv::InputArray input, cv::OutputArray output) {
  // build a lookup table mapping the pixel values [0, 255] to
  // their adjusted gamma values
  if (CV_MAT_TYPE(input.type()) == CV_8U ||
      CV_MAT_TYPE(input.type()) == CV_16U) {
    double minVal, maxVal;
    cv::minMaxIdx(input, &minVal, &maxVal);
    input.getMat().convertTo(output, input.type(), 255.0 / (maxVal - minVal),
                             -minVal);
  } else if (CV_MAT_TYPE(input.type()) == CV_32F) {
    // get histogram density probability in order to cut values under above
    // edges limits (here 5-95%)... usefull for HDR pixel errors cancellation
    cv::Mat out;
    rescaleGrayLevelMat(input.getMat(), out, 0.05f);
    out.convertTo(output, CV_8U, 1, 0);
  }
}

bool point_comparator(const cv::Point2i &a, const cv::Point2i &b) {
  if (a.y != b.y)
    return a.y < b.y;
  else
    return (a.x < b.x);
}

cv::Mat drawCheckerboardCorners(const cv::Mat &result,
                                const BoardObservation &obs, int thickness,
                                bool triangular_board) {
  cv::Mat output = result;
  drawCheckerboardCorners(output, obs, thickness, triangular_board);
  return output;
}

void drawCheckerboardCornersOnly(cv::Mat &result, const BoardObservation &obs,
                                 int thickness, bool triangular_board) {
  const int line_max = 10;
  static const cv::Scalar line_colors[line_max] = {
      //            R, G, B
      cv::Scalar(0, 0, 255),  cv::Scalar(0, 128, 255), cv::Scalar(0, 224, 224),
      cv::Scalar(0, 255, 0),  cv::Scalar(224, 224, 0), cv::Scalar(255, 128, 0),
      cv::Scalar(255, 0, 0),  cv::Scalar(255, 0, 128), cv::Scalar(224, 0, 224),
      cv::Scalar(128, 0, 255)};
  float rs = 1.0f + thickness;
  const cv::Mat &checkerboard = obs.board;
  const std::vector<cv::Point2f> &corners = obs.corner_locations;
  for (int r = 0; r < checkerboard.rows; ++r) {
    cv::Scalar color = line_colors[r % line_max];
    for (int c = 0; c < checkerboard.cols; ++c) {
      int id1, id2;
      bool both = true;
      if (c + 1 < checkerboard.cols) {
        id1 = checkerboard.at<int>(r, c);
        id2 = checkerboard.at<int>(r, c + 1);
        if (id1 != -1 && id2 != -1) {
          subpixel_marker(result, corners[id1], rs, color, thickness);
          subpixel_marker(result, corners[id2], rs, color, thickness);
        }
      } else
        both = false;
      if (r + 1 < checkerboard.rows) {
        id1 = checkerboard.at<int>(r, c);
        id2 = checkerboard.at<int>(r + 1, c);
        if (id1 != -1 && id2 != -1) {
          subpixel_marker(result, corners[id1], rs, color, thickness);
          subpixel_marker(result, corners[id2], rs, color, thickness);
        }
      } else
        both = false;
      if (both && triangular_board) {
        id1 = checkerboard.at<int>(r, c + 1);
        id2 = checkerboard.at<int>(r + 1, c);
#if 0
                if (id1 != -1 && id2 != -1)
                {
                    const cv::Point2f &s1 = corners[id1], &s2 = corners[id2];
                   // cv::line(result, cv::Point(s1.x * ss, s1.y * ss), cv::Point(s2.x * ss, s2.y * ss), color, thickness, CV_AA, shift);
                }
#endif
      }
    }
  }
}

void drawCheckerboardCorners(cv::Mat &result, const BoardObservation &obs,
                             int thickness, bool triangular_board) {
  const int line_max = 10;
  static const cv::Scalar line_colors[line_max] = {
      //            R, G, B
      cv::Scalar(0, 0, 255),  cv::Scalar(0, 128, 255), cv::Scalar(0, 224, 224),
      cv::Scalar(0, 255, 0),  cv::Scalar(224, 224, 0), cv::Scalar(255, 128, 0),
      cv::Scalar(255, 0, 0),  cv::Scalar(255, 0, 128), cv::Scalar(224, 0, 224),
      cv::Scalar(128, 0, 255)};
  float rs = 1.0f + thickness;
  const cv::Mat &checkerboard = obs.board;
  const std::vector<cv::Point2f> &corners = obs.corner_locations;
  //    printf("Drawing board %d:\n", obs.board_id);
  //    std::cout << checkerboard << std::endl;
  bool index_shown = false;
  for (int r = 0; r < checkerboard.rows; ++r) {
    cv::Scalar color = line_colors[r % line_max];
    if (!obs.indexed) {
      color[0] = std::max(128, int(color[0]));
      color[1] = std::max(128, int(color[1]));
      color[2] = std::max(128, int(color[2]));
    }
    for (int c = 0; c < checkerboard.cols; ++c) {
      int id1, id2;
      bool both = true;
      if (c + 1 < checkerboard.cols) {
        id1 = checkerboard.at<int>(r, c);
        id2 = checkerboard.at<int>(r, c + 1);
        if (id1 != -1 && id2 != -1) {
          subpixel_line(result, corners[id1], corners[id2], color, thickness);
          if (obs.indexed) {
            subpixel_marker(result, corners[id1], rs, color, thickness);
            subpixel_marker(result, corners[id2], rs, color, thickness);
            if (!index_shown) {
              cv::putText(result, std::to_string(obs.board_id),
                          cv::Point(int(corners[id1].x), int(corners[id1].y)),
                          cv::FONT_HERSHEY_PLAIN, 1.0, cv::Scalar(255, 0, 255));
              index_shown = true;
            }
          }
        }
      } else
        both = false;
      if (r + 1 < checkerboard.rows) {
        id1 = checkerboard.at<int>(r, c);
        id2 = checkerboard.at<int>(r + 1, c);
        if (id1 != -1 && id2 != -1) {
          subpixel_line(result, corners[id1], corners[id2], color, thickness);
          if (obs.indexed) {
            subpixel_marker(result, corners[id1], rs, color, thickness);
            subpixel_marker(result, corners[id2], rs, color, thickness);
            if (!index_shown) {
              cv::putText(result, std::to_string(obs.board_id),
                          cv::Point(int(corners[id1].x), int(corners[id1].y)),
                          cv::FONT_HERSHEY_PLAIN, 1.0, cv::Scalar(255, 0, 255));
              index_shown = true;
            }
          }
        }
      } else
        both = false;
      if (both && triangular_board) {
        id1 = checkerboard.at<int>(r, c + 1);
        id2 = checkerboard.at<int>(r + 1, c);
        if (id1 != -1 && id2 != -1)
          subpixel_line(result, corners[id1], corners[id2], color, thickness);
      }
    }
  }
  if (!obs.indexed) {
    cv::Scalar color(0, 0, 0);
    for (int r = 0; r < checkerboard.rows; ++r) {
      for (int c = 0; c < checkerboard.cols; ++c) {
        int id1, id2;
        if (c + 1 < checkerboard.cols) {
          id1 = checkerboard.at<int>(r, c);
          id2 = checkerboard.at<int>(r, c + 1);
          if (id1 != -1 && id2 != -1) {
            subpixel_marker(result, corners[id1], rs, color, thickness);
            subpixel_marker(result, corners[id2], rs, color, thickness);
          }
        }
        if (r + 1 < checkerboard.rows) {
          id1 = checkerboard.at<int>(r, c);
          id2 = checkerboard.at<int>(r + 1, c);
          if (id1 != -1 && id2 != -1) {
            subpixel_marker(result, corners[id1], rs, color, thickness);
            subpixel_marker(result, corners[id2], rs, color, thickness);
          }
        }
      }
    }
  }
}
}
}
