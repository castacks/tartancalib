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
#ifndef DELTILLE_APPS_TARGET_DETECTOR_H
#define DELTILLE_APPS_TARGET_DETECTOR_H

//
// simple wrapper demonstrating how to use the deltille library
//

#include <vector>
#include <string>
#include <stdexcept>

#include <opencv2/core.hpp>

#include <deltille/GridDetectorContext.h>
#include <deltille/TaggedBoardIndexer.h>


/**
 */
struct CalibrationCorner {
  int boardId;   //< the board id
  int pointId;   //< point index
  int isOrdered; //< 1 if ordered

  double x; //< x-coordinate (column)
  double y; //< y-coordinate (row)

  CalibrationCorner(double x_ = -1.0, double y_ = -1.0, int b_id = -1,
                    int p_id = -1, int o = 0)
      : boardId(b_id), pointId(p_id), isOrdered(o), x(x_), y(y_) {}

  bool isValid() const { return CalibrationCorner::IsValid(*this); }

  static inline bool IsValid(const CalibrationCorner &c) {
    return c.x >= 0.0 && c.y >= 0.0;
  }

private:
  double _padding = 0.0;

  friend std::ostream &operator<<(std::ostream &os,
                                  const CalibrationCorner &c) {
    os << c.boardId << "," << c.pointId << "," << c.isOrdered << "," << c.x
       << "," << c.y;
    return os;
  }
};

template <class SaddlePointType, class FloatImageType = double,
          typename ImageType = uint8_t>
static inline int
FindBoardsHelper(const cv::Mat &image, const cv::Size &board_size,
                 std::vector<orp::calibration::BoardObservation> &boards) {
  boards.clear();
  orp::calibration::GridDetectorContext<SaddlePointType, ImageType,
                                        FloatImageType>gd(image);

  return gd.findBoards(board_size, boards, false);
}

/**
 * \return the number of detected corners
 */
template <class SaddlePointType, class FloatImageType = double>
static inline int
FindBoards(const cv::Mat &I, const cv::Size &board_size,
           std::vector<orp::calibration::BoardObservation> &boards) {
  switch (I.type()) {
  case cv::DataType<uint8_t>::type:
    return FindBoardsHelper<SaddlePointType, FloatImageType, uint8_t>(
        I, board_size, boards);
  case cv::DataType<uint16_t>::type:
    return FindBoardsHelper<SaddlePointType, FloatImageType, uint16_t>(
        I, board_size, boards);
  case cv::DataType<uint32_t>::type:
    return FindBoardsHelper<SaddlePointType, FloatImageType, uint32_t>(
        I, board_size, boards);
  case cv::DataType<float>::type:
    return FindBoardsHelper<SaddlePointType, FloatImageType, float>(
        I, board_size, boards);
  case cv::DataType<double>::type:
    return FindBoardsHelper<SaddlePointType, FloatImageType, double>(
        I, board_size, boards);
  default:
    throw std::runtime_error("unsupported image type/depth");
  }
}

class TargetDetector {
  using CornerVector = std::vector<CalibrationCorner>;

public:
  /**
   */
  TargetDetector(const std::string &target_dsc_fn, int board_width = 22,
                 int board_height = 22)
      : _indexer() {
    if (!target_dsc_fn.empty()) {
      _indexer.addBoardDefinitions(target_dsc_fn);
    } else {
      _indexer.chessboard_col = board_width;
      _indexer.chessboard_row = board_height;
    }

    if (_indexer.board_defs.empty()) {
      throw std::invalid_argument("invalid target *.dsc file '" + target_dsc_fn +
                             "'" + " no boards were found");
    }
  }

  /**
   */
  void run(const cv::Mat &src, CornerVector &corners,
           cv::Mat *debug_image = nullptr) {
    using namespace orp::calibration;

    if (debug_image) {
      cv::cvtColor(src, *debug_image, cv::COLOR_GRAY2BGR);
      _indexer.setDebugImage(*debug_image);
    }

    const cv::Size board_size(_indexer.chessboard_col, _indexer.chessboard_row);

    std::vector<orp::calibration::BoardObservation> boards;
    if (is_triangular()) {
      FindBoards<MonkeySaddlePointSpherical>(src, board_size, boards);
      _indexer.fixTriangleBoards(src, boards);
    } else {
      FindBoards<SaddlePoint>(src, board_size, boards);
      _indexer.fixCheckerBoards(src, boards);
    }

    corners.clear();
    corners.reserve(boards.size() * board_size.area());
    for (auto &b : boards) {
      for (size_t i = 0; i < b.corner_locations.size(); ++i) {
        const auto &c = b.corner_locations[i];
        corners.emplace_back(c.x, c.y, b.board_id, int(i), b.indexed);
      }
    }

    if (debug_image) {
      auto thickness = std::max(
          1, int(0.5 +
                 sqrt(debug_image->size().area()) /
                     orp::calibration::detector_params.working_resolution));

      for (const auto &b : boards) {
        auto is_triangular = _indexer.hasDefinitions() && b.indexed &&
                             _indexer.board_defs[b.board_id].triangular;
        orp::calibration::drawCheckerboardCorners(*debug_image, b, thickness,
                                                  is_triangular);
      }
    }
  }

private:
  orp::calibration::TaggedBoardIndexer _indexer;

  inline bool is_triangular() const {
    return !_indexer.board_defs.empty() &&
           _indexer.board_defs.front().triangular;
  }
};



#endif
