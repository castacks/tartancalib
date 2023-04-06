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
#ifndef INCLUDE_DELTILLE_TAGGED_BOARD_INDEXDER_H
#define INCLUDE_DELTILLE_TAGGED_BOARD_INDEXDER_H


#include <opencv2/opencv.hpp>
#include <vector>

#include <deltille/DetectorTools.h>

#include <fstream>
#include <iostream>
#include <map>
#include <memory>

namespace orp {
namespace calibration {
typedef std::map<std::string, std::shared_ptr<TagFamily>> AprilTagFamilies;

struct BoardDefinition {
  typedef std::map<int, cv::Point2i> BoardTagMap;
  typedef BoardTagMap::iterator BoardTagMapIterator;

  int id;
  int target_id;
  std::string family;
  float border_bits;
  int cols, rows;
  double size;
  bool triangular;
  BoardTagMap tag_locations; // code -> location of lower left corner on the
                             // checkerboard
  cv::Mat corner_locations;  // matrix of Vec3d -> 3D locations of all corners
  std::shared_ptr<TagFamily> detector;
};

const AprilTags::TagCodes &tagFamilyNameToCodes(const std::string &family);

void readBoardDefinitions(std::istream &in, std::vector<BoardDefinition> &defs,
                          AprilTagFamilies &detectors);
bool updateBoardsWithCalibratedTargetFile(std::istream &in,
                                          std::vector<BoardDefinition> &defs);

void writeBoardObservations(const char *filename,
                            const std::vector<BoardDefinition> &defs,
                            const std::vector<BoardObservation> &boards);

// interface with Hyowon's calibration
void writeBoardDefinitionsHH(const char *filename,
                             const std::vector<BoardDefinition> &defs);
void writeBoardObservationsHH(const char *filename,
                              const std::vector<BoardDefinition> &defs,
                              const std::vector<BoardObservation> &boards,
                              bool output_unindexed_boards = false);

template <typename Point3d, typename Point2d>
void get2Dto3DCorrespondences(BoardDefinition &def, BoardObservation &board,
                              std::vector<Point3d> &objectPoints,
                              std::vector<Point2d> &imagePoints) {
  std::vector<cv::Point2f> &obs_corners = board.corner_locations;
  objectPoints.reserve(obs_corners.size());
  imagePoints.reserve(obs_corners.size());
  for (int j = 0; j < obs_corners.size(); j++) {
    if (obs_corners[j].x != -1 && obs_corners[j].y != -1) {
      int r = j / def.cols, c = j % def.cols;
      cv::Vec3d &pt = def.corner_locations.at<cv::Vec3d>(r, c);
      objectPoints.push_back(Point3d(pt[0], pt[1], pt[2]));
      imagePoints.push_back(Point2d(obs_corners[j].x, obs_corners[j].y));
    }
  }
}

int fixFullCheckerBoardOrientations(const cv::Mat &img,
                                    const cv::Size &board_size,
                                    std::vector<BoardObservation> &boards);
bool fixFullCheckerBoardOrientation(const cv::Mat &img,
                                    const cv::Size &board_size,
                                    BoardObservation &board);

struct TaggedBoardIndexer {
  typedef std::map<int, std::pair<int, cv::Point2i>> TagToBoardMap;

  AprilTagFamilies detectors;
  std::vector<BoardDefinition> board_defs;
  TagToBoardMap tag_to_board_map;
  cv::Mat dbg;
  int chessboard_row;
  int chessboard_col;

  TaggedBoardIndexer() {}

  void init() {
    detectors.clear();
    board_defs.clear();
    tag_to_board_map.clear();
    dbg = cv::Mat();
    chessboard_row = 13;
    chessboard_col = 13;
  }

  bool hasDefinitions() const { return board_defs.size() > 0; }

  void addBoardDefinitions(const std::string &filename) {
    std::ifstream f(filename);
    if (f.good())
      readBoardDefinitions(f, board_defs, detectors);
    updateBoardDefinitions();
  }

  bool readCalibratedTargetFile(const std::string &filename) {
    std::ifstream f(filename);
    if (!f.good())
      return false;

    return updateBoardsWithCalibratedTargetFile(f, board_defs);
  }

  void updateBoardDefinitions() {
    int maxr = 0, maxc = 0;
    for (size_t i = 0; i < board_defs.size(); i++) {
      if (maxr < board_defs[i].rows)
        maxr = board_defs[i].rows;
      if (maxc < board_defs[i].cols)
        maxc = board_defs[i].cols;
    }
    chessboard_row = maxr;
    chessboard_col = maxc;
    tag_to_board_map.clear();

    for (size_t b = 0; b < board_defs.size(); ++b) {
      int board_tag_offset = board_defs[b].detector->tagFamilyOffset;
      for (auto &&tl : board_defs[b].tag_locations)
        tag_to_board_map[tl.first + board_tag_offset] =
            std::make_pair(int(b), tl.second);
    }
  }

  void writeObservationsHH(const std::string &filename,
                           const std::vector<BoardObservation> &boards,
                           bool output_unindexed_boards = false) const {
    writeBoardObservationsHH(filename.c_str(), board_defs, boards,
                             output_unindexed_boards);
  }

  void writedDefinitionsHH(const std::string &filename) const {
    writeBoardDefinitionsHH(filename.c_str(), board_defs);
  }

  void setDebugImage(cv::Mat &dbg_image) { dbg = dbg_image; }

  void fixCheckerBoards(const cv::Mat &img,
                        std::vector<BoardObservation> &boards);
  void fixTriangleBoards(const cv::Mat &img,
                         std::vector<BoardObservation> &boards);
  bool detectDeltaTag(const cv::Mat &img, BoardObservation &obs, int r, int c,
                      bool lower, TagDetection &det);

private:
  std::vector<int> rotation_histo;
  std::vector<int> board_id_histo;
  std::map<int, int> offsetx_histo;
  std::map<int, int> offsety_histo;
};

}; // calibration
}; // orp

#endif
