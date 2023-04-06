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
 * DetectorTools.h
 *
 *  Created on: Oct 12, 2016
 *      Author: mperdoch
 */
#pragma once
#ifndef INCLUDE_CHECKERBOARD_DETECTOR_DETECTORTOOLS_H_
#define INCLUDE_CHECKERBOARD_DETECTOR_DETECTORTOOLS_H_

#include <algorithm>
#include <opencv2/opencv.hpp>
#include <vector>

#include <deltille/TagFamily.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

// #define DEBUG_INDEXING 1            // levels 1..6
// #define DEBUG_REINDEXING 1          // levels 1..3
// #define DEBUG_TIMING

namespace orp {
namespace calibration {
struct BoardObservation {
  // seen in the image
  cv::Mat board;
  std::vector<cv::Point2f> corner_locations;

  // board index in the board definitions
  int board_id;

  // are corners ordered in any way?
  bool indexed;

  // detected april tags
  std::vector<TagDetection> tags;

  static bool orderById(const BoardObservation &a, const BoardObservation &b);
};

struct SaddleClusterDesc {
  typedef struct SaddlePoint PointType;

  int id;
  std::vector<int> idxs;
  double cx, cy;
  double a1, a2;
  SaddleClusterDesc() : id(-1), cx(0.0), cy(0.0) {}

  static bool sortByClusterSize(const SaddleClusterDesc &a,
                                const SaddleClusterDesc &b) {
    return a.idxs.size() < b.idxs.size();
  }

  static bool clusterUnused(const SaddleClusterDesc &c) { return c.id < 0; }

  static bool isTriangular() { return false; }

  template <typename PointType>
  void computeClusterMeans(const std::vector<PointType> &pts) {
    int sigp = 0, sigm = 0;
    std::vector<int> sgn(idxs.size());
    for (size_t j = 0; j < idxs.size(); j++) {
      double cost = (cos(pts[idxs[j]].a1) * cos(pts[idxs[0]].a1) +
                     sin(pts[idxs[j]].a1) * sin(pts[idxs[0]].a1));
      if (cost < 0) {
        sgn[j] = -1;
        sigm++;
      } else {
        sgn[j] = 1;
        sigp++;
      }
    }
    double sign = sigp < sigm ? -1 : 1;
    // reset means
    a1 = 0, a2 = 0;
    cx = 0;
    cy = 0;
    double c1 = 0, c2 = 0;

    for (size_t j = 0; j < idxs.size(); j++) {
      const PointType &pt = pts[idxs[j]];
      cx += pt.x;
      cy += pt.y;

      a1 += sign * sgn[j] * sin(pt.a1);
      c1 += sign * sgn[j] * cos(pt.a1);

      a2 += sin(pt.a2);
      c2 += cos(pt.a2);
    }
    a1 = atan(a1 / c1);
    a2 = atan(a2 / c2);

    cx /= idxs.size();
    cy /= idxs.size();
  }
};

struct MonkeySaddleClusterDesc {
  typedef struct MonkeySaddlePoint PointType;
  int id;
  std::vector<int> idxs;
  double cx, cy;
  std::complex<double> a1, a2, a3;

  MonkeySaddleClusterDesc() : id(-1), cx(0.0), cy(0.0) {}

  static bool sortByClusterSize(const MonkeySaddleClusterDesc &a,
                                const MonkeySaddleClusterDesc &b) {
    return a.idxs.size() < b.idxs.size();
  }

  static bool clusterUnused(const MonkeySaddleClusterDesc &c) {
    return c.id < 0;
  }

  static bool isTriangular() { return true; }

  template <typename PointType>
  void computeClusterMeans(const std::vector<PointType> &pts) {
    std::complex<double> c1 = 0, c2 = 0, c3 = 0;
    a1 = 0;
    a2 = 0;
    a3 = 0;
    cx = 0;
    cy = 0;
    for (size_t j = 0; j < idxs.size(); j++) {
      const PointType &pt = pts[idxs[j]];
      cx += pt.x;
      cy += pt.y;

      a1 += sin(pt.a1);
      c1 += cos(pt.a1);

      a2 += sin(pt.a2);
      c2 += cos(pt.a2);

      a3 += sin(pt.a3);
      c3 += cos(pt.a3);
    }
    a1 = atan(a1 / c1);
    a2 = atan(a2 / c2);
    a3 = atan(a3 / c3);

    cx /= idxs.size();
    cy /= idxs.size();
  }
};

struct Quad {
  int i00_, i01_, i10_, i11_;
  Quad(int i00, int i01, int i10, int i11)
      : i00_(i00), i01_(i01), i10_(i10), i11_(i11) {}
};

struct CoordIndex {
  int r, c;
  int id;
  CoordIndex(int r_, int c_, int id) : r(r_), c(c_), id(id) {}
};

void stretchIntensities(cv::InputArray input, cv::OutputArray output);

template <typename FloatPoint>
void subpixel_marker(cv::Mat &img, const FloatPoint &pt, float sz,
                     const cv::Scalar &color, int thickness = 1,
                     int lineType = CV_AA) {
  const int shift = 16;
  const float ss = (1 << shift);
  cv::rectangle(img, cv::Point(int((pt.x - sz) * ss), int((pt.y - sz) * ss)),
                cv::Point(int((pt.x + sz) * ss), int((pt.y + sz) * ss)), color,
                thickness, lineType, shift);
}

template <typename FloatPoint>
void subpixel_line(cv::Mat &img, const FloatPoint &pt1, const FloatPoint &pt2,
                   const cv::Scalar &color, int thickness = 1,
                   int lineType = CV_AA) {
  const int shift = 16;
  const float ss = (1 << shift);
  cv::line(img, cv::Point(int(pt1.x * ss), int(pt1.y * ss)),
           cv::Point(int(pt2.x * ss), int(pt2.y * ss)), color, thickness,
           lineType, shift);
}

void drawCheckerboardCorners(cv::Mat &result, const BoardObservation &obs,
                             int thickness = 1, bool triangular_board = false);
cv::Mat drawCheckerboardCorners(const cv::Mat &result,
                                const BoardObservation &obs, int thickness = 1,
                                bool triangular_board = false);
void drawCheckerboardCornersOnly(cv::Mat &result, const BoardObservation &obs,
                                 int thickness, bool triangular_board = false);

double sqr(double a);
double mod(double a, double b);
double mod2pi(double vin);
float mod2pi(float vin);
void hessianResponse(const cv::Mat &inputImage, cv::Mat &outputImage);
void adjustGamma(cv::InputArray input, cv::OutputArray output,
                 double gamma = 1.0);
void stretchIntensities(cv::InputArray input, cv::OutputArray output);
int findRoot(std::vector<int> &cluster_ids, int id);

bool point_comparator(const cv::Point2i &a, const cv::Point2i &b);

//#ifdef _WIN32
#ifndef __linux__
template <class ForwardIterator, class T>
void iota(ForwardIterator first, ForwardIterator last, T value) {
  while (first != last) {
    *first++ = value;
    ++value;
  }
}
#endif

template <typename PointType>
double distance2(const PointType &a, const PointType &b) {
  const double x = a.x - b.x, y = a.y - b.y;
  return x * x + y * y;
}

template <typename Float>
inline void solveCubicPolynomial(const Float *a, std::complex<Float> roots[3]) {
  // compute roots of cubic equation:
  // a4 * x^3 + a3 * x ^2 + a2 * x + a1 = 0
  //
  // syms a4 a3 a2 a1
  // roots([a4,a3,a2,a1])
  //
  // (((a1/(2*a4) + a3^3/(27*a4^3) - (a2*a3)/(6*a4^2))^2 + (- a3^2/(9*a4^2) +
  // a2/(3*a4))^3)^(1/2) - a3^3/(27*a4^3) - a1/(2*a4) + (a2*a3)/(6*a4^2))^(1/3)
  // - a3/(3*a4) - (- a3^2/(9*a4^2) + a2/(3*a4))/(((a1/(2*a4) + a3^3/(27*a4^3) -
  // (a2*a3)/(6*a4^2))^2 + (a2/(3*a4) - a3^2/(9*a4^2))^3)^(1/2) - a3^3/(27*a4^3)
  // - a1/(2*a4) + (a2*a3)/(6*a4^2))^(1/3)
  // (- a3^2/(9*a4^2) + a2/(3*a4))/(2*(((a1/(2*a4) + a3^3/(27*a4^3) -
  // (a2*a3)/(6*a4^2))^2 + (a2/(3*a4) - a3^2/(9*a4^2))^3)^(1/2) - a3^3/(27*a4^3)
  // - a1/(2*a4) + (a2*a3)/(6*a4^2))^(1/3)) - (3^(1/2)*((- a3^2/(9*a4^2) +
  // a2/(3*a4))/(((a1/(2*a4) + a3^3/(27*a4^3) - (a2*a3)/(6*a4^2))^2 + (a2/(3*a4)
  // - a3^2/(9*a4^2))^3)^(1/2) - a3^3/(27*a4^3) - a1/(2*a4) +
  // (a2*a3)/(6*a4^2))^(1/3) + (((a1/(2*a4) + a3^3/(27*a4^3) -
  // (a2*a3)/(6*a4^2))^2 + (- a3^2/(9*a4^2) + a2/(3*a4))^3)^(1/2) -
  // a3^3/(27*a4^3) - a1/(2*a4) + (a2*a3)/(6*a4^2))^(1/3))*1i)/2 - a3/(3*a4) -
  // (((a1/(2*a4) + a3^3/(27*a4^3) - (a2*a3)/(6*a4^2))^2 + (- a3^2/(9*a4^2) +
  // a2/(3*a4))^3)^(1/2) - a3^3/(27*a4^3) - a1/(2*a4) +
  // (a2*a3)/(6*a4^2))^(1/3)/2
  // (- a3^2/(9*a4^2) + a2/(3*a4))/(2*(((a1/(2*a4) + a3^3/(27*a4^3) -
  // (a2*a3)/(6*a4^2))^2 + (a2/(3*a4) - a3^2/(9*a4^2))^3)^(1/2) - a3^3/(27*a4^3)
  // - a1/(2*a4) + (a2*a3)/(6*a4^2))^(1/3)) + (3^(1/2)*((- a3^2/(9*a4^2) +
  // a2/(3*a4))/(((a1/(2*a4) + a3^3/(27*a4^3) - (a2*a3)/(6*a4^2))^2 + (a2/(3*a4)
  // - a3^2/(9*a4^2))^3)^(1/2) - a3^3/(27*a4^3) - a1/(2*a4) +
  // (a2*a3)/(6*a4^2))^(1/3) + (((a1/(2*a4) + a3^3/(27*a4^3) -
  // (a2*a3)/(6*a4^2))^2 + (- a3^2/(9*a4^2) + a2/(3*a4))^3)^(1/2) -
  // a3^3/(27*a4^3) - a1/(2*a4) + (a2*a3)/(6*a4^2))^(1/3))*1i)/2 - a3/(3*a4) -
  // (((a1/(2*a4) + a3^3/(27*a4^3) - (a2*a3)/(6*a4^2))^2 + (- a3^2/(9*a4^2) +
  // a2/(3*a4))^3)^(1/2) - a3^3/(27*a4^3) - a1/(2*a4) +
  // (a2*a3)/(6*a4^2))^(1/3)/2

  // p = - a3 / (3 * a4)
  // r =   a2 / (3 * a4)
  // p2 = p * p      //  a3^2/(9 *a4^2)
  // p3 = p * p * p  // -a3^3/(27*a4^3)
  //
  // q = p3 + (a2 * a3 - 3.0 * a4 * a1) / (6.0 * a4 * a4)     =               p3
  // + a2 * a3 / (6.0 * a4 * a4) - a1 / (2.0 * a4)
  // q2 = q * q                                               =              (p3
  // + a2 * a3 / (6.0 * a4 * a4) - a1 / (2.0 * a4) )^2
  //
  // x      =  p -
  // x2,3   =  p + {q + [q2 + (r-p2)^3]^1/2}^1/3  + {q - [q2 +
  // (r-p2)^3]^1/2}^1/3

  //
  // (((a1/(2*a4) + a3^3/(27*a4^3) - (a2*a3)/(6*a4^2))^2 + (a2/(3*a4) -
  // a3^2/(9*a4^2))^3)^(1/2) - (a1/(2*a4) + a3^3/(27*a4^3) -
  // (a2*a3)/(6*a4^2)))^(1/3))
  // ((---------------------R-------------------------^2 +
  // ------------Q--------------^3)^(1/2) - ----------------- R
  // ---------------------------)^(1/3)
  // (----------------------------------------------D-----------------------------------)^(1/2)
  // -                   R                            )^(1/3)

  // convenience mapping
  const Float &a1 = a[0], &a2 = a[1], &a3 = a[2], &a4 = a[3];

  Float p = -a3 / (3 * a4);
  Float r = a2 / (3 * a4);
  Float R = p * p * p + (a2 * a3 - 3.0 * a4 * a1) / (6.0 * a4 * a4);
  Float R2 = R * R;
  Float Q = r - p * p;
  Float Q3 = Q * Q * Q;
  Float D = R2 + Q3;

  std::complex<Float> S, T;
  if (D < 0) {
    // conjugate complex roots
    S = pow(R + sqrt(std::complex<Float>(D)), 1.0 / 3.0);
    T = pow(R - sqrt(std::complex<Float>(D)), 1.0 / 3.0);
  } else if (D >= 0) {
    // real or multiple roots
    S = cbrt(R + sqrt(D));
    T = cbrt(R - sqrt(D));
  }
  roots[0] = S + T + p;
  roots[1] =
      -(S + T) / 2.0 + p + std::complex<double>(0, sqrt(3.0)) * (S - T) / 2.0;
  roots[2] =
      -(S + T) / 2.0 + p - std::complex<double>(0, sqrt(3.0)) * (S - T) / 2.0;
}

template <typename ElementType, typename IndexType>
void partial_sort_indexes(const ElementType *v, IndexType *idx, size_t len,
                          size_t up_to) {
// initialize original index locations
#ifndef __linux__
  iota(idx, idx + len, 0);
#else
  std::iota(idx, idx + len, 0);
#endif

  // sort indexes based on comparing values in v
  std::partial_sort(idx, idx + up_to, idx + len,
                    [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
}

template <typename PointType> bool PointIsInf(const PointType &pt) {
  return (std::isinf(pt.x) || std::isinf(pt.y));
}

template <typename PointType> bool PointIsValid(const PointType &pt) {
  return (pt.x != -1 || pt.y != -1);
}

template <typename PointType, typename ClusterDesc>
void clusterPoints2(const std::vector<PointType> &pts,
                    const cv::Size &input_size, std::vector<int> &cluster_ids,
                    std::vector<ClusterDesc> &cluster_stats, int &num_clusters,
                    double threshold = 2.0) {
  // printf("Clustering %zd points\n", pts.size());
  cluster_ids.resize(pts.size());
  // check if threshold is higher than pixel grid resolution
  CV_Assert(threshold >= 1.0);
  double threshold2 = threshold * threshold;
  num_clusters = 0;

  // cluster links
  cv::Mat cluster(input_size, CV_32S);
  cluster.setTo(cv::Scalar(-1));

  cluster_stats.clear();
  cluster_stats.resize(pts.size());

  for (size_t i = 0; i < pts.size(); i++) {
    const PointType &pt = pts[i];
    int &loc_cluster = cluster.at<int>(int(round(pt.y)), int(round(pt.x)));
    if (loc_cluster == -1) {
      // a new cluster at this location (with point i being it's root)
      loc_cluster = int(i);
      cluster_stats[i].id = int(i);
      num_clusters++;
    }
    // assign point to cluster
    cluster_ids[i] = loc_cluster;
    // push point to cluster struct
    ClusterDesc &c = cluster_stats[loc_cluster];
    c.idxs.push_back(int(i));
  }

  int num_dist = 0;
  for (size_t i = 0; i < pts.size(); i++) {
    const PointType &pt1 = pts[i];
    // find distances to already clustered points
    int cluster_id = findRoot(cluster_ids, int(i));
    ClusterDesc &cur = cluster_stats[cluster_id];

    // const int lbx = std::max(0, int(pt1.x - threshold)), lby = std::max(0,
    // int(pt1.y - threshold)),
    //          ubx = std::min(input_size.width - 1, int(pt1.x + threshold +
    //          .5)), uby = std::min(input_size.height - 1, int(pt1.y +
    //          threshold + .5));
    const int lbx = std::max(0, int(pt1.x)), lby = std::max(0, int(pt1.y)),
              ubx = std::min(input_size.width - 1, int(pt1.x + .5)),
              uby = std::min(input_size.height - 1, int(pt1.y + .5));

    for (int y = lby; y <= uby; y++) {
      for (int x = lbx; x <= ubx; x++) {
        // go through relevant neighbouring clusters and check if they need to
        // be merged in based on the distance to this point
        int other_id = findRoot(cluster_ids, cluster.at<int>(y, x));
        if (other_id == -1)
          continue;
        if (other_id == cluster_id)
          continue;

        ClusterDesc &other = cluster_stats[other_id];
        for (size_t idx = 0; idx < other.idxs.size(); idx++) {
          const PointType &pt2 = pts[other.idxs[idx]];
          double dist = distance2(pt1, pt2);
          num_dist++;
          // elongated cluster? bail out soon
          // if (dist > 5*threshold) break;
          if (dist < threshold2) {
            // merge cluster other_id and cluster_id...
            cluster_ids[other.id] = cluster_id;
            other.id = -1;
            for (size_t k = 0; k < other.idxs.size(); k++) {
              cur.idxs.push_back(other.idxs[k]);
              cluster_ids[other.idxs[k]] = cluster_id; // bug fixed by Hyowon !
            }
            other.idxs.clear();
            num_clusters--;
          }
        }
      }
    }
  }
  // printf("Remains %d clusters\n", num_clusters);
  // renumber clusters
  num_clusters = 0;
  for (size_t i = 0; i < cluster_stats.size(); i++)
    if (cluster_stats[i].id != -1)
      cluster_stats[i].id = num_clusters++;
  // collapse trees and fix cluster ids to roots
  for (size_t i = 0; i < pts.size(); i++) {
    findRoot(cluster_ids, int(i));
    // int c = findRoot(cluster_ids, int(i));
    // if (c != cluster_ids[i])
    //  ;
    // printf("cluster root: %d, id: %d\n", c, cluster_ids[i]);
  }
  // renumber to match cluster_stats
  for (size_t i = 0; i < pts.size(); i++) {
    // if (cluster_stats[cluster_ids[i]].id < 0 ||
    //    cluster_stats[cluster_ids[i]].id > (int) pts.size())
    //  ;
    // printf("invalid cluster id, cluster %d, stats id: %d\n", cluster_ids[i],
    // cluster_stats[cluster_ids[i]].id);

    cluster_ids[i] = cluster_stats[cluster_ids[i]].id;
  }
  // squeeze cluster stats
  typename std::vector<ClusterDesc>::iterator last = std::remove_if(
      cluster_stats.begin(), cluster_stats.end(), ClusterDesc::clusterUnused);
  cluster_stats.resize(last - cluster_stats.begin());
  // printf("Total clusters: %d, distances: %d\n", num_clusters, num_dist);

  // combine points in clusters...
  for (size_t i = 0; i < cluster_stats.size(); i++)
    cluster_stats[i].computeClusterMeans(pts);
}

#if DEBUG_INDEXING
extern cv::Mat DEBUG;
#endif
}
}

#endif /* INCLUDE_CHECKERBOARD_DETECTOR_DETECTORTOOLS_H_ */
