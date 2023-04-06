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
 * GridDetector.h
 *
 *  Created on: Dec 1, 2016
 *      Author: mperdoch
 */
#pragma once
#ifndef INCLUDE_CHECKERBOARD_DETECTOR_GRIDDETECTORCONTEXT_H_
#define INCLUDE_CHECKERBOARD_DETECTOR_GRIDDETECTORCONTEXT_H_

#include <unordered_map>

#include <deltille/DetectorTools.h>
#include <deltille/PolynomialFit.h>
#include <deltille/PolynomialSaddleDetectorContext.h>

namespace orp {
namespace calibration {

template <typename SaddlePointType, typename InputImageType = uint8_t,
          typename FloatImageType = float>
struct GridDetectorContext {

  typedef typename std::vector<SaddlePointType> SaddlePointVector;

private:
  // shared parameters of the image
  int width, height;

  // precomputed and shared values for each image, only needs to be computed
  // once
  SaddlePointVector pts;
  int num_pts;
  cv::Mat dist_u, dist_v, vec_u, vec_v, dist;
  cv::Mat NN, issame, isdiff, polarity, input;
  cv::Mat active;

  // indexes of points usable for the initial point in quad selection
  std::vector<int> keypoints;

  // auxiliary parameters
  bool deltilleGrid;
  int max_nn;

  std::unordered_map<uint64_t, bool> edge_cache;

  PolynomialSaddleDetectorContext<SaddlePointType, InputImageType,
                                  FloatImageType>
      detector;

public:
  GridDetectorContext(const cv::Mat &image) : detector(image) {
    input = detector.input_lowres;
    width = input.cols;
    height = input.rows;
    deltilleGrid = SaddlePointType::isTriangular;
  }

  template <typename UnaryPredicateOnIndex>
  std::vector<int> getClosestNNs(int idx, int num,
                                 UnaryPredicateOnIndex fn) const {
    std::vector<int> closenn(num);
    const int *idxs = NN.ptr<int>(idx);
    const uint8_t *mask = active.ptr<const uint8_t>(0);
    int cnt = 0;
    for (int i = 0; i < max_nn; ++i) {
      int nn = idxs[i];
      if (mask[nn] != 0 && nn != idx && fn(nn))
        closenn[cnt++] = nn;
      if (cnt == num)
        break;
    }
    // remove any slack space...
    closenn.resize(cnt);
    return closenn;
  }

  std::vector<int> getActiveIndices() const {
    std::vector<int> idxs(pts.size());
    const uint8_t *mask = active.ptr<const uint8_t>(0);
    int cnt = 0;
    for (size_t i = 0; i < pts.size(); ++i)
      if (mask[i] != 0)
        idxs[cnt++] = int(i);
    // remove any slack space...
    idxs.resize(cnt);
    return idxs;
  }

  bool precomputePolaritiesAndNN(const SaddlePointVector &points) {
    // cleanup any previous state...
    pts.clear();
    edge_cache.clear();
    pts = points;

    // recreate all precomputed data
    num_pts = int(pts.size());
    dist_u.create(num_pts, num_pts, CV_64F);
    dist_v.create(num_pts, num_pts, CV_64F);
    vec_u.create(num_pts, num_pts, CV_64F);
    vec_v.create(num_pts, num_pts, CV_64F);
    dist.create(num_pts, num_pts, CV_64F);
    NN.create(num_pts, num_pts, CV_32S);
    issame.create(num_pts, num_pts, CV_8U);
    memset(issame.ptr<uint8_t>(), 0, num_pts * num_pts * sizeof(uint8_t));
    active.create(1, num_pts, CV_8U);
    active.setTo(1);

    polarity.create(num_pts, SaddlePointType::NumPolarities,
                    SaddlePointType::polarityStorage);
    if (!deltilleGrid)
      isdiff.create(num_pts, num_pts, CV_8U);

    for (int i = 0; i < num_pts; ++i)
      pts[i].computePolarities(
          polarity.ptr<typename SaddlePointType::PolarityStorageType>(i));

    // limit the number of closest neighbours taken into account
    max_nn = std::min(detector_params.max_nearest_neighbours, int(pts.size()));
    int max_close_nn = std::min(num_pts, deltilleGrid ? 7 : 16);
    for (int i = 0; i < num_pts; ++i) {
      double *d = dist.ptr<double>(i), *u = vec_u.ptr<double>(i),
             *v = vec_v.ptr<double>(i), *dst_u = dist_u.ptr<double>(i),
             *dst_v = dist_v.ptr<double>(i);
      uint8_t *sm = issame.ptr<uint8_t>(i);
      // deltille grid does not have different polarity crossings
      uint8_t *df = deltilleGrid ? sm : isdiff.ptr<uint8_t>(i);
      for (size_t j = 0; j < pts.size(); ++j, ++u, ++v, ++dst_u, ++dst_v, ++d) {
        const double du = pts[j].x - pts[i].x;
        const double dv = pts[j].y - pts[i].y;
        *dst_u = du;
        *dst_v = dv;
        const double dst = sqrt(du * du + dv * dv);
        // normalize to unit vector
        *u = du / (dst + std::numeric_limits<double>::epsilon());
        *v = dv / (dst + std::numeric_limits<double>::epsilon());
        *d = dst;
      }
      // after each row... find nearest neighbours
      partial_sort_indexes(dist.ptr<double>(i), NN.ptr<int>(i), pts.size(),
                           max_nn);
      // only compare polarities with the spatially close buddies
      int *nn = NN.ptr<int>(i);
      for (int j = 0; j < max_nn; ++j) {
        // polarity constraints...
        SaddlePointType::comparePolarities(
            polarity.ptr<typename SaddlePointType::PolarityStorageType>(i),
            polarity.ptr<typename SaddlePointType::PolarityStorageType>(nn[j]),
            sm[nn[j]], df[nn[j]]);
        sm[nn[j]] = (i != nn[j] && sm[nn[j]] == 1) ? 1 : 0;
      }
    }
    if (!deltilleGrid) {
      isdiff = isdiff.mul(isdiff.t());
      // select perspective sample points, we want them to have at least
      // num_good same and num_good different neighbors
      int num_good = 3;
      keypoints.clear();
      keypoints.reserve(num_pts);
      for (int j = 0; j < num_pts; ++j) {
        int *idxs = NN.ptr<int>(j);
        uint8_t *sm = issame.ptr<uint8_t>(j);
        uint8_t *df = isdiff.ptr<uint8_t>(j);
        int good1 = 0, good2 = 0;
        for (int i = 1; i < max_close_nn; ++i) {
          const int id = idxs[i];
          good1 += sm[id];
          good2 += df[id];
          if (good1 >= num_good && good2 >= num_good) {
            keypoints.push_back(j);
            break;
          }
        }
      }
    }
    return true;
  }

  bool checkTriangleConsistency(const SaddlePointType &pt0, int i00, int i11,
                                int pt, double &mean) const {
    const double *dst_u = dist_u.ptr<double>(i00);
    const double *dst_v = dist_v.ptr<double>(i00);
    int cnt = 0;
    mean = 0;
    // yes... go, back to image check consistency of the edge
    double minI = DBL_MAX, maxI = -DBL_MAX;
    for (double triy = 0.1; triy <= 0.9; triy += 0.05) {
      for (double trix = 0.1; trix <= (0.9 - triy); trix += 0.05) {
        double x = pt0.x + trix * dst_u[pt] + triy * dst_u[i11],
               y = pt0.y + trix * dst_v[pt] + triy * dst_v[i11];
        int x0 = int(x), y0 = int(y);
        if (x0 < 0 || y0 < 0 || x0 > width - 2 || y0 > height - 2)
          continue;
        double xw = x - x0, yw = y - y0;
        const InputImageType *row0 = input.ptr<InputImageType>(y0),
                             *row1 = input.ptr<InputImageType>(y0 + 1);
        double I = (1.0 - xw) * (1.0 - yw) * row0[x0] +
                   (xw) * (1.0 - yw) * row0[x0 + 1] +
                   (1.0 - xw) * (yw)*row1[x0] + (xw) * (yw)*row1[x0 + 1];

        mean += I;
        ++cnt;
        if (I < minI)
          minI = I;
        if (I > maxI)
          maxI = I;
      }
    }
    mean /= cnt;
// check quality of the triangle
#if DEBUG_INDEXING >= 5
    cv::cvtColor(input, DEBUG, CV_GRAY2RGB);
    cv::line(
        DEBUG, cv::Point(pt0.x * 65536, pt0.y * 65536),
        cv::Point((pt0.x + dst_u[i11]) * 65536, (pt0.y + dst_v[i11]) * 65536),
        cv::Scalar(0, 0, 255), 1, CV_AA, 16);
    cv::line(
        DEBUG, cv::Point(pt0.x * 65536, pt0.y * 65536),
        cv::Point((pt0.x + dst_u[pt]) * 65536, (pt0.y + dst_v[pt]) * 65536),
        cv::Scalar(0, 255, 0), 1, CV_AA, 16);
    cv::line(
        DEBUG,
        cv::Point((pt0.x + dst_u[i11]) * 65536, (pt0.y + dst_v[i11]) * 65536),
        cv::Point((pt0.x + dst_u[pt]) * 65536, (pt0.y + dst_v[pt]) * 65536),
        cv::Scalar(255, 0, 0), 1, CV_AA, 16);
    if (maxI - minI < detector_params.rectangle_consistency_threshold)
      printf(
          "Accepting quad candidate: minI: %.3lf, maxI: %.3lf, thresh: %.3lf\n",
          minI, maxI, detector_params.rectangle_consistency_threshold);
    else
      printf(
          "Rejecting quad candidate: minI: %.3lf, maxI: %.3lf, thresh: %.3lf\n",
          minI, maxI, detector_params.rectangle_consistency_threshold);
    cv::imshow("output", DEBUG);
#if DEBUG_INDEXING == 6
    cv::waitKey(0);
#else
    cv::waitKey(1);
#endif
#endif
    return (maxI - minI < detector_params.rectangle_consistency_threshold);
  }

  std::vector<bool> computeShadowMask(const std::vector<int> &closenn,
                                      const double *u, const double *v) const {
    // compute shadowing mask
    const double shadow_threshold = cos(detector_params.shadow_angle);
    std::vector<bool> shadow_mask(closenn.size(), false);
    for (size_t r = 0; r < closenn.size(); ++r)
      for (size_t c = r + 1; c < closenn.size(); ++c)
        if (u[closenn[r]] * u[closenn[c]] + v[closenn[r]] * v[closenn[c]] >
            shadow_threshold)
          // remove column indexed nn
          shadow_mask[c] = true;
    return shadow_mask;
  }

  bool initialQuadSelection(const std::vector<int> &idxs,
                            std::vector<Quad> &quads) const {
    if (idxs.empty())
      return false;
    int i00 = idxs[int(double(rand() * (idxs.size())) / (1.0 + RAND_MAX))];
    double polarity_threshold =
        cos(detector_params.rectangle_polarity_angle_threshold);

    const uint8_t *sm = issame.ptr<uint8_t>(i00);
    const uint8_t *df = isdiff.ptr<uint8_t>(i00);
    const std::vector<int> closenn =
        getClosestNNs(i00, max_nn, [sm, df](int nn) -> bool {
          return (sm[nn] + df[nn]) > 0;
        });
    const double *u = vec_u.ptr<double>(i00), *v = vec_v.ptr<double>(i00);
    std::vector<bool> shadow_mask = computeShadowMask(closenn, u, v);

    // get non-shadowed close nn of i00 and split them into same and diff
    // polarity
    std::vector<int> same, diff;
    same.reserve(closenn.size());
    diff.reserve(closenn.size());
    for (size_t c = 0; c < closenn.size(); ++c) {
      if (!shadow_mask[c]) {
        int cnn = closenn[c];
        if (sm[cnn] > 0)
          same.push_back(cnn);
        if (df[cnn] > 0)
          diff.push_back(cnn);
      }
    }

    if (same.empty() || diff.empty())
      return false;

    // looking for the best diagonal candidate
    for (size_t i = 0; i < same.size(); ++i) {
      int i11 = same[i];

      // get all non-diagonal candidates and check if they have correct polarity
      // (different from i00 and i11)
      std::vector<int> candidates;
      candidates.reserve(2);
      double mean;
      for (size_t j = 0; j < diff.size(); ++j) {
        int didx = diff[j];
        // if candidate has proper polarity, go back to image check consistency
        // of the triangle
        if (isdiff.at<uint8_t>(i11, didx) &&
            checkTriangleConsistency(pts[i00], i00, i11, didx, mean))
          candidates.push_back(didx);
      }

      if (candidates.size() != 2 ||
          issame.at<uint8_t>(candidates[0], candidates[1]) == 0)
        // not 2 candidates or they do not have the same polarity
        continue;
      // ok we have two candidates to finish the quad, check if they are in
      // expected location
      double cost1 = u[candidates[0]] * vec_u.at<double>(candidates[1], i11) +
                     v[candidates[0]] * vec_v.at<double>(candidates[1], i11),
             cost2 = u[candidates[1]] * vec_u.at<double>(candidates[0], i11) +
                     v[candidates[1]] * vec_v.at<double>(candidates[0], i11);
      if (cost1 < polarity_threshold || cost2 < polarity_threshold)
        continue;
      cost1 = u[i11] * v[candidates[0]] - v[i11] * u[candidates[0]];
      cost2 = u[i11] * v[candidates[1]] - v[i11] * u[candidates[1]];
      int i10, i01;
      if (cost1 > 0 && cost2 < 0) {
        i10 = candidates[0];
        i01 = candidates[1]; // axis changed by Hyowon
      } else if (cost1 < 0 && cost2 > 0) {
        i10 = candidates[1];
        i01 = candidates[0]; // axis changed by Hyowon
      } else
        continue;
      // final check, do angles play well together to form a quad?
      double a1 = atan2(vec_v.at<double>(i00, i01), vec_u.at<double>(i00, i01));
      double a2 = atan2(vec_v.at<double>(i01, i11), vec_u.at<double>(i01, i11));
      double a3 = atan2(vec_v.at<double>(i11, i10), vec_u.at<double>(i11, i10));
      double a4 = atan2(vec_v.at<double>(i10, i00), vec_u.at<double>(i10, i00));
      double cost =
          mod2pi(a2 - a1) + mod2pi(a3 - a2) + mod2pi(a4 - a3) + mod2pi(a1 - a4);

      if (fabs(cost - 2 * M_PI) > 1.0e-5 && fabs(cost + 2 * M_PI) > 1.0e-5)
        continue;

      quads.push_back(Quad(i00, i01, i10, i11));
    }
    return true;
  }

  bool initialDeltilleQuadSelection(const std::vector<int> &idxs,
                                    std::vector<Quad> &quads) const {
    if (idxs.empty())
      return false;

    int i00 = idxs[int(double(rand() * (idxs.size())) / (1.0 + RAND_MAX))];
    // get point's stats
    const uint8_t *sm = issame.ptr<uint8_t>(i00);
    const std::vector<int> &closenn = getClosestNNs(
        i00, max_nn, [&sm](int idx) -> bool { return sm[idx] > 0; });

    const double tri_angle_threshold =
        sin(detector_params.rectangle_polarity_angle_threshold);
    const double triangle_min_angle_threshold =
        cos(detector_params.rectangle_polarity_angle_threshold);
    const double *u = vec_u.ptr<double>(i00), *v = vec_v.ptr<double>(i00);
    std::vector<bool> shadow_mask = computeShadowMask(closenn, u, v);

    // get non-shadowed close nn of i00 of same polarity
    std::vector<int> same;
    same.reserve(closenn.size());
    for (size_t c = 0; c < closenn.size(); ++c)
      if (!shadow_mask[c]) {
        int cnn = closenn[c];
        if (sm[cnn] > 0)
          same.push_back(cnn);
      }

    if (same.empty()) {
#if DEBUG_INDEXING == 5
      printf("Rejecting initial point %d: no same neighbours found\n", i00);
#endif
      return false;
    }
    // looking for the best diagonal quad candidate out of the same polarity
    // close neighbors
    for (size_t i = 0; i < same.size(); ++i) {
      int i11 = same[i];
      std::vector<int> candidates;
      const uint8_t *sm2 = issame.ptr<uint8_t>(i11);

      for (size_t c = 0; c < closenn.size(); ++c) {
        int s2 = closenn[c];
        if (sm2[s2] > 0) {
          double cost1 = u[i11] * vec_v.at<double>(i11, s2) -
                         v[i11] * vec_u.at<double>(i11, s2);
          double cost2 = u[s2] * vec_v.at<double>(i11, s2) -
                         v[s2] * vec_u.at<double>(i11, s2);
          double cost3 = u[i11] * v[s2] - v[i11] * u[s2];
          if (fabs(cost1) > tri_angle_threshold &&
              fabs(cost2) > tri_angle_threshold &&
              fabs(cost3) > tri_angle_threshold) {
            // check for sign consistency
            if ((cost1 < 0 && cost2 < 0 && cost3 < 0) ||
                (cost1 > 0 && cost2 > 0 && cost3 > 0))
              candidates.push_back(s2);
          }
        }
      }
      if (candidates.size() < 2) {
#if DEBUG_INDEXING == 5
        printf("Rejecting diagonal %d, %d: no suitable quad candidates.\n", i00,
               i11);
#endif
        continue;
      }

      // get best two candidates for of diagonal points, check homogeneity
      std::vector<double> mean(candidates.size());
      int good = 0;
      for (size_t j = 0; j < candidates.size(); ++j)
        if (checkTriangleConsistency(pts[i00], i00, i11, candidates[j],
                                     mean[good])) {
          candidates[good] = candidates[j];
          good++;
        }

      // triangles need to be homogeneous and different enough from each other
      if (good == 2 &&
          fabs(mean[0] - mean[1]) >=
              detector_params.triangle_consistency_threshold) {
        double cost1 = u[candidates[0]] * vec_u.at<double>(candidates[1], i11) +
                       v[candidates[0]] * vec_v.at<double>(candidates[1], i11),
               cost2 = u[candidates[1]] * vec_u.at<double>(candidates[0], i11) +
                       v[candidates[1]] * vec_v.at<double>(candidates[0], i11);
        // do not allow too sharp (<30deg) triangles for initial quad
        if (cost1 < triangle_min_angle_threshold ||
            cost2 < triangle_min_angle_threshold)
          continue;
        cost1 = u[i11] * v[candidates[0]] - v[i11] * u[candidates[0]];
        cost2 = u[i11] * v[candidates[1]] - v[i11] * u[candidates[1]];
        // order the quad properly
        if (cost1 > 0 && cost2 < 0)
          quads.push_back(Quad(i00, i11, candidates[1], candidates[0]));
        else if (cost1 < 0 && cost2 > 0)
          quads.push_back(Quad(i00, i11, candidates[0], candidates[1]));
      }
    }
    return true;
  }

  std::vector<bool> computeShadowMaskIgnoreVisitedAndCheckPolarity(
      const std::vector<int> &closenn, const std::vector<bool> &visited,
      const uint8_t *pol, const double *u, const double *v) const {
    const double shadow_threshold = cos(detector_params.shadow_angle);
    std::vector<bool> shadow_mask(closenn.size(), true);
    for (size_t idx = 0; idx < closenn.size(); ++idx) {
      int n = closenn[idx];
      // ignore visited and ignore same polarity
      if (visited[n] || (!pol[n]))
        shadow_mask[idx] = false;
    }
    // compute shadowing mask
    for (size_t r = 0; r < closenn.size(); ++r)
      for (size_t c = r + 1; c < closenn.size(); ++c)
        if (shadow_mask[r] &&
            u[closenn[r]] * u[closenn[c]] + v[closenn[r]] * v[closenn[c]] >
                shadow_threshold)
          // remove column indexed nn
          shadow_mask[c] = false;
    return shadow_mask;
  }

  uint64_t computeEdgeHash(int n, int id0, int id1) {
    return uint64_t(n) + num_pts * (uint64_t(id0) + num_pts * uint64_t(id1));
  }

  bool verifyEdge(int n, int id0, int id1) {
    uint64_t hash = computeEdgeHash(n, id0, id1);
    bool result = false;
    if (edge_cache.count(hash) != 0) {
#if DEBUG_INDEXING > 2
      printf("Cache hit: %d, %d, %d (%lu) = %s\n", n, id0, id1, hash,
             edge_cache.at(hash) ? "true" : "false");
#endif
      result = edge_cache.at(hash);
    }
#if DEBUG_INDEXING <= 2
    else {
#endif
      cv::Mat u = (cv::Mat_<double>(3, 1) << pts[n].x, pts[id0].x, pts[id1].x);
      cv::Mat v = (cv::Mat_<double>(3, 1) << pts[n].y, pts[id0].y, pts[id1].y);
      cv::Mat A = cv::Mat::ones(3, 3, CV_64F);
      u.copyTo(A.col(0));
      v.copyTo(A.col(1));
      cv::Mat f, b = -u.mul(u) - v.mul(v);
      cv::solve(A, b, f, cv::DECOMP_SVD);

      double x0 = -f.at<double>(0) / 2.0;
      double y0 = -f.at<double>(1) / 2.0;
      double r0 = sqrt(x0 * x0 + y0 * y0 - f.at<double>(2));
      double theta1 = atan2(pts[n].y - y0, pts[n].x - x0);
      double theta2 = atan2(pts[id0].y - y0, pts[id0].x - x0);
      double da = mod(theta1 - theta2, 2.0 * M_PI);
      double dtheta = mod(2.0 * da, 2.0 * M_PI) - da;

      // yes... go, back to image check consistency of the edge
      double minDiffI = DBL_MAX, maxDiffI = DBL_MIN;
      double sgnI = 0;
      int cnt = 0;
      double minI0 = DBL_MAX, maxI0 = 0;
      double minI1 = DBL_MAX, maxI1 = 0;
      double d = dist.at<double>(id0, n);

      double delta_r_cur = std::max(detector_params.delta_r_min,
                                    std::min(detector_params.delta_r_max,
                                             d * detector_params.delta_r_rel));
      for (double x = detector_params.edge_stability_interval_min;
           x <= detector_params.edge_stability_interval_max;
           x += detector_params.edge_stability_interval_step) {
        double a = theta2 + dtheta * x, ca = cos(a), sa = sin(a);
        double u0 = x0 + (r0 - delta_r_cur) * ca,
               u1 = x0 + (r0 + delta_r_cur) * ca;
        double v0 = y0 + (r0 - delta_r_cur) * sa,
               v1 = y0 + (r0 + delta_r_cur) * sa;

        ++cnt;

        // bug fixed by Hyowon
        if (u0 < 0 || u0 > width - 1 || v0 < 0 || v0 > height - 1 || u1 < 0 ||
            u1 > width - 1 || v1 < 0 || v1 > height - 1)
          break;

        // interpolate one side of the edge
        int xi = int(u0), yi = int(v0);
        double xw = u0 - xi, yw = v0 - yi;

        double I0 = (1.0 - xw) * (1.0 - yw) * input.at<InputImageType>(yi, xi) +
                    xw * (1.0 - yw) * input.at<InputImageType>(yi, xi + 1) +
                    (1.0 - xw) * yw * input.at<InputImageType>(yi + 1, xi) +
                    xw * yw * input.at<InputImageType>(yi + 1, xi + 1);
        if (I0 < minI0)
          minI0 = I0;
        if (I0 > maxI0)
          maxI0 = I0;
        // and other side of the edge
        xi = int(u1), yi = int(v1);
        xw = u1 - xi, yw = v1 - yi;
        double I1 = (1.0 - xw) * (1.0 - yw) * input.at<InputImageType>(yi, xi) +
                    xw * (1.0 - yw) * input.at<InputImageType>(yi, xi + 1) +
                    (1.0 - xw) * yw * input.at<InputImageType>(yi + 1, xi) +
                    xw * yw * input.at<InputImageType>(yi + 1, xi + 1);
        if (I1 < minI1)
          minI1 = I1;
        if (I1 > maxI1)
          maxI1 = I1;

        double DiffI = I0 - I1;

        sgnI += (DiffI < 0) ? -1 : 1;
        DiffI = fabs(DiffI);
        if (DiffI < minDiffI)
          minDiffI = DiffI;
        if (DiffI > maxDiffI)
          maxDiffI = DiffI;
      }
      result =
          (minDiffI >
               detector_params.edge_stability_threshold && // absolute threshold
                                                           // on edge contrast
           minDiffI > maxDiffI / 2 && // the minimum gradient is at least half
                                      // of maximum gradient
           std::abs(sgnI) == cnt &&   // there are no sign flips
           (maxI0 - minI0) < detector_params.edge_consistency_threshold &&
           (maxI1 - minI1) < detector_params.edge_consistency_threshold);
#if DEBUG_INDEXING <= 2
    }
#endif
    if (result) {
#if DEBUG_INDEXING > 2
      printf("Accepted %d from: id0: %d, id1: %d, minI: %f, maxI: %f, cnt: %d, "
             "diffI0: %f, diffI1: %f\n",
             n, id0, id1, minDiffI, maxDiffI, cnt, maxI0 - minI0,
             maxI1 - minI1);
#endif
#if DEBUG_INDEXING > 2
      for (double x = detector_params.edge_stability_interval_min;
           x <= detector_params.edge_stability_interval_max;
           x += detector_params.edge_stability_interval_step) {
        double a = theta2 + dtheta * x, ca = cos(a), sa = sin(a);
        double u0 = x0 + (r0 - delta_r_cur) * ca,
               u1 = x0 + (r0 + delta_r_cur) * ca;
        double v0 = y0 + (r0 - delta_r_cur) * sa,
               v1 = y0 + (r0 + delta_r_cur) * sa;
        // bug fixed by Hyowon
        if (u0 < 0 || u0 > width - 1 || v0 < 0 || v0 > height - 1 || u1 < 0 ||
            u1 > width - 1 || v1 < 0 || v1 > height - 1)
          break;

        // interpolate
        int xi = u0, yi = v0;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(0, 255, 255);
        xi = u1, yi = v1;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(0, 255, 255);
      }
#if DEBUG_INDEXING == 4
      imshow("output", DEBUG);
      cv::waitKey(0);
#endif
      for (double x = detector_params.edge_stability_interval_min;
           x <= detector_params.edge_stability_interval_max;
           x += detector_params.edge_stability_interval_step) {
        double a = theta2 + dtheta * x, ca = cos(a), sa = sin(a);
        double u0 = x0 + (r0 - delta_r_cur) * ca,
               u1 = x0 + (r0 + delta_r_cur) * ca;
        double v0 = y0 + (r0 - delta_r_cur) * sa,
               v1 = y0 + (r0 + delta_r_cur) * sa;
        if (u0 < 0 || u0 > width - 1 || v0 < 0 || v0 > height - 1 || u1 < 0 ||
            u1 > width - 1 || v1 < 0 || v1 > height - 1)
          break;

        // interpolate
        int xi = u0, yi = v0;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(0, 0, 255);
        xi = u1, yi = v1;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(0, 0, 255);
      }
      printf("INS: %d, %d, %d <- true\n", n, id0, id1);
#endif
      edge_cache.insert({hash, true});
      return true;
    }
#if DEBUG_INDEXING > 2
    else {
      printf("Rejected %d from: id0: %d, id1: %d, minI: %f, maxI: %f, cnt: %d, "
             "diffI0: %f, diffI1: %f\n",
             n, id0, id1, minDiffI, maxDiffI, cnt, maxI0 - minI0,
             maxI1 - minI1);

#if DEBUG_INDEXING > 2
      for (double x = detector_params.edge_stability_interval_min;
           x <= detector_params.edge_stability_interval_max;
           x += detector_params.edge_stability_interval_step) {
        double a = theta2 + dtheta * x, ca = cos(a), sa = sin(a);
        double u0 = x0 + (r0 - delta_r_cur) * ca,
               u1 = x0 + (r0 + delta_r_cur) * ca;
        double v0 = y0 + (r0 - delta_r_cur) * sa,
               v1 = y0 + (r0 + delta_r_cur) * sa;
        // bug fixed by Hyowon
        if (u0 < 0 || u0 > width - 1 || v0 < 0 || v0 > height - 1 || u1 < 0 ||
            u1 > width - 1 || v1 < 0 || v1 > height - 1)
          break;

        // interpolate
        int xi = u0, yi = v0;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(255, 255, 0);
        xi = u1, yi = v1;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(255, 255, 0);
      }
#if DEBUG_INDEXING == 4
      imshow("output", DEBUG);
      cv::waitKey(0);
#endif
      for (double x = detector_params.edge_stability_interval_min;
           x <= detector_params.edge_stability_interval_max;
           x += detector_params.edge_stability_interval_step) {
        double a = theta2 + dtheta * x, ca = cos(a), sa = sin(a);
        double u0 = x0 + (r0 - delta_r_cur) * ca,
               u1 = x0 + (r0 + delta_r_cur) * ca;
        double v0 = y0 + (r0 - delta_r_cur) * sa,
               v1 = y0 + (r0 + delta_r_cur) * sa;
        if (u0 < 0 || u0 > width - 1 || v0 < 0 || v0 > height - 1 || u1 < 0 ||
            u1 > width - 1 || v1 < 0 || v1 > height - 1)
          break;

        // interpolate
        int xi = u0, yi = v0;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(255, 0, 0);
        xi = u1, yi = v1;
        DEBUG.at<cv::Vec3b>(yi, xi) = cv::Vec3b(255, 0, 0);
      }
      printf("INS: %d, %d, %d <- false\n", n, id0, id1);
#endif
    }
#endif
    edge_cache.insert({hash, false});
    return false;
  }

  bool expandGrid(std::vector<CoordIndex> &coords, std::vector<bool> &visited,
                  cv::Mat &idxmap, int offr, int offc, int ri, int ci,
                  std::vector<int> &id) {
    double vec0[2] = {vec_u.at<double>(id[1], id[0]),
                      vec_v.at<double>(id[1], id[0])};
    // get valid nn indices
    const std::vector<int> closenn =
        getClosestNNs(id[0], max_nn, [](int idx) -> bool {
          UNUSED(idx);
          return true;
        });
    const double *vu = vec_u.ptr<double>(id[0]), *vv = vec_v.ptr<double>(id[0]);
    // go through nearest neighbors and find suitable candidates
    std::vector<bool> id_mask = computeShadowMaskIgnoreVisitedAndCheckPolarity(
        closenn, visited, deltilleGrid ? issame.ptr(id[0]) : isdiff.ptr(id[0]),
        vu, vv);
    const double triangle_polarity_threshold =
        cos(detector_params.triangle_polarity_angle_threshold);

    for (size_t idfi = 0; idfi < id_mask.size(); ++idfi) {
      if (id_mask[idfi]) {
        int n = closenn[idfi];
        double cost = vec0[0] * vu[n] + vec0[1] * vv[n];
        if (cost < triangle_polarity_threshold)
          continue;

        if (verifyEdge(n, id[0], id[1])) {
          coords.push_back(CoordIndex(ri - offr, ci - offc, n));
          idxmap.at<int>(ri, ci) = n;
          visited[n] = true;
          return true;
        }
      }
    }
    return false;
  }

  void initGridSearch(const Quad &quad, std::vector<CoordIndex> &coords,
                      std::vector<bool> &visited) {
    visited[quad.i00_] = true;
    visited[quad.i01_] = true;
    visited[quad.i10_] = true;
    visited[quad.i11_] = true;

    coords.push_back(CoordIndex(0, 0, quad.i00_));
    coords.push_back(CoordIndex(0, 1, quad.i01_));
    coords.push_back(CoordIndex(1, 0, quad.i10_));
    coords.push_back(CoordIndex(deltilleGrid ? -1 : 1, 1, quad.i11_));
  }

  void initGridIndexMap(const std::vector<CoordIndex> &coords, cv::Mat &idxmap,
                        int &offr, int &offc) {
    int minr = INT_MAX, maxr = INT_MIN, minc = INT_MAX, maxc = INT_MIN;
    for (size_t i = 0; i < coords.size(); i++) {
      if (coords[i].r < minr)
        minr = coords[i].r;
      if (coords[i].c < minc)
        minc = coords[i].c;
      if (coords[i].r > maxr)
        maxr = coords[i].r;
      if (coords[i].c > maxc)
        maxc = coords[i].c;
    }
    offr = int(1 - minr);
    offc = int(1 - minc);
    idxmap.create(int(maxr - minr + 3), int(maxc - minc + 3), CV_32S);
    idxmap.setTo(cv::Scalar(-1));
    // column indices into coords
    for (size_t i = 0; i < coords.size(); ++i)
      idxmap.at<int>(offr + coords[i].r, offc + coords[i].c) = coords[i].id;
  }

  bool growGrid(std::vector<CoordIndex> &coords, std::vector<bool> &visited,
                cv::Mat &idxmap, int offr, int offc) {
    int added = 0;
    std::vector<CoordIndex> old_coords(coords);
    for (size_t i = 0; i < old_coords.size(); ++i) {
      int ci = offc + old_coords[i].c;
      int ri = offr + old_coords[i].r;
      int idx = idxmap.at<int>(ri, ci);
      if (idx != -1) {
        std::vector<int> id(2);
        id[0] = idx;
        // has left... check if there are two occupied fields left of the target
        // field...
        if (idxmap.at<int>(ri, ci + 1) == -1 &&
            idxmap.at<int>(ri, ci - 1) != -1) {
          id[1] = idxmap.at<int>(ri, ci - 1);
          if (expandGrid(coords, visited, idxmap, offr, offc, ri, ci + 1, id))
            added++;
        }
        // has right... check if there are two occupied fields right of the
        // target field...
        if (idxmap.at<int>(ri, ci - 1) == -1 &&
            idxmap.at<int>(ri, ci + 1) != -1) {
          id[1] = idxmap.at<int>(ri, ci + 1);
          if (expandGrid(coords, visited, idxmap, offr, offc, ri, ci - 1, id))
            added++;
        }
        // has up
        if (idxmap.at<int>(ri + 1, ci) == -1 &&
            idxmap.at<int>(ri - 1, ci) != -1) {
          id[1] = idxmap.at<int>(ri - 1, ci);
          if (expandGrid(coords, visited, idxmap, offr, offc, ri + 1, ci, id))
            added++;
        }
        // has down
        if (idxmap.at<int>(ri - 1, ci) == -1 &&
            idxmap.at<int>(ri + 1, ci) != -1) {
          id[1] = idxmap.at<int>(ri + 1, ci);
          if (expandGrid(coords, visited, idxmap, offr, offc, ri - 1, ci, id))
            added++;
        }
        if (deltilleGrid) {
          // has down left
          if (idxmap.at<int>(ri - 1, ci + 1) == -1 &&
              idxmap.at<int>(ri + 1, ci - 1) != -1) {
            id[1] = idxmap.at<int>(ri + 1, ci - 1);
            if (expandGrid(coords, visited, idxmap, offr, offc, ri - 1, ci + 1,
                           id))
              added++;
          }
          // has up right
          if (idxmap.at<int>(ri + 1, ci - 1) == -1 &&
              idxmap.at<int>(ri - 1, ci + 1) != -1) {
            id[1] = idxmap.at<int>(ri - 1, ci + 1);
            if (expandGrid(coords, visited, idxmap, offr, offc, ri + 1, ci - 1,
                           id))
              added++;
          }
        }
      }
    }
    // has the grid changed?
    return added != 0;
  }

  bool findBestGrid(const cv::Size &board_size, cv::Mat &best_grid) {
    cv::Scalar total = cv::sum(active);
    if (total[0] < 16)
      return false;

    if (board_size.area() == 0)
      return false;

#if DEBUG_INDEXING > 1
    cv::Mat bestDebug, origDebug;

    cv::cvtColor(input, DEBUG, CV_GRAY2RGB);
    int cnt = 0;
    for (size_t i = 0; i < pts.size(); i++) {
      if (active.at<uint8_t>(i) != 0) {
        cv::Point2f pt;
        pt.x = pts[i].x;
        pt.y = pts[i].y;
        ++cnt;
        cv::rectangle(DEBUG, cv::Point((pt.x - 3) * 65536, (pt.y - 3) * 65536),
                      cv::Point((pt.x + 3) * 65536, (pt.y + 3) * 65536),
                      cv::Scalar(255, 0, 0), 1, CV_AA, 16);
      }
    }
    printf("Remaining %d active points\n", cnt);
    // show remaining detections
    cv::imshow("output", DEBUG);
    cv::waitKey(0);

    origDebug = DEBUG.clone();
#endif

    best_grid.create(board_size, CV_32S);
    best_grid.setTo(cv::Scalar(-1));

    // now the sequential part...
    int iter = 0, trial = 0, nbest = 0;
    cv::Mat idxmap;

    std::vector<Quad> axis; //    srand(time(NULL));
    const std::vector<int> idxs = getActiveIndices();
    while (trial < detector_params.grid_search_quad_growing_trials &&
           iter < detector_params.grid_search_max_iterations) {
      ++iter;
      axis.clear();

      bool have_quad = deltilleGrid ? initialDeltilleQuadSelection(idxs, axis)
                                    : initialQuadSelection(keypoints, axis);

      if (!have_quad)
        continue;

      if (!axis.empty())
        trial++;

      for (size_t k = 0; k < axis.size(); ++k) {
        std::vector<CoordIndex> coords;
        std::vector<bool> visited(num_pts, false);
        Quad &quad = axis[k];
        initGridSearch(quad, coords, visited);

#if DEBUG_INDEXING > 1
        cv::line(DEBUG,
                 cv::Point(pts[quad.i00_].x * 65536, pts[quad.i00_].y * 65536),
                 cv::Point(pts[quad.i10_].x * 65536, pts[quad.i10_].y * 65536),
                 cv::Scalar(0, 0, 255), 1, CV_AA, 16);
        cv::line(DEBUG,
                 cv::Point(pts[quad.i10_].x * 65536, pts[quad.i10_].y * 65536),
                 cv::Point(pts[quad.i01_].x * 65536, pts[quad.i01_].y * 65536),
                 cv::Scalar(0, 255, 0), 1, CV_AA, 16);
        cv::line(DEBUG,
                 cv::Point(pts[quad.i01_].x * 65536, pts[quad.i01_].y * 65536),
                 cv::Point(pts[quad.i11_].x * 65536, pts[quad.i11_].y * 65536),
                 cv::Scalar(255, 0, 0), 1, CV_AA, 16);
        cv::line(DEBUG,
                 cv::Point(pts[quad.i11_].x * 65536, pts[quad.i11_].y * 65536),
                 cv::Point(pts[quad.i00_].x * 65536, pts[quad.i00_].y * 65536),
                 cv::Scalar(255, 0, 255), 1, CV_AA, 16);
#endif
        while (true) {
          int offr, offc;
          initGridIndexMap(coords, idxmap, offr, offc);
          // std::cout << "IdxMap:" << std::endl << idxmap << std::endl;
          bool grown = growGrid(coords, visited, idxmap, offr, offc);
#if DEBUG_INDEXING > 2
          cv::imshow("output", DEBUG);
          cv::waitKey(0);
#endif
          if (!grown)
            break;
        }

#ifdef DEBUG_INDEXING
// std::cout << "Final IDX map:" << std::endl << idxmap << std::endl;
#endif
        int nnz = 0;
        for (int r = 0; r < idxmap.rows; ++r) {
          int *row = idxmap.ptr<int>(r);
          for (int c = 0; c < idxmap.cols; ++c)
            if (row[c] >= 0)
              ++nnz;
        }
        if (nnz > 0.3 * board_size.area())
          ++trial;
#ifdef DEBUG_INDEXING
        printf("Iters: %d, trials: %d, found: %d points, best: %d\n", iter,
               trial, nnz, nbest);
#endif
        if (nnz > nbest) {
#ifdef DEBUG_INDEXING
          std::cout << "Best IDX map:" << std::endl << idxmap << std::endl;
#endif
          nbest = nnz;
          idxmap.copyTo(best_grid);
          if (nbest == board_size.area())
            break;
#if DEBUG_INDEXING == 2
          cv::imshow("output", DEBUG);
          cv::waitKey(0);
#endif
#if DEBUG_INDEXING > 1
          bestDebug = DEBUG.clone();
#endif
        }
#if DEBUG_INDEXING > 1
        origDebug.copyTo(DEBUG);
#endif
      }
    }

    return true;
  }

  int findBoards(const cv::Size &board_size,
                 std::vector<BoardObservation> &boards, bool best_only) {
    // do all the preprocessing here
    SaddlePointVector refclust;

    // find initial saddles on lowres image...
    detector.findSaddles(refclust);
    precomputePolaritiesAndNN(refclust);

    while (true) {
      BoardObservation obs;
#if DEBUG_TIMING
      auto t3 = high_resolution_clock::now();
#endif
      findBestGrid(board_size, obs.board);
      obs.indexed = false;
      // squeeze all corner locations on this board...
      std::vector<cv::Point2f> &corners = obs.corner_locations;
      int cnt = 0;
      for (int r = 0; r < obs.board.rows; ++r)
        for (int c = 0; c < obs.board.cols; ++c) {
          int &ptidx = obs.board.at<int>(r, c);
          if (ptidx > -1) {
            ++cnt;
            SaddlePointType &pt = refclust[ptidx];
            // inactivate this point location
            active.at<uint8_t>(ptidx) = 0;
            ptidx = int(corners.size());
            corners.push_back(cv::Point2f(float(pt.x), float(pt.y)));
          } else
            corners.push_back(
                cv::Point2f(std::numeric_limits<float>::infinity(),
                            std::numeric_limits<float>::infinity()));
        }
#if DEBUG_INDEXING > 2
      cv::cvtColor(input, DEBUG, CV_GRAY2RGB);
      orp::calibration::drawCheckerboardCorners(DEBUG, obs, 1.5, true);
      imshow("output", DEBUG);
      cv::waitKey(0);
#endif

#ifdef DEBUG_TIMING
      auto t4 = high_resolution_clock::now();
      printf("Indexing took: %.3f ms , found %d corners\n",
             duration_cast<microseconds>(t4 - t3).count() / 1000.0, cnt);
#endif
      // found less than 10 corners on a board... stop here
      if (cnt < 9)
        break;

      boards.push_back(obs);
#ifdef DEBUG_INDEXING
//            std::cout << "FINAL board map:" << std::endl << obs.board<<
//            std::endl;
#endif

      if (cnt < 20 || best_only)
        break;
      // auto rend = std::remove_if(refclust.begin(), refclust.end(),
      // PointIsInf<SaddlePointType>);
      // refclust.resize(rend - refclust.begin());
    }
    detector.finalizeSaddles(boards);
    return int(boards.size());
  }
};
}
}

#endif /* INCLUDE_CHECKERBOARD_DETECTOR_GRIDDETECTORCONTEXT_H_ */
