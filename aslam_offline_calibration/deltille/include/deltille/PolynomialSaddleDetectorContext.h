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
 * PolynomialSaddleDetector.h
 *
 *  Created on: Dec 1, 2016
 *      Author: mperdoch
 */
#pragma once
#ifndef INCLUDE_CHECKERBOARD_DETECTOR_POLYNOMIALSADDLEDETECTORCONTEXT_H_
#define INCLUDE_CHECKERBOARD_DETECTOR_POLYNOMIALSADDLEDETECTORCONTEXT_H_

#include <deltille/DetectorTools.h>
#include <deltille/PolynomialFit.h>

#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>

using namespace std::chrono;

namespace orp {
namespace calibration {

template <typename SaddlePointType, typename InputImageType = uint8_t,
          typename FloatImageType = float>
struct PolynomialSaddleDetectorContext {
private:
  // low res and full size polynomial fits...
  PolynomialFit<SaddlePointType> lowresFitting;
  PolynomialFit<SaddlePointType> fullFitting;

  double scaling;
  int num_iterations;

public:
  // potentially resized original image
  cv::Mat input_lowres;

  // low res and full size of the input image (smoothed)
  cv::Mat lowres;
  cv::Mat full;

public:
  PolynomialSaddleDetectorContext(const cv::Mat &img) {
#ifdef DEBUG_TIMING
    auto t0 = high_resolution_clock::now();
#endif
    preprocessImage(img);
#ifdef DEBUG_TIMING
    auto t1 = high_resolution_clock::now();
    printf("Filtering took: %.3f ms\n",
           duration_cast<microseconds>(t1 - t0).count() / 1000.0);
#endif
    num_iterations = SaddlePointType::isTriangular ? 5 : 20;
  }

  int findSaddles(std::vector<SaddlePointType> &refclust) {
#ifdef DEBUG_TIMING
    auto t0 = high_resolution_clock::now();
#endif
    std::vector<cv::Point> locations;

    getInitialSaddleLocations(input_lowres, locations);
#ifdef DEBUG_TIMING
    auto t1 = high_resolution_clock::now();
    printf("Preprocessing took: %.3f ms\n",
           duration_cast<microseconds>(t1 - t0).count() / 1000.0);
#endif
    std::vector<SaddlePointType> refined;
    lowresFitting.saddleSubpixelRefinement(lowres, locations, refined,
                                           num_iterations, true);

    if (!SaddlePointType::isTriangular) {
      // early prefilter outliers
      double smax = -DBL_MAX;
      for (size_t i = 0; i < refined.size(); ++i)
        if (smax < refined[i].s)
          smax = refined[i].s;

      const double min_angle_width = 15.0 * M_PI / 180;
      const double max_angle_width = 75.0 * M_PI / 180;
      smax *= detector_params.rectangular_saddle_threshold;
      for (size_t i = 0; i < refined.size(); ++i) {
        if (refined[i].s < smax || std::abs(refined[i].a2) < min_angle_width ||
            std::abs(refined[i].a2) > max_angle_width)
          refined[i].x = refined[i].y = std::numeric_limits<double>::infinity();
      }
    }
    // squeeze out diverged points
    refined.resize(std::remove_if(refined.begin(), refined.end(),
                                  PointIsInf<SaddlePointType>) -
                   refined.begin());

#ifdef DEBUG_TIMING
    auto t2 = high_resolution_clock::now();
    printf("Refinement took: %.3f ms\n",
           duration_cast<microseconds>(t2 - t1).count() / 1000.0);
#endif
    int num_clusters = 0;
    std::vector<int> cluster_ids;

    std::vector<typename SaddlePointType::ClusterDescType> cluster_stats;
    clusterPoints2(refined, lowres.size(), cluster_ids, cluster_stats,
                   num_clusters, 1.0);
#ifdef DEBUG_TIMING
    auto t3 = high_resolution_clock::now();
    printf("Clustering took: %.3f ms\n",
           duration_cast<microseconds>(t3 - t2).count() / 1000.0);
#endif
    refclust.resize(cluster_stats.size());
    for (size_t i = 0; i < cluster_stats.size(); i++)
      refclust[i] = cluster_stats[i];

#ifdef DEBUG_INDEXING
    cv::cvtColor(input_lowres, DEBUG, CV_GRAY2RGB);
    for (size_t i = 0; i < refclust.size(); i++) {
      cv::Point2f pt;
      pt.x = refclust[i].x;
      pt.y = refclust[i].y;
      cv::rectangle(DEBUG, cv::Point((pt.x - 3) * 65536, (pt.y - 3) * 65536),
                    cv::Point((pt.x + 3) * 65536, (pt.y + 3) * 65536),
                    cv::Scalar(255, 0, 0), 1, CV_AA, 16);
      if (!SaddlePointType::isTriangular)
        refclust[i].plotPolarities(DEBUG, 3.0 / scaling);
    }
#endif

#ifdef DEBUG_TIMING
    auto t4 = high_resolution_clock::now();
#endif
    if (SaddlePointType::isTriangular) {
      // a bit hackier monkey saddle second filter, that checks if points would
      // converge to the same location with larger scale...
      PolynomialFit<SaddlePointType> tempFitting;
      tempFitting.initSaddleFitting(
          lowresFitting.getHalfKernelSize() +
          detector_params.deltille_stability_kernel_size_increase);

      std::vector<SaddlePointType> tempclust;
      cv::Mat temp;
      cv::filter2D(input_lowres, temp, cv::DataType<FloatImageType>::depth,
                   tempFitting.getSmoothingKernel());
      tempFitting.saddleSubpixelRefinement(temp, refclust, tempclust);
      for (size_t i = 0; i < refclust.size(); i++) {
        double dx = refclust[i].x - tempclust[i].x,
               dy = refclust[i].y - tempclust[i].y;
        if (dx * dx + dy * dy > detector_params.deltille_stability_threshold) {
          refclust[i].x = std::numeric_limits<double>::infinity();
          refclust[i].y = std::numeric_limits<double>::infinity();
        }
      }
    }
    // squeeze out any remaining unwanted points...
    refclust.resize(std::remove_if(refclust.begin(), refclust.end(),
                                   PointIsInf<SaddlePointType>) -
                    refclust.begin());
#ifdef DEBUG_INDEXING
    // show final detections
    for (size_t i = 0; i < refclust.size(); i++) {
      SaddlePointType &s1 = refclust[i];
      int sz = scaling * 3;
      cv::rectangle(DEBUG, cv::Point((s1.x - sz) * 65536, (s1.y - sz) * 65536),
                    cv::Point((s1.x + sz) * 65536, (s1.y + sz) * 65536),
                    cv::Scalar(0, 255, 0), 3.0 / scaling, CV_AA, 16);
      s1.plotPolarities(DEBUG, 3.0 / scaling);
      cv::putText(DEBUG, std::to_string(i), cv::Point((s1.x + sz), (s1.y + sz)),
                  cv::FONT_HERSHEY_PLAIN, 0.8, cv::Scalar(255, 0, 255));
    }
    imshow("detections", DEBUG);
#endif
#ifdef DEBUG_TIMING
    auto t5 = high_resolution_clock::now();
#endif
    return int(refclust.size());
  }

  int finalizeSaddles(std::vector<BoardObservation> &boards) {
#ifdef DEBUG_TIMING
    auto t1 = high_resolution_clock::now();
#endif
    // finally, refine at full scale..    .
    for (size_t i = 0; i < boards.size(); ++i) {
      BoardObservation &obs = boards[i];
      std::vector<cv::Point2f> fullScaleLocations(obs.corner_locations);
      for (auto &pt : fullScaleLocations) {
        pt.x = float((pt.x + 0.5) * scaling - 0.5);
        pt.y = float((pt.y + 0.5) * scaling - 0.5);
      }
      std::vector<SaddlePointType> fullScalePoints(fullScaleLocations.size());
      // spend a bit more time on final refinement... to make sure we do not
      // throw away good points after upscaling
      fullFitting.saddleSubpixelRefinement(
          full, fullScaleLocations, fullScalePoints, num_iterations * 2, false);
      for (size_t c = 0; c < obs.corner_locations.size(); ++c) {
        SaddlePointType &pt = fullScalePoints[c];
        if (PointIsInf(pt)) {
#ifdef DEBUG_INDEXING
          if (!PointIsInf(fullScaleLocations[c])) {
            printf("Refinement failed for: %d, board #%zu, reason: %d, x: "
                   "%.3lf -> %.3lf, y: %.3lf -> %.3lf, det: %.3lf\n",
                   obs.board.at<int>(c / obs.board.cols, c % obs.board.cols), i,
                   int(pt.s), fullScaleLocations[c].x, std::abs(pt.a1),
                   fullScaleLocations[c].y, std::abs(pt.a2), pt.det);
          }
#endif
          obs.corner_locations[c].x = -1.0f;
          obs.corner_locations[c].y = -1.0f;
          obs.board.at<int>(int(c / obs.board.cols), int(c % obs.board.cols)) =
              -1;
        } else {
          obs.corner_locations[c].x = float(pt.x);
          obs.corner_locations[c].y = float(pt.y);
        }
      }
    }
#ifdef DEBUG_TIMING
    auto t2 = high_resolution_clock::now();
    printf("Final refinement took: %.3f ms\n",
           duration_cast<microseconds>(t2 - t1).count() / 1000.0);
#endif
    return 0;
  }

private:
  void stretchIntensities(cv::InputArray input, cv::OutputArray output) {
    // build a lookup table mapping the pixel values [0, 255] to
    // their adjusted gamma values
    double minVal, maxVal;
    cv::minMaxIdx(input, &minVal, &maxVal);
    input.getMat().convertTo(output, input.type(), 255.0 / (maxVal - minVal),
                             -minVal);
  }

  void preprocessImage(cv::InputArray img) {
    cv::Mat gray_img;

    if (img.getMat().channels() == 3)
      // convert to a single channel first, keep original intensity range for
      // now...
      cvtColor(img, gray_img, CV_RGB2GRAY);
    else
      gray_img = img.getMat();

    scaling = 1.0;
    double res = std::max(gray_img.rows, gray_img.cols);
    if (res > detector_params.working_resolution)
      scaling = int(res / detector_params.working_resolution * 10) / 10.0;

    // presmooth full image with the full size kernel to get a smooth one
    int fullres_half_kernel_size =
        int(fmin(7, detector_params.half_kernel_size * scaling + 0.5));
    lowresFitting.initSaddleFitting(detector_params.half_kernel_size);
    fullFitting.initSaddleFitting(fullres_half_kernel_size);

    cv::filter2D(gray_img, full, cv::DataType<FloatImageType>::depth,
                 fullFitting.getSmoothingKernel());

    // resize to lowres, no additional smoothing should be necessary
    cv::resize(gray_img, input_lowres, cv::Size(), 1.0 / scaling, 1.0 / scaling,
               cv::INTER_CUBIC);
    // for initial detection adjust the intensity range to (0-255) for proper
    // thresholding
    stretchIntensities(input_lowres, input_lowres);
    cv::filter2D(input_lowres, lowres, cv::DataType<FloatImageType>::depth,
                 lowresFitting.getSmoothingKernel());
  }

  void getInitialSaddleLocations(const cv::Mat &input,
                                 std::vector<cv::Point> &locations) {
    cv::Mat gauss_img;
    input.convertTo(gauss_img, cv::DataType<FloatImageType>::depth);

    cv::GaussianBlur(gauss_img, gauss_img, cv::Size(7, 7), 1.5, 1.5);
    cv::Mat hessian_img;
    hessianResponse(gauss_img, hessian_img);
    double mn = 0, mx = 0;
    cv::minMaxIdx(hessian_img, &mn, &mx, NULL, NULL);
    cv::Mat mask = hessian_img < mn * 0.01;
    cv::dilate(mask, mask, cv::Mat());
    // remove saddles around the border
    for (int i = 0; i < 2; i++) {
      mask.row(i) = cv::Scalar::all(0);
      mask.row(mask.rows - 1 - i) = cv::Scalar::all(0);
      mask.col(i) = cv::Scalar::all(0);
      mask.col(mask.cols - 1 - i) = cv::Scalar::all(0);
    }
    // get locations and sort them by y and then x
    cv::Mat locs;
    cv::findNonZero(mask, locs);
    if (locs.total() > 0) {
      locations.resize(locs.rows);
      cv::Point *tmp = locs.ptr<cv::Point>();
      std::copy(tmp, tmp + locs.rows, locations.begin());
    }
    std::sort(locations.begin(), locations.end(), point_comparator);
  }
};
}
}

#endif /* INCLUDE_CHECKERBOARD_DETECTOR_POLYNOMIALSADDLEDETECTORCONTEXT_H_ */
