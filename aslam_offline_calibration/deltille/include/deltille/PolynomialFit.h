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
 * PolynomialFit.h
 *
 *  Created on: Nov 29, 2016
 *      Author: mperdoch
 */
#pragma once
#ifndef INCLUDE_CHECKERBOARD_DETECTOR_POLYNOMIALFIT_H_
#define INCLUDE_CHECKERBOARD_DETECTOR_POLYNOMIALFIT_H_

#include <deltille/DetectorParams.h>
#include <deltille/DetectorTools.h>
#include <deltille/utils.h>

#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>

using namespace std::chrono;

namespace orp {
namespace calibration {

struct SaddlePoint {
  typedef SaddleClusterDesc ClusterDescType;
  typedef double PolarityStorageType;
  const static int polarityStorage = CV_64FC1;
  const static bool isTriangular = false;
  enum {
    PolaritySine1 = 0,
    PolarityCosine1 = 1,
    PolaritySine2 = 2,
    PolarityCosine2 = 3,
    NumPolarities // keep this one at the end...
  };

public:
  SaddlePoint() {}
  SaddlePoint(double x, double y) : x(x), y(y) {}
  SaddlePoint(const ClusterDescType &other)
      : x(other.cx), y(other.cy), a1(other.a1), a2(other.a2) {}

  void computePolarities(double *p) const {
    double a_1 = a1 + a2;
    double a_2 = a1 - a2;
    p[PolaritySine1] = sin(a_1);
    p[PolarityCosine1] = cos(a_1);
    p[PolaritySine2] = sin(a_2);
    p[PolarityCosine2] = cos(a_2);
  }

  static void comparePolarities(const PolarityStorageType *p1,
                                const PolarityStorageType *p2, uint8_t &same,
                                uint8_t &diff) {
    const double rectangular_polarity_threshold =
        cos(detector_params.rectangle_polarity_angle_threshold);
    double cost1 = fabs(p1[PolaritySine1] * p2[PolaritySine2] +
                        p1[PolarityCosine1] * p2[PolarityCosine2]);
    double cost2 = fabs(p1[PolaritySine1] * p2[PolaritySine1] +
                        p1[PolarityCosine1] * p2[PolarityCosine1]);
    double cost3 = fabs(p1[PolaritySine2] * p2[PolaritySine2] +
                        p1[PolarityCosine2] * p2[PolarityCosine2]);
    diff = (cost1 > rectangular_polarity_threshold) ? 1 : 0;
    same = (cost2 > rectangular_polarity_threshold &&
            cost3 > rectangular_polarity_threshold)
               ? 1
               : 0;
  }
#if DEBUG_INDEXING > 1
  static void debugComparePolarities(const PolarityStorageType *p1,
                                     const PolarityStorageType *p2, int id1,
                                     int id2, double dx, double dy) {
    printf("pair: %d, %d, dx: %f, dy: %f\n", id1, id2, dx, dy);
    double cost1 = fabs(p1[PolaritySine1] * p2[PolaritySine2] +
                        p1[PolarityCosine1] * p2[PolarityCosine2]);
    double cost2 = fabs(p1[PolaritySine1] * p2[PolaritySine1] +
                        p1[PolarityCosine1] * p2[PolarityCosine1]);
    double cost3 = fabs(p1[PolaritySine2] * p2[PolaritySine2] +
                        p1[PolarityCosine2] * p2[PolarityCosine2]);
    printf("  cost1: %f, cost2: %f, cost3: %f\n", cost1, cost2, cost3);
  }
#endif
  void plotPolarities(cv::Mat &im, double scaling = 1.0) {
    int sz = int(scaling * 3);
    double r = 5.0;
    double a_1 = a1 + a2;
    double a_2 = a1 - a2;
    subpixel_line(im, cv::Point2d(x, y),
                  cv::Point2d(x + cos(a_1) * r, y + sin(a_1) * r),
                  cv::Scalar(255, 0, 0), sz);
    subpixel_line(im, cv::Point2d(x, y),
                  cv::Point2d(x + cos(a_2) * r, y + sin(a_2) * r),
                  cv::Scalar(255, 0, 0), sz);
  }

public:
  double x, y;
  PolarityStorageType a1, a2;
  double s, det;
};

struct MonkeySaddlePoint {
  typedef MonkeySaddleClusterDesc ClusterDescType;
  typedef std::complex<double> PolarityStorageType;
  const static int polarityStorage = CV_64FC2;
  const static bool isTriangular = true;
  enum {
    PolaritySine1 = 0,
    PolarityCosine1 = 1,
    PolaritySine2 = 2,
    PolarityCosine2 = 3,
    PolaritySine3 = 4,
    PolarityCosine3 = 5,
    NumPolarities // keep this one at the end...
  };

public:
  MonkeySaddlePoint() {}
  MonkeySaddlePoint(double x, double y) : x(x), y(y) {}
  MonkeySaddlePoint(const ClusterDescType &other)
      : x(other.cx), y(other.cy), a1(other.a1), a2(other.a2), a3(other.a3) {}

  void computePolarities(PolarityStorageType *p) const {
    // printf("pol: a1: %f%+fi, a2: %f%+fi, a3: %f%+fi\n", a1.real(), a1.imag(),
    // a2.real(), a2.imag(), a3.real(), a3.imag());
    p[PolaritySine1] = sin(a1);
    p[PolarityCosine1] = cos(a1);
    p[PolaritySine2] = sin(a2);
    p[PolarityCosine2] = cos(a2);
    p[PolaritySine3] = sin(a3);
    p[PolarityCosine3] = cos(a3);
  }

  static void comparePolarities(const PolarityStorageType *p1,
                                const PolarityStorageType *p2, uint8_t &same,
                                uint8_t &dummy) {
    const double triangle_polarity_threshold =
        cos(10.0 / 180 * M_PI); // cos(detector_params.triangle_polarity_angle);

    PolarityStorageType cost1 = p1[PolaritySine1] * p2[PolaritySine1] +
                                p1[PolarityCosine1] * p2[PolarityCosine1];
    PolarityStorageType cost2 = p1[PolaritySine2] * p2[PolaritySine2] +
                                p1[PolarityCosine2] * p2[PolarityCosine2];
    PolarityStorageType cost3 = p1[PolaritySine3] * p2[PolaritySine3] +
                                p1[PolarityCosine3] * p2[PolarityCosine3];
    // 1-3, 2-1, 3-2
    bool sm = std::abs(cost1) > triangle_polarity_threshold &&
              std::abs(cost2) > triangle_polarity_threshold &&
              std::abs(cost3) > triangle_polarity_threshold;
    same = sm ? 1 : 0;
    UNUSED(dummy);
  }
#if DEBUG_INDEXING > 1
  static void debugComparePolarities(const PolarityStorageType *p1,
                                     const PolarityStorageType *p2, int id1,
                                     int id2, double dx, double dy) {
    // const double triangle_polarity_threshold = cos(10.0/180*M_PI);
    // //cos(detector_params.triangle_polarity_angle);

    PolarityStorageType cost1 = p1[PolaritySine1] * p2[PolaritySine1] +
                                p1[PolarityCosine1] * p2[PolarityCosine1];
    PolarityStorageType cost2 = p1[PolaritySine2] * p2[PolaritySine2] +
                                p1[PolarityCosine2] * p2[PolarityCosine2];
    PolarityStorageType cost3 = p1[PolaritySine3] * p2[PolaritySine3] +
                                p1[PolarityCosine3] * p2[PolarityCosine3];
    printf("pair: %d, %d, dx: %f, dy: %f\n", id1, id2, dx, dy);
    for (int i = 0; i < 6; ++i)
      printf("  p[%d]: %f%+fi   <-> %f%+fi\n", i, p1[i].real(), p1[i].imag(),
             p2[i].real(), p2[i].imag());
    printf("  cost1: %f%+fi, cost2: %f%+fi, cost3: %f%+fi\n", cost1.real(),
           cost1.imag(), cost2.real(), cost2.imag(), cost3.real(),
           cost3.imag());
  }
#endif
  void plotPolarities(cv::Mat &im, double scaling = 1.0) {
    int sz = int(scaling);
    double r = 25.0;
    subpixel_line(im, cv::Point2d(x, y), cv::Point2d(x + std::abs(cos(a1)) * r,
                                                     y + std::abs(sin(a1)) * r),
                  cv::Scalar(255, 0, 0), sz);
    r = 20.0;
    subpixel_line(im, cv::Point2d(x, y), cv::Point2d(x + std::abs(cos(a2)) * r,
                                                     y + std::abs(sin(a2)) * r),
                  cv::Scalar(255, 255, 0), sz);
    r = 15.0;
    subpixel_line(im, cv::Point2d(x, y), cv::Point2d(x + std::abs(cos(a3)) * r,
                                                     y + std::abs(sin(a3)) * r),
                  cv::Scalar(255, 255, 0), sz);
  }

public:
  double x, y;
  PolarityStorageType a1, a2, a3;
  double s, det;
};

struct MonkeySaddlePointSpherical {
  typedef MonkeySaddleClusterDesc ClusterDescType;
  typedef cv::Vec3d PolarityStorageType;
  const static int polarityStorage = CV_64FC3;
  const static bool isTriangular = true;
  enum {
    Polarity1 = 0,
    Polarity2 = 1,
    Polarity3 = 2,
    NumPolarities // keep this one at the end...
  };

public:
  MonkeySaddlePointSpherical() {}
  MonkeySaddlePointSpherical(double x, double y) : x(x), y(y) {}
  MonkeySaddlePointSpherical(const ClusterDescType &other)
      : x(other.cx), y(other.cy), a1(other.a1), a2(other.a2), a3(other.a3) {}

  void complexToSpherical(const std::complex<double> &angle,
                          PolarityStorageType &vec) const {
    const double ci = cos(angle.imag());
    vec[0] = ci * cos(angle.real());
    vec[1] = ci * sin(angle.real());
    vec[2] = sin(angle.imag());
  }

  void computePolarities(PolarityStorageType *p) const {
    // convert complex angles to vectors on unit sphere
    complexToSpherical(a1, p[0]);
    complexToSpherical(a2, p[1]);
    complexToSpherical(a3, p[2]);
  }

  static bool comparePolaritiesUnderRotation(const PolarityStorageType *p1,
                                             const PolarityStorageType *p2,
                                             int rotation) {
    const double triangle_polarity_threshold =
        cos(10.0 / 180 * M_PI); // cos(detector_params.triangle_polarity_angle);
    return std::abs(p1[0].dot(p2[rotation])) > triangle_polarity_threshold &&
           std::abs(p1[1].dot(p2[(rotation + 1) % 3])) >
               triangle_polarity_threshold &&
           std::abs(p1[2].dot(p2[(rotation + 2) % 3])) >
               triangle_polarity_threshold;
  }

  static void comparePolarities(const PolarityStorageType *p1,
                                const PolarityStorageType *p2, uint8_t &same,
                                uint8_t &dummy) {
    bool sm = comparePolaritiesUnderRotation(p1, p2, 0) ||
              comparePolaritiesUnderRotation(p1, p2, 1) ||
              comparePolaritiesUnderRotation(p1, p2, 2);
    same = sm ? 1 : 0;
    UNUSED(dummy);
  }
#if DEBUG_INDEXING > 1
  static void debugComparePolarities(const PolarityStorageType *p1,
                                     const PolarityStorageType *p2, int id1,
                                     int id2, double dx, double dy) {
    std::cout << p1[0] << " vs. " << p2[0] << std::endl;
    std::cout << p1[1] << " vs. " << p2[1] << std::endl;
    std::cout << p1[2] << " vs. " << p2[2] << std::endl;
  }
#endif
  void plotPolarities(cv::Mat &im, double scaling = 1.0) {
    PolarityStorageType p[3];
    computePolarities(p);
    int sz = int(scaling);
    double r, z;
    r = 25.0;
    z = 1.0 - p[0].val[2];
    subpixel_line(im, cv::Point2d(x, y),
                  cv::Point2d(x + p[0].val[0] / z * r, y + p[0].val[1] / z * r),
                  cv::Scalar(255, 0, 0), sz);
    r = 20.0;
    z = 1.0 - p[1].val[2];
    subpixel_line(im, cv::Point2d(x, y),
                  cv::Point2d(x + p[1].val[0] / z * r, y + p[1].val[1] / z * r),
                  cv::Scalar(255, 255, 0), sz);
    r = 15.0;
    z = 1.0 - p[2].val[2];
    subpixel_line(im, cv::Point2d(x, y),
                  cv::Point2d(x + p[2].val[0] / z * r, y + p[2].val[1] / z * r),
                  cv::Scalar(255, 255, 0), sz);
  }

public:
  double x, y;
  std::complex<double> a1, a2, a3; // these can be possibly complex numbers ...
  double s, det;
};

const double separable_kernel_coefs[][31] = {
    {},
    {},
    {0.10805848, 0.23395099, 0.31932424, 0.23395099, 0.10805848},
    {0.05745145, 0.13032769, 0.19602977, 0.23725367, 0.19602977, 0.13032769,
     0.05745145},
    {0.03399957, 0.08105924, 0.12712243, 0.16627552, 0.18958185, 0.16627552,
     0.12712243, 0.08105924, 0.03399957},
    {0.02259464, 0.05455135, 0.08703301, 0.11746379, 0.14276789, 0.15724555,
     0.14276789, 0.11746379, 0.08703301, 0.05455135, 0.02259464},
    {0.01579519, 0.03862299, 0.06283421, 0.08618732, 0.10751836, 0.12494305,
     0.13461367, 0.12494305, 0.10751836, 0.08618732, 0.06283421, 0.03862299,
     0.01579519},
    {0.01143871, 0.02864804, 0.04706062, 0.06534370, 0.08271709, 0.09832237,
     0.11088212, 0.11768435, 0.11088212, 0.09832237, 0.08271709, 0.06534370,
     0.04706062, 0.02864804, 0.01143871},
    {0.00864721, 0.02192319, 0.03635533, 0.05094765, 0.06507559, 0.07835904,
     0.09016019, 0.09953920, 0.10451709, 0.09953920, 0.09016019, 0.07835904,
     0.06507559, 0.05094765, 0.03635533, 0.02192319, 0.00864721},
    {0.00677711, 0.01714213, 0.02875466, 0.04058941, 0.05231192, 0.06349388,
     0.07389996, 0.08307035, 0.09028023, 0.09404233, 0.09028023, 0.08307035,
     0.07389996, 0.06349388, 0.05231192, 0.04058941, 0.02875466, 0.01714213,
     0.00677711},
    {0.00539676, 0.01380479, 0.02324506, 0.03301030, 0.04274016, 0.05222798,
     0.06124167, 0.06956222, 0.07684212, 0.08251187, 0.08542734, 0.08251187,
     0.07684212, 0.06956222, 0.06124167, 0.05222798, 0.04274016, 0.03301030,
     0.02324506, 0.01380479, 0.00539676},
    {0.00436407, 0.01128483, 0.01907035, 0.02724371, 0.03548930, 0.04358317,
     0.05139063, 0.05877813, 0.06555273, 0.07144232, 0.07599151, 0.07830131,
     0.07599151, 0.07144232, 0.06555273, 0.05877813, 0.05139063, 0.04358317,
     0.03548930, 0.02724371, 0.01907035, 0.01128483, 0.00436407},
    {0.00357401, 0.00933387, 0.01591871, 0.02283374, 0.02984252,
     0.03679677, 0.04360144, 0.05011413, 0.05625260, 0.06185059,
     0.06668946, 0.07039971, 0.07226263, 0.07039971, 0.06668946,
     0.06185059, 0.05625260, 0.05011413, 0.04360144, 0.03679677,
     0.02984252, 0.02283374, 0.01591871, 0.00933387, 0.00357401},
    {0.00300967, 0.00785704, 0.01343839, 0.01934213, 0.02538384, 0.03141969,
     0.03734759, 0.04311305, 0.04861305, 0.05377237, 0.05845760, 0.06248655,
     0.06555555, 0.06708126, 0.06555555, 0.06248655, 0.05845760, 0.05377237,
     0.04861305, 0.04311305, 0.03734759, 0.03141969, 0.02538384, 0.01934213,
     0.01343839, 0.00785704, 0.00300967},
    {0.00255589, 0.00668959, 0.01144620, 0.01656341, 0.02180479, 0.02706700,
     0.03229665, 0.03739725, 0.04232875, 0.04702181, 0.05140609, 0.05537219,
     0.05876676, 0.06133728, 0.06260387, 0.06133728, 0.05876676, 0.05537219,
     0.05140609, 0.04702181, 0.04232875, 0.03739725, 0.03229665, 0.02706700,
     0.02180479, 0.01656341, 0.01144620, 0.00668959, 0.00255589},
    {0.00218643, 0.00573808, 0.00987911, 0.01432013, 0.01890347, 0.02353501,
     0.02814209, 0.03268025, 0.03710248, 0.04135805, 0.04539827, 0.04915900,
     0.05254896, 0.05543801, 0.05761402, 0.05867763, 0.05761402, 0.05543801,
     0.05254896, 0.04915900, 0.04539827, 0.04135805, 0.03710248, 0.03268025,
     0.02814209, 0.02353501, 0.01890347, 0.01432013, 0.00987911, 0.00573808,
     0.00218643}};

template <typename ImageType>
bool interpolatePatch(double x, double y, int window_half_size,
                      const cv::Mat &input, const cv::Mat &mask,
                      cv::Mat &b_vec) {
  if (x > window_half_size + 1 && x < input.cols - (window_half_size + 2) &&
      y > window_half_size + 1 && y < input.rows - (window_half_size + 2)) {
    int x0 = int(x);
    int y0 = int(y);
    double xw = x - x0;
    double yw = y - y0;

    // precompute bilinear interpolation weights
    double w00 = (1.0 - xw) * (1.0 - yw), w01 = xw * (1.0 - yw),
           w10 = (1.0 - xw) * yw, w11 = xw * yw;

    // fit to local neighborhood = b vector...
    const uint8_t *v = mask.ptr<const uint8_t>(0);
    double *m = b_vec.ptr<double>(0);
    double mn = DBL_MAX;
    double mx = DBL_MIN;
    for (int wy = -window_half_size; wy <= window_half_size; wy++) {
      const ImageType *im00 = input.ptr<ImageType>(y0 + wy),
                      *im10 = input.ptr<ImageType>(y0 + wy + 1);
      for (int wx = -window_half_size; wx <= window_half_size; wx++) {
        if (*v > 0) {
          const int col0 = x0 + wx;
          const int col1 = col0 + 1;
          double val = im00[col0] * w00 + im00[col1] * w01 + im10[col0] * w10 +
                       im10[col1] * w11;
          *(m++) = val;
          if (mn > val)
            mn = val;
          if (mx < val)
            mx = val;
        }
        v++;
      }
    }
    if (mx - mn > 1.0 / 255)
      return true;
  }
  return false;
}

template <typename SaddlePointType> struct PolynomialFit {
  int initSaddleFitting(int half_kernel_size);

  template <typename LocationsPointType, typename SmoothedImageType = double>
  void saddleSubpixelRefinement(const cv::Mat &smoothed_input,
                                const std::vector<LocationsPointType> &initial,
                                std::vector<SaddlePoint> &refined,
                                int max_iterations = 3,
                                bool tight_convergence = true) {
    cv::Mat b(invAtAAt.cols, 1, CV_64FC1);

    double convergence_region = window_half_size;
    if (tight_convergence)
      convergence_region = 1.0;

    refined.resize(initial.size());
    for (size_t idx = 0; idx < initial.size(); idx++)
      refined[idx] = SaddlePoint(initial[idx].x, initial[idx].y);

    for (size_t idx = 0; idx < refined.size(); idx++) {
      SaddlePoint &pt = refined[idx];
      bool point_diverged = true;
      for (int it = 0; it < max_iterations; it++) {
        if (interpolatePatch<SmoothedImageType>(pt.x, pt.y, window_half_size,
                                                smoothed_input, mask, b)) {
          // fit quadric to surface by solving LSQ
          cv::Mat p = invAtAAt * b;

          // k5, k4, k3, k2, k1, k0
          // 0 , 1 , 2 , 3 , 4 , 5
          double *r = p.ptr<double>(0);
          pt.det = 4.0 * r[0] * r[1] - r[2] * r[2]; // 4.0 * k5 * k4 - k3 * k3

          // check if it is still a saddle point
          if (pt.det > 0)
            break;

          // compute the new location
          double dx = (-2.0 * r[1] * r[4] + r[2] * r[3]) /
                      pt.det; // - 2 * k4 * k1 +     k3 * k2
          double dy = (r[2] * r[4] - 2.0 * r[0] * r[3]) /
                      pt.det; //       k3 * k1 - 2 * k5 * k2
          pt.x += dx;
          pt.y += dy;

          if (detector_params.spatial_convergence_threshold > fabs(dx) &&
              detector_params.spatial_convergence_threshold > fabs(dy)) {
            double k4mk5 = r[1] - r[0];
            pt.s = sqrt(r[2] * r[2] + k4mk5 * k4mk5);
            pt.a1 = atan2(-r[2], k4mk5) / 2.0;
            pt.a2 = acos((r[1] + r[0]) / pt.s) / 2.0;
            // converged
            point_diverged = false;
            break;
          }
          // check for divergence due to departure out of convergence region or
          // point type change
          if (fabs(pt.x - initial[idx].x) > convergence_region ||
              fabs(pt.y - initial[idx].y) > convergence_region)
            break;
        } else
          break;
      }
      if (point_diverged) {
        pt.x = pt.y = std::numeric_limits<double>::infinity();
        ++diverged;
      }
    }
  }

  template <typename MonkeySaddlePointType, typename LocationsPointType,
            typename SmoothedImageType = double>
  void saddleSubpixelRefinement(const cv::Mat &smoothed_input,
                                const std::vector<LocationsPointType> &initial,
                                std::vector<MonkeySaddlePointType> &refined,
                                int max_iterations = 10,
                                bool tight_convergence = true) {
    cv::Mat b(invAtAAt.cols, 1, CV_64FC1);

    refined.resize(initial.size());
    for (size_t idx = 0; idx < initial.size(); idx++)
      refined[idx] = MonkeySaddlePointType(initial[idx].x, initial[idx].y);

    double convergence_region = window_half_size;
    if (tight_convergence)
      convergence_region = 1.0;

    std::complex<double> roots[3];
    for (size_t idx = 0; idx < refined.size(); idx++) {
      MonkeySaddlePointType &pt = refined[idx];
      bool point_diverged = true;
      int divergence_reason = 3; //  max iterations reached..
      UNUSED(divergence_reason);
      for (int it = 0; it < max_iterations; it++) {
        if (interpolatePatch<SmoothedImageType>(pt.x, pt.y, window_half_size,
                                                smoothed_input, mask, b)) {
          // fit cubic to surface by solving LSQ
          cv::Mat p = invAtAAt * b;
          const double *a =
              p.ptr<double>(0) - 1; // 1 based indexing for convenience

          // now use second derivatives to find the location of critical point
          //
          // f_xx = 6 * a1 * dx + 2 * a2 * dy + 2 * a5 = 0
          // f_xy = 2 * a2 * dx + 2 * a3 * dy +     a6 = 0
          // f_yy = 2 * a3 * dx + 6 * a4 * dy + 2 * a7 = 0
          //
          // A  = [ 3*a1   a2 ]           b = [ -a5 ]
          //      [ 2*a2 2*a3 ]               [ -a6 ]
          //      [   a3 3*a4 ]               [ -a7 ]
          //
          // this is over determined system, solve it by: inv(A'A) = A'b
          //
          const double AtA_a = 9 * a[1] * a[1] + 4 * a[2] * a[2] + a[3] * a[3];
          const double AtA_b =
              3 * a[1] * a[2] + 4 * a[2] * a[3] + 3 * a[3] * a[4];
          const double AtA_d = a[2] * a[2] + 4 * a[3] * a[3] + 9 * a[4] * a[4];

          // inverse...
          const double det = AtA_a * AtA_d - AtA_b * AtA_b;

          const double Atb_a = -3 * a[1] * a[5] - 2 * a[2] * a[6] - a[3] * a[7];
          const double Atb_b = -a[2] * a[5] - 2 * a[3] * a[6] - 3 * a[4] * a[7];

          const double dx = (AtA_d * Atb_a - AtA_b * Atb_b) / det;
          const double dy = (-AtA_b * Atb_a + AtA_a * Atb_b) / det;

          // journal of physics tests for umbilic point types
          // J = 3 * (a1 * a3 + a2 * a4) - (a2 * a2 + a3 * a3)
          pt.det =
              3 * (a[1] * a[3] + a[2] * a[4]) - (a[2] * a[2] + a[3] * a[3]);

          // new location
          pt.x += dx;
          pt.y += dy;

// keep dx, dy for debugging in case the point diverges
#ifdef DEBUG_INDEXING
          pt.a3 = std::complex<double>(dx, dy);
#endif
          // check if it is still a monkey saddle point
          if (pt.det >= 0) {
            divergence_reason = 1; // does not have proper shape
            break;
          }

          if (detector_params.spatial_convergence_threshold > fabs(dx) &&
              detector_params.spatial_convergence_threshold > fabs(dy)) {
            // recover angles as roots of cubic equation (it assumes 0 indexed
            // array, thus pass a+1):
            solveCubicPolynomial(a + 1, roots);
            pt.a1 = atan(roots[0]);
            pt.a2 = atan(roots[1]);
            pt.a3 = atan(roots[2]);
            point_diverged = false;
            break;
          }
          // check for divergence due to departure out of convergence region
          if (fabs(pt.x - initial[idx].x) > convergence_region ||
              fabs(pt.y - initial[idx].y) > convergence_region) {
            divergence_reason = 2; // departed from convergence region
            break;
          }
        } else
          break;
      }
      if (point_diverged) {
#ifdef DEBUG_INDEXING
        // keep some divergence debugging info...
        if (divergence_reason == 3) {
          pt.a1 = pt.a3.real();
          pt.a2 = pt.a3.imag();
        } else {
          pt.a1 = pt.x;
          pt.a2 = pt.y;
        }
        pt.s = divergence_reason;
#endif
        pt.x = pt.y = std::numeric_limits<double>::infinity();
        ++diverged;
      }
    }
  }

  cv::Mat getSmoothingKernel() { return smoothingKernel; }

  int getHalfKernelSize() { return window_half_size; }

private:
  int window_half_size;
  int diverged = 0;
  int iterations = 0;
  cv::Mat smoothingKernel;
  cv::Mat invAtAAt;
  cv::Mat mask;

  int initConeSmoothingKernel() {
    int window_size = window_half_size * 2 + 1;
    smoothingKernel.create(window_size, window_size,
                           cv::DataType<float>::depth);
    mask.create(window_size, window_size, CV_8UC1);
    double maxVal = window_half_size + 1, sum = 0;
    float *w = smoothingKernel.ptr<float>(0);
    // cone kernel
    int nnz = 0;
    for (int y = -window_half_size; y <= window_half_size; y++)
      for (int x = -window_half_size; x <= window_half_size; x++) {
        *w = float(maxVal - sqrt(x * x + y * y));
        if (*w > 0)
          nnz++;
        else
          *w = 0;
        sum += *w;
        w++;
      }
    // scale kernel
    smoothingKernel /= sum;
    return nnz;
  }
};

template <>
int PolynomialFit<SaddlePoint>::initSaddleFitting(int half_kernel_size) {
  window_half_size = half_kernel_size;
  int nnz = initConeSmoothingKernel();
  cv::Mat A(nnz, 6, CV_64FC1);

  double *a = A.ptr<double>(0);
  float *w = smoothingKernel.ptr<float>(0);
  uint8_t *m = mask.ptr<uint8_t>(0);
  for (int y = -half_kernel_size; y <= half_kernel_size; y++)
    for (int x = -half_kernel_size; x <= half_kernel_size; x++) {
      if (*w > 0) {
        *m = 255;
        a[0] = x * x;
        a[1] = y * y;
        a[2] = x * y;
        a[3] = y;
        a[4] = x;
        a[5] = 1;
        a += 6;
      } else
        *m = 0;
      w++;
      m++;
    }
  // compute pseudoinverse of AtA
  cv::invert(A.t() * A, invAtAAt, cv::DECOMP_SVD);
  invAtAAt *= A.t();
  return nnz;
}

template <>
int PolynomialFit<MonkeySaddlePoint>::initSaddleFitting(int half_kernel_size) {
  window_half_size = half_kernel_size;
  int nnz = initConeSmoothingKernel();
  cv::Mat A(nnz, 10, CV_64FC1);
  double *a = A.ptr<double>(0);
  uint8_t *m = mask.ptr<uint8_t>(0);
  float *w = smoothingKernel.ptr<float>(0);

  for (int y = -half_kernel_size; y <= half_kernel_size; y++)
    for (int x = -half_kernel_size; x <= half_kernel_size; x++) {
      if (*w > 0) {
        *m = 255;
        a[0] = x * x * x;
        a[1] = x * x * y;
        a[2] = x * y * y;
        a[3] = y * y * y;
        a[4] = x * x;
        a[5] = x * y;
        a[6] = y * y;
        a[7] = x;
        a[8] = y;
        a[9] = 1;
        a += 10;
      } else
        *m = 0;
      w++;
      m++;
    }
  // compute pseudoinverse of AtA
  cv::invert(A.t() * A, invAtAAt, cv::DECOMP_SVD);
  invAtAAt *= A.t();
  return nnz;
}

template <>
int PolynomialFit<MonkeySaddlePointSpherical>::initSaddleFitting(
    int half_kernel_size) {
  window_half_size = half_kernel_size;
  int nnz = initConeSmoothingKernel();
  cv::Mat A(nnz, 10, CV_64FC1);
  double *a = A.ptr<double>(0);
  uint8_t *m = mask.ptr<uint8_t>(0);
  float *w = smoothingKernel.ptr<float>(0);

  for (int y = -half_kernel_size; y <= half_kernel_size; y++)
    for (int x = -half_kernel_size; x <= half_kernel_size; x++) {
      if (*w > 0) {
        *m = 255;
        a[0] = x * x * x;
        a[1] = x * x * y;
        a[2] = x * y * y;
        a[3] = y * y * y;
        a[4] = x * x;
        a[5] = x * y;
        a[6] = y * y;
        a[7] = x;
        a[8] = y;
        a[9] = 1;
        a += 10;
      } else
        *m = 0;
      w++;
      m++;
    }
  // compute pseudoinverse of AtA
  cv::invert(A.t() * A, invAtAAt, cv::DECOMP_SVD);
  invAtAAt *= A.t();
  return nnz;
}
}
}

#endif // INCLUDE_CHECKERBOARD_DETECTOR_POLYNOMIALFIT_H_
