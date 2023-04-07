#include <deltille/PolynomialFit.h>
namespace orp {
namespace calibration {

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

} // namespace orp
} // namespace calibration
