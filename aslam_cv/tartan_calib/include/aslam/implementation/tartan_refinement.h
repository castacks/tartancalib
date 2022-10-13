// Copyright 2019 ETH Zürich, Thomas Schöps
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include "libvis/eigen.h"
#include "libvis/image.h"
// #include <libvis/libvis.h>

namespace vis {


template <typename T, typename Derived, typename C>
bool RefineFeatureByTartan(
    int num_samples,
    vector<vector<Eigen::VectorXd>> samples_pointcloud,
    const Image<T>& image,
    int window_half_size,
    const MatrixBase<Derived>& position,
    const Mat4d& local_pattern_tr_pixel,
    const Eigen::VectorXd corner_3d,
    Vec2f* out_position,
    float* final_cost,
    bool debug,
    C* original_camera,
    cv::Mat cv_image
    ) {
    constexpr int kDim = 6;
    Matrix<float, kDim, kDim> H;
    Matrix<float, kDim, 1> b;
    float lambda = -1;

    Eigen::MatrixXd original_cam_params, final_params;
    original_camera->getParameters(original_cam_params,true,false,false);
    const Eigen::MatrixXd  original_cam_params_const = original_cam_params;
    final_params = original_cam_params;


    // SM_INFO_STREAM("Cost "<<init_cost);
    // Initialize lambda?

    bool applied_update = false;
    for (int meta_iteration = 0; meta_iteration < 30; ++ meta_iteration) {
      float init_cost = 0;
      original_camera->setParameters(final_params,true,false,false);

      ComputeCostandJacobian(
      num_samples,
      samples_pointcloud,
      image,
      window_half_size,
      position,
      local_pattern_tr_pixel,
      corner_3d,
      out_position,
      &init_cost,
      debug,
      original_camera,
      cv_image,
      &H,
      &b
      );

      lambda = -1;
      if (lambda < 0) {
      lambda = 0.001f * (1.f / kDim) * H.diagonal().sum();
      }

      for (int lm_iteration = 0; lm_iteration < 2; ++lm_iteration)
      {
        Matrix<float, kDim, kDim> H_LM;
        H_LM.triangularView<Eigen::Upper>() = H.triangularView<Eigen::Upper>();
        // SM_INFO_STREAM("Matrix: "<<H);
        // H_LM.diagonal().array() += lambda;
        
        // Solve for the update.
        // NOTE: Not sure if using double is helpful here
        
        Eigen::Matrix<double, kDim, 1> x = H_LM.cast<double>().selfadjointView<Eigen::Upper>().ldlt().solve(b.cast<double>()).cast<double>();
        // SM_INFO_STREAM("Suggested change: "<<x);
        // cv::waitKey(500);
        // SM_INFO_STREAM("Suggested changes: "<<x);
        const Eigen::MatrixXd& test_params = final_params-x;
        original_camera->setParameters(test_params,true,false,false);   

        float test_cost = 0;
        ComputeCost(
        num_samples,
        samples_pointcloud,
        image,
        window_half_size,
        position,
        local_pattern_tr_pixel,
        corner_3d,
        out_position,
        &test_cost,
        debug,
        original_camera,
        cv_image,
        &H,
        &b
        ); 
        if (test_cost < init_cost)
        {
          // init_cost = test_cost;
          final_params = test_params;
          // lambda *= 0.9f;
          break;
        }
        else
        {
          // lambda *= 1.1f;
        }
        original_camera->setParameters(final_params,true,false,false);

      }
    }
    // get the final predicted point
    Eigen::VectorXd out_position_eig;
    original_camera->setParameters(final_params,true,false,false);
    original_camera->vsEuclideanToKeypoint(corner_3d,out_position_eig);
    *out_position = Vec2f(out_position_eig(0),out_position_eig(1));
    

    // setting parameters back to the initial value
    original_camera->setParameters(original_cam_params_const,true,false,false);   

      // ComputeCostandJacobian(
      // num_samples,
      // samples_pointcloud,
      // image,
      // window_half_size,
      // position,
      // local_pattern_tr_pixel,
      // corner_3d,
      // out_position,
      // &out_cost,
      // debug,
      // original_camera,
      // cv_image,
      // &H,
      // &b
      // );
      // const Eigen::MatrixXd& updated_cam_params = cam_params - x.cast<double>() ;
      // local_camera->setParameters(updated_cam_params,true,false,false);
      // SM_INFO_STREAM("ORIGINAL PARAMETERS: " << cam_params);
      // debugging

    // local_camera->vsEuclideanToKeypoint(corner_3d,corner_imageframe);
    // int red = static_cast<int>((1-static_cast<float>(it)/static_cast<float>(num_iterations))*255.);
    // int green =  static_cast<int>((static_cast<float>(it)/static_cast<float>(num_iterations))*255.);
    // cv::circle(cv_image, cv::Point2f(corner_imageframe(0),corner_imageframe(1)),0, cv::Scalar(0,green,red),2);           
    // cv::circle(cv_image, cv::Point2f(negative_sample_image_frame(0),negative_sample_image_frame(1)),0, cv::Scalar(0,255,0),2);    
    // cv::imwrite("test_"+std::to_string(it)+".png",cv_image);



//       // TEST: use pseudoinverse to ensure Gauge fixing
//       H_LM.template triangularView<Eigen::Lower>() = H_LM.template triangularView<Eigen::Upper>().transpose();
//       Eigen::Matrix<float, kDim, 1> x =
//           H_LM.cast<double>().completeOrthogonalDecomposition().solve(b.cast<double>()).cast<float>();
    //   if (kDebug) {
    //     // LOG(INFO) << "  x in LM iteration " << lm_iteration << ": " << x.transpose();
    //   }
      
    //   // Test whether the update improves the cost.
    //   Mat3f test_pixel_tr_pattern_samples;
    //   test_pixel_tr_pattern_samples(0, 0) = pixel_tr_pattern_samples(0, 0) - x(0);
    //   test_pixel_tr_pattern_samples(0, 1) = pixel_tr_pattern_samples(0, 1) - x(1);
    //   test_pixel_tr_pattern_samples(0, 2) = pixel_tr_pattern_samples(0, 2) - x(2);
    //   test_pixel_tr_pattern_samples(1, 0) = pixel_tr_pattern_samples(1, 0) - x(3);
    //   test_pixel_tr_pattern_samples(1, 1) = pixel_tr_pattern_samples(1, 1) - x(4);
    //   test_pixel_tr_pattern_samples(1, 2) = pixel_tr_pattern_samples(1, 2) - x(5);
    //   test_pixel_tr_pattern_samples(2, 0) = pixel_tr_pattern_samples(2, 0) - x(6);
    //   test_pixel_tr_pattern_samples(2, 1) = pixel_tr_pattern_samples(2, 1) - x(7);
    //   test_pixel_tr_pattern_samples(2, 2) = pixel_tr_pattern_samples(2, 2);
    //   float test_cost;
    //   if (!CostFunction::ComputeCornerRefinementCost(test_pixel_tr_pattern_samples, image, num_samples, pattern_samples, &test_cost)) {
    //     if (kDebug) {
    //     //   LOG(INFO) << "  CostFunction::ComputeCornerRefinementCost() failed, aborting.";
    //     }
    //     return false;
    //   }
      
    //   if (kDebug) {
    //     // LOG(INFO) << "  test_cost: " << test_cost << ", cost: " << cost;
    //   }
      
    //   if (test_cost < cost) {
    //     if (final_cost) {
    //       *final_cost = test_cost;
    //     }
    //     last_step_squared_norm = Vec2f(x(2), x(5)).squaredNorm();  // using the translation only
    //     pixel_tr_pattern_samples = test_pixel_tr_pattern_samples;
    //     lambda *= 0.5f;
    //     applied_update = true;
    //     break;
    //   } else {
    //     lambda *= 2.f;
    //   }
    // }
    
    // if (applied_update) {
    //   // Since the element at (2, 2) is always 1, we can directly assign the
    //   // translation values instead of computing:
    //   // *out_position = (pixel_tr_pattern_samples * Vec3f(0, 0, 1)).hnormalized();
    //   out_position->x() = pixel_tr_pattern_samples(0, 2);
    //   out_position->y() = pixel_tr_pattern_samples(1, 2);
    //   if (kDebug) {
    //     // LOG(INFO) << "New position: " << out_position->transpose();
    //   }
    // } else {
    //   // Cannot find an update that improves the cost. Treat this as converged.
    //   if (kDebug) {
    //     // LOG(INFO) << "Cannot find an update to improve the cost. Returning convergence (iteration " << iteration << ").";
    //   }
    //   return true;
    // }
    
    // // Check for divergence.
    // if (fabs(original_position.x() - out_position->x()) >= window_half_size ||
    //     fabs(original_position.y() - out_position->y()) >= window_half_size) {
    //   // The result is probably not the originally intended corner,
    //   // since it is not within the original search window.
    //   if (debug || kDebug) {
    //     // LOG(WARNING) << "Corner refinement failed because the refined position left the original window";
    //   }
    //   return false;


    // TODO: Why was this commented out? For parity with the CUDA version?
//     // Check for convergence.
//     if (x.squaredNorm() < numeric_limits<float>::epsilon()) {
//       return true;
//     }
  }








// The positions are specified in pixel-center origin convention.
template <typename T, typename Derived, typename C>
bool ComputeCostandJacobian(
    int num_samples,
    vector<vector<Eigen::VectorXd>> samples_pointcloud,
    const Image<T>& image,
    int window_half_size,
    const MatrixBase<Derived>& position,
    const Mat4d& local_pattern_tr_pixel,
    const Eigen::VectorXd corner_3d,
    Vec2f* out_position,
    float* out_cost,
    bool debug,
    C* original_camera,
    cv::Mat cv_image,
    Matrix<float, 6, 6>* H,
    Matrix<float, 6, 1>* b
    ) {
      H->triangularView<Eigen::Upper>().setZero();
      b->setZero();
      // doing a deep copy with the original cam, we can't change the original camera
      C* local_camera = original_camera;
      *out_cost = 0;
      // SM_INFO_STREAM("1");
      // Eigen::MatrixXd original_cam_params;
      // original_camera->getParameters(original_cam_params,true,true,true);
      // SM_INFO_STREAM("Original_params: "<< original_cam_params);
      // local_camera->setParameters(original_cam_params,true,true,true);
      // SM_INFO_STREAM("1");
      int num_iterations = 2;
      double learning_rate = 0.000000001;
      Eigen::MatrixXd original_cam_params;
      original_camera->getParameters(original_cam_params,true,true,true);
      int num_intr_params = original_cam_params.cols();
      Eigen::VectorXd corner_imageframe;
      double cost = 0;

      // cv::imwrite("grad_debug.png",grad_x);
      // SM_INFO_STREAM("Wrote grad debug image");
      // cv::waitKey(20000);

      // // debug effort
      // cv::Mat grad_x(cv::Size(1224,1028),CV_32FC1);
      // float intensity;
      // Eigen::Matrix<float, 1, 2> gradient;
      // Vec2f sample_coord;
      // for (int i = 0; i< 1028; i++)
      // {
      //   SM_INFO_STREAM("ROW "<<i);
      //   for (int j = 0; j < 1224; j++)
      //   {
      //       sample_coord = Vec2f(i,j);
      //       image.InterpolateBilinearWithJacobian(sample_coord, &intensity, &gradient);
      //       grad_x.at<float>(i,j) = gradient.coeff(0,0);
      //   }
      // }
      // cv::imwrite("grad_debug.png",grad_x);
      // cv::waitKey(20000);


        Eigen::Matrix<float,1,6> jacobian_intrinsics = Eigen::MatrixXf::Zero(1,6); 
        // step 1: create the first pointcloud
        for (int num_sample = 0; num_sample < num_samples; num_sample++)
        {
          Eigen::VectorXd positive_sample_image_frame, negative_sample_image_frame;
          Eigen::MatrixXd positive_outJacobian_projection, negative_outJacobian_projection;
          vector<Eigen::VectorXd> sample_couple = samples_pointcloud[num_sample]; 
          // getting jacobian of keypoint w.r.t. intrinsics
          local_camera->euclideanToKeypointIntrinsicsJacobian(sample_couple[0],positive_outJacobian_projection,true,false,false);
          local_camera->euclideanToKeypointIntrinsicsJacobian(sample_couple[1],negative_outJacobian_projection,true,false,false);
          // reproject 3D points into image plane
          local_camera->vsEuclideanToKeypoint(sample_couple[0],positive_sample_image_frame);
          local_camera->vsEuclideanToKeypoint(sample_couple[1],negative_sample_image_frame);

          Vec2f sample_pos = Vec2f(positive_sample_image_frame(1),positive_sample_image_frame(0));
          Vec2f sample_neg = Vec2f(negative_sample_image_frame(1),negative_sample_image_frame(0));
          bool in_image = (image.ContainsPixelCenterConv(sample_pos) && image.ContainsPixelCenterConv(sample_neg));

          if (in_image)
          {
            // cv::circle(cv_image, cv::Point2f(positive_sample_image_frame(0),positive_sample_image_frame(1)),0, cv::Scalar(255,0,0),2); 
            // cv::circle(cv_image, cv::Point2f(negative_sample_image_frame(0),negative_sample_image_frame(1)),0, cv::Scalar(255,0,0),2); 

            // gradient of positive sample
            float intensity_pos;
            Eigen::Matrix<float, 1, 2> gradient_pos;
            image.InterpolateBilinearWithJacobian(sample_pos, &intensity_pos, &gradient_pos);
            float grad_y = gradient_pos.coeff(0,0);
            gradient_pos(0,0) = gradient_pos(0,1);
            gradient_pos(0,1) = grad_y;

            // gradient of negative sample
            float intensity_neg;
            Eigen::Matrix<float, 1, 2> gradient_neg;
            image.InterpolateBilinearWithJacobian(sample_neg, &intensity_neg, &gradient_neg);
            grad_y = gradient_neg.coeff(0,0);
            gradient_neg(0,0) = gradient_neg(0,1);
            gradient_neg(0,1) = grad_y;

            double residual = static_cast<double>(intensity_pos - intensity_neg);
            *out_cost += static_cast<float>(residual*residual);

            jacobian_intrinsics += gradient_pos*(positive_outJacobian_projection.cast<float>());
            jacobian_intrinsics -= gradient_neg*(negative_outJacobian_projection.cast<float>());

             // Accumulate update equation coefficients.
            float weight = 1.0;
            Matrix<float, 1, 6> jacobian_weighted = weight * jacobian_intrinsics;
            H->triangularView<Eigen::Upper>() += jacobian_intrinsics.transpose() * jacobian_weighted;
            *b += static_cast<float>(residual) * jacobian_weighted;

            // SM_INFO_STREAM("Position pos: "<<positive_sample_image_frame);
            // SM_INFO_STREAM("Negative pos: "<<negative_sample_image_frame);
            // SM_INFO_STREAM("Raw Jacobian positive :"<<positive_outJacobian_projection);
            // SM_INFO_STREAM("Img jac positive :"<<gradient_pos);

            // SM_INFO_STREAM("Raw Jacobian negative :"<<negative_outJacobian_projection);
            // SM_INFO_STREAM("Img jac negative :"<<gradient_neg);
            // cv::waitKey(1000);
          }
        // *out_cost = static_cast<float>(cost);
        // Eigen::MatrixXd cam_params;
        // local_camera->getParameters(cam_params,true,false,false);
        

        // Eigen::MatrixXd jacobian_corrected = learning_rate*jacobian_intrinsics.cast<double>();
        // // SM_INFO_STREAM("Jacobian : "<<jacobian_corrected);
        // // SM_INFO_STREAM("Cost: "<<cost);
        // const Eigen::MatrixXd& updated_cam_params = cam_params - jacobian_corrected ;
        // local_camera->setParameters(updated_cam_params,true,false,false);
        // // SM_INFO_STREAM("ORIGINAL PARAMETERS: " << cam_params);
        // // debugging

        // local_camera->vsEuclideanToKeypoint(corner_3d,corner_imageframe);
        // int red = static_cast<int>((1-static_cast<float>(it)/static_cast<float>(num_iterations))*255.);
        // int green =  static_cast<int>((static_cast<float>(it)/static_cast<float>(num_iterations))*255.);
        // cv::circle(cv_image, cv::Point2f(corner_imageframe(0),corner_imageframe(1)),0, cv::Scalar(0,green,red),2);           
        // cv::circle(cv_image, cv::Point2f(negative_sample_image_frame(0),negative_sample_image_frame(1)),0, cv::Scalar(0,255,0),2);    
        // cv::imwrite("test_"+std::to_string(it)+".png",cv_image);
      }
      // SM_INFO_STREAM("Jacobian: "<<jacobian_intrinsics);
    
    }


// The positions are specified in pixel-center origin convention.
template <typename T, typename Derived, typename C>
bool ComputeCost(
    int num_samples,
    vector<vector<Eigen::VectorXd>> samples_pointcloud,
    const Image<T>& image,
    int window_half_size,
    const MatrixBase<Derived>& position,
    const Mat4d& local_pattern_tr_pixel,
    const Eigen::VectorXd corner_3d,
    Vec2f* out_position,
    float* out_cost,
    bool debug,
    C* original_camera,
    cv::Mat cv_image,
    Matrix<float, 6, 6>* H,
    Matrix<float, 6, 1>* b
    ) {
      // doing a deep copy with the original cam, we can't change the original camera
      C* local_camera = original_camera;
      *out_cost = 0;
      // SM_INFO_STREAM("1");
      // Eigen::MatrixXd original_cam_params;
      // original_camera->getParameters(original_cam_params,true,true,true);
      // SM_INFO_STREAM("Original_params: "<< original_cam_params);
      // local_camera->setParameters(original_cam_params,true,true,true);
      // SM_INFO_STREAM("1");
      Eigen::MatrixXd original_cam_params;
      original_camera->getParameters(original_cam_params,true,true,true);
      int num_intr_params = original_cam_params.cols();
      Eigen::VectorXd corner_imageframe;

      // //debug effort
      // cv::Mat grad_x(cv::Size(1028,1224),CV_32FC1);
      // float intensity;
      // Eigen::Matrix<float, 1, 2> gradient;
      // Vec2f sample_coord;
      // for (int i = 0; i< 1028; i++)
      // {
      //   SM_INFO_STREAM("ROW "<<i);
      //   for (int j = 0; j < 1224; j++)
      //   {
      //       sample_coord = Vec2f(i,j);
      //       image.InterpolateBilinearWithJacobian(sample_coord, &intensity, &gradient);
      //       grad_x.at<float>(i,j) = gradient.coeff(0,1);
      //   }
      // }
      // cv::imwrite("grad_debug.png",grad_x);

        // step 1: create the first pointcloud
        for (int num_sample = 0; num_sample < num_samples; num_sample++)
        {
          Eigen::VectorXd positive_sample_image_frame, negative_sample_image_frame;
          Eigen::MatrixXd positive_outJacobian_projection, negative_outJacobian_projection;
          vector<Eigen::VectorXd> sample_couple = samples_pointcloud[num_sample]; 
          // getting jacobian of keypoint w.r.t. intrinsics
          // local_camera->euclideanToKeypointIntrinsicsJacobian(sample_couple[0],positive_outJacobian_projection,true,false,false);
          // local_camera->euclideanToKeypointIntrinsicsJacobian(sample_couple[1],negative_outJacobian_projection,true,false,false);
        
          // reproject 3D points into image plane
          local_camera->vsEuclideanToKeypoint(sample_couple[0],positive_sample_image_frame);
          local_camera->vsEuclideanToKeypoint(sample_couple[1],negative_sample_image_frame);

          Vec2f sample_pos = Vec2f(positive_sample_image_frame(1),positive_sample_image_frame(0));
          Vec2f sample_neg = Vec2f(negative_sample_image_frame(1),negative_sample_image_frame(0));
          bool in_image = (image.ContainsPixelCenterConv(sample_pos) && image.ContainsPixelCenterConv(sample_neg));

          if (in_image)
          {
            // cv::circle(cv_image, cv::Point2f(positive_sample_image_frame(0),positive_sample_image_frame(1)),0, cv::Scalar(255,0,0),2); 
            // cv::circle(cv_image, cv::Point2f(negative_sample_image_frame(0),negative_sample_image_frame(1)),0, cv::Scalar(255,0,0),2); 

            // gradient of positive sample
            float intensity_pos;
            Eigen::Matrix<float, 1, 2> gradient_pos;
            image.InterpolateBilinearWithJacobian(sample_pos, &intensity_pos, &gradient_pos);
            // float grad_y = gradient_pos.coeff(0,0);
            // gradient_pos(0,0) = gradient_pos(0,1);
            // gradient_pos(0,1) = grad_y;

            // gradient of negative sample
            float intensity_neg;
            Eigen::Matrix<float, 1, 2> gradient_neg;
            image.InterpolateBilinearWithJacobian(sample_neg, &intensity_neg, &gradient_neg);
            // grad_y = gradient_neg.coeff(0,0);
            // gradient_neg(0,0) = gradient_neg(0,1);
            // gradient_neg(0,1) = grad_y;

            double residual = static_cast<double>(intensity_pos - intensity_neg);
            *out_cost += static_cast<float>(residual*residual);

            // SM_INFO_STREAM("Raw Jacobian positive :"<<positive_outJacobian_projection);
            // SM_INFO_STREAM("Img jac positive :"<<gradient_pos);

            // SM_INFO_STREAM("Raw Jacobian positive :"<<negative_outJacobian_projection);
            // SM_INFO_STREAM("Img jac positive :"<<gradien_neg);
          }
          else
          {
            *out_cost = numeric_limits<float>::infinity();
          }
        // *out_cost = static_cast<float>(cost);
        // Eigen::MatrixXd cam_params;
        // local_camera->getParameters(cam_params,true,false,false);
        

        // Eigen::MatrixXd jacobian_corrected = learning_rate*jacobian_intrinsics.cast<double>();
        // // SM_INFO_STREAM("Jacobian : "<<jacobian_corrected);
        // // SM_INFO_STREAM("Cost: "<<cost);
        // const Eigen::MatrixXd& updated_cam_params = cam_params - jacobian_corrected ;
        // local_camera->setParameters(updated_cam_params,true,false,false);
        // // SM_INFO_STREAM("ORIGINAL PARAMETERS: " << cam_params);
        // // debugging

        local_camera->vsEuclideanToKeypoint(corner_3d,corner_imageframe);
        // int red = static_cast<int>((1-static_cast<float>(it)/static_cast<float>(num_iterations))*255.);
        // int green =  static_cast<int>((static_cast<float>(it)/static_cast<float>(num_iterations))*255.);
        cv::circle(cv_image, cv::Point2f(corner_imageframe(0),corner_imageframe(1)),0, cv::Scalar(255,0,0),2);           
        // cv::circle(cv_image, cv::Point2f(negative_sample_image_frame(0),negative_sample_image_frame(1)),0, cv::Scalar(0,255,0),2);    
        // cv::imwrite("test_"+std::to_string(it)+".png",cv_image);
      }
     
    
    }


}