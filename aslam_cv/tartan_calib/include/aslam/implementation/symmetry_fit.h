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

template <typename C>
Vec2f TargetToImageFrame(
    Eigen::Vector4d target_location,
    Eigen::Matrix4d T_target_to_euc,
    C* camera)
{
  Eigen::VectorXd distorted_pixel_location;
  Eigen::VectorXd sample_3D = T_target_to_euc*target_location; 

  Eigen::Vector3d sample_input;
  sample_input(0) = sample_3D(0);
  sample_input(1) = sample_3D(1);
  sample_input(2) = sample_3D(2);

  camera->vsEuclideanToKeypoint(sample_input,distorted_pixel_location);
  return Vec2f(distorted_pixel_location(0),distorted_pixel_location(1));
}

template <typename C>
bool IsValid(
    Eigen::Vector4d target_location,
    Eigen::Matrix4d T_target_to_euc,
    C* camera,
    cv::Mat cv_img,
    double threshold
    )
{
  Eigen::VectorXd distorted_pixel_location;
  Eigen::VectorXd sample_3D = T_target_to_euc*target_location; 

  Eigen::Vector3d sample_input;
  sample_input(0) = sample_3D(0);
  sample_input(1) = sample_3D(1);
  sample_input(2) = sample_3D(2);

  camera->vsEuclideanToKeypoint(sample_input,distorted_pixel_location);
  // return Vec2f(distorted_pixel_location(0),distorted_pixel_location(1));
  bool valid = false;
  if (distorted_pixel_location(0) > threshold && distorted_pixel_location(0)< (cv_img.cols-threshold))
  {
    if (distorted_pixel_location(1) > threshold && distorted_pixel_location(1)< (cv_img.rows-threshold))
    {
      valid = true; 
    }
  }
  return valid;
}

template <typename T, typename C, size_t rows, size_t cols>
bool DebugScreen(
    Eigen::Vector4d start_target_frame,
    int num_samples_symmetry,
    int num_meta_samples_axis,
    vector<vector<Eigen::VectorXd>> samples_targetframe,
    Eigen::Vector4d (&meta_locations)[rows][cols],
    const Image<T>& image,
    Vec2f* out_position,
    float* final_cost,
    C* camera,
    cv::Mat cv_image,
    Eigen::Matrix4d T_target_to_euc,
    double* mean_sym
    )
    {
      cv::Mat intensity_img(cv::Size(num_meta_samples_axis,num_meta_samples_axis),CV_32FC1);
      for (int i = 0; i< num_meta_samples_axis; i++)
      {
        for (int j = 0; j< num_meta_samples_axis; j++)
        {
          intensity_img.at<float>(i,j) = static_cast<float>(TargetToIntensity(image, T_target_to_euc, camera, start_target_frame +  meta_locations[i][j],1.0));
        }
      }
      cv::imwrite("intensity_test.png",intensity_img);

      cv::Mat symmetry_img(cv::Size(num_meta_samples_axis,num_meta_samples_axis),CV_32FC1);

      for (int i = 0; i< num_meta_samples_axis; i++)
      {
        for (int j = 0; j< num_meta_samples_axis; j++)
        {
          Eigen::Vector2d Jacobian;
          Jacobian.setZero();

          symmetry_img.at<float>(i,j) = static_cast<float>(TargetToSymmetry(image, T_target_to_euc, camera, start_target_frame +  meta_locations[i][j],samples_targetframe,num_samples_symmetry,cv_image,1.0));
          // symmetry_img.at<float>(i,j) = static_cast<float>(TargetToJacobian(image, T_target_to_euc, camera, start_target_frame +  meta_locations[i][j],samples_targetframe,num_samples_symmetry,cv_image,&Jacobian));

        }
      }

      double minVal; 
      double maxVal; 
      cv::Point minLoc; 
      cv::Point maxLoc;

      minMaxLoc( symmetry_img, &minVal, &maxVal, &minLoc, &maxLoc );
      symmetry_img *= 1.0/maxVal;
      symmetry_img = 1.0 - symmetry_img;
      cv::Scalar tempVal = cv::mean( symmetry_img );
      *mean_sym = static_cast<double>(tempVal.val[0]);

      // cv::waitKey(20000);
      cv::imwrite("symmetry_test.png",symmetry_img*255);
      // SM_INFO_STREAM("Best target: "<<start_target_frame +  meta_locations[minLoc.x][minLoc.y]);
      *out_position = TargetToImageFrame(start_target_frame +  meta_locations[minLoc.x][minLoc.y], T_target_to_euc, camera);
    }  

template <typename T, typename C, size_t rows, size_t cols>
bool FitSymmetry(
    Eigen::Vector4d& start_target_frame,
    int num_samples_symmetry,
    int num_meta_samples_axis,
    vector<vector<Eigen::VectorXd>> samples_targetframe,
    Eigen::Vector4d (&meta_locations)[rows][cols],
    const Image<T>& image,
    Vec2f* out_position,
    float* final_cost,
    C* camera,
    cv::Mat cv_image,
    Eigen::Matrix4d T_target_to_euc,
    double scaling_factor
    ) {
      // SM_INFO_STREAM("Target initial position: "<<start_target_frame);
      Eigen::Vector4d intermediate_target_frame = start_target_frame;
      Eigen::Vector4d test_target_frame;
      // DebugScreen(start_target_frame, num_samples_symmetry,num_meta_samples_axis,samples_targetframe,meta_locations,image,out_position,final_cost,camera,cv_image,T_target_to_euc);
      constexpr int kDim = 2;
      Matrix<double, kDim, kDim> H;
      Matrix<double, kDim, 1> b;
      double lambda = -1;
      double last_step_squared_norm = numeric_limits<float>::infinity();
      Eigen::Vector2d Jacobian;
      // constexpr int num_it = 10;
      // for (int i = 0; i< num_it; i++)
      // {
      //   Eigen::Vector2d Jacobian;
      //   double cost = TargetToJacobian(image, T_target_to_euc, camera, intermediate_target_frame,samples_targetframe,num_samples_symmetry,cv_image,&Jacobian,&H, &b, scaling_factor);
      //   double jac_norm = Jacobian.norm();
      //   Eigen::Vector2d Jac_norm = Jacobian/jac_norm * 0.0005;

      //   intermediate_target_frame(0) -= Jac_norm(0);
      //   intermediate_target_frame(1) -= Jac_norm(1);
      //   // SM_INFO_STREAM("Refined target position: "<<out_pos);
      // }
      

      // *out_position = TargetToImageFrame(intermediate_target_frame,T_target_to_euc,camera);

      constexpr int kMaxIterationCount = 30;
      double norm_threshold = 0.005;

      // if (IsValid(start_target_frame,T_target_to_euc,camera,cv_image,10.0))
      if(true)
      {

        for (int iteration = 0; iteration < kMaxIterationCount; ++ iteration) {
          // lambda = -1;

        
        //   // *out_position = TargetToImageFrame(start_target_frame, T_target_to_euc, camera);

        //   // if (kDebug) {
        //   //   SM_INFO_STREAM("cost: " << cost);
        //   // }
        //   // if (final_cost) {
        //   //   *final_cost = cost;
        //   // }
      
        // // Initialize lambda?
        // SM_INFO_STREAM("Meta iteration: "<<iteration);
        if (lambda < 0) {
          // lambda = 0.001f * (1.f / kDim) * H.diagonal().sum();
          lambda = 1000;
        }

      
        bool applied_update = false;
        for (int lm_iteration = 0; lm_iteration < 10; ++ lm_iteration) {
          // SM_INFO_STREAM("Minor iteration: "<<lm_iteration);

          // this fill in H and b with the jacobian parts
          double cost = TargetToJacobian(image, T_target_to_euc, camera, intermediate_target_frame,samples_targetframe,num_samples_symmetry,cv_image,&Jacobian,&H, &b,scaling_factor);

          Matrix<double, kDim, kDim> H_LM;
          H_LM.triangularView<Eigen::Upper>() = H.triangularView<Eigen::Upper>();
          H_LM.diagonal().array() += lambda;
        

          Eigen::Matrix<double, kDim, 1> x = H_LM.cast<double>().selfadjointView<Eigen::Upper>().ldlt().solve(b.cast<double>());

          if (x.norm()>norm_threshold && x.norm() >0.000000001)
          {
            x *= (norm_threshold/x.norm());
          }
          

          test_target_frame = intermediate_target_frame;
          test_target_frame(0) -= x(0,0);
          test_target_frame(1) -= x(1,0);

          // if (IsValid(test_target_frame,T_target_to_euc,camera,cv_image,10.0))
          if (true)
          {
            double test_cost = TargetToSymmetry(image, T_target_to_euc, camera, test_target_frame,samples_targetframe,num_samples_symmetry,cv_image,scaling_factor);

            if (test_cost < cost) {
              if (final_cost) {
                *final_cost = test_cost;
              }
              intermediate_target_frame = test_target_frame;
              lambda *= 0.5f;
              applied_update = true;
              break;
            } else {
              lambda *= 2.f;
            }

          }
        

      }
      }
    }
   *out_position = TargetToImageFrame(intermediate_target_frame,T_target_to_euc,camera);
   start_target_frame = intermediate_target_frame;
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
    // }
  
  }




template <typename T, typename C>
double TargetToIntensity(
    const Image<T>& image,
    Eigen::Matrix4d T_target_to_euc,
    C* camera,
    Eigen::VectorXd target_location,
    double scaling_factor
    )
    {
      Eigen::VectorXd distorted_pixel_location;
      Eigen::VectorXd target_location_H = target_location;
      // target_location_H  << 0.0 , 1.0;
      Eigen::VectorXd sample_3D = T_target_to_euc*target_location_H;
      camera->vsEuclideanToKeypoint(sample_3D,distorted_pixel_location);
      Vec2f sample_pos = Vec2f(distorted_pixel_location(1),distorted_pixel_location(0));    
      return image.InterpolateBilinear(scaling_factor*sample_pos);
    }



template <typename T, typename C>
double TargetToSymmetry(
    const Image<T>& image,
    Eigen::Matrix4d T_target_to_euc,
    C* camera,
    Eigen::VectorXd target_location,
    vector<vector<Eigen::VectorXd>> samples_targetframe,
    int num_samples_symmetry,
    cv::Mat cv_img,
    double scaling_factor
    )
    {
      double C_symmetry = 0;
      double intensity_pos, intensity_neg;
      Eigen::VectorXd distorted_pixel_location;
      Eigen::VectorXd target_location_pos, target_location_neg, pos_3D, neg_3D, image_pos, image_neg;
      for (int i = 0; i < num_samples_symmetry ; i++)
      {
        vector<Eigen::VectorXd> sample_pair = samples_targetframe[i];
        
        target_location_pos = target_location;
        target_location_pos(0) += sample_pair[0](0);
        target_location_pos(1) += sample_pair[0](1);
        
        target_location_neg = target_location;
        target_location_neg(0) += sample_pair[1](0);
        target_location_neg(1) += sample_pair[1](1);
        
        pos_3D =  T_target_to_euc*target_location_pos;
        neg_3D =  T_target_to_euc*target_location_neg;

        Eigen::Vector3d pos_3D_input, neg_3D_input;

        pos_3D_input(0) = pos_3D(0);
        pos_3D_input(1) = pos_3D(1);
        pos_3D_input(2) = pos_3D(2);

        neg_3D_input(0) = neg_3D(0);
        neg_3D_input(1) = neg_3D(1);
        neg_3D_input(2) = neg_3D(2);


        camera->vsEuclideanToKeypoint(pos_3D_input,image_pos);
        camera->vsEuclideanToKeypoint(neg_3D_input,image_neg);
        
        Vec2f sample_pos = scaling_factor*Vec2f(image_pos(1),image_pos(0));   
        Vec2f sample_neg = scaling_factor*Vec2f(image_neg(1),image_neg(0)); 

        intensity_pos = image.InterpolateBilinear(sample_pos);
        intensity_neg = image.InterpolateBilinear(sample_neg);

        C_symmetry += (intensity_pos-intensity_neg)*(intensity_pos-intensity_neg);
        // cv::circle(cv_img, cv::Point2f(sample_pos(1),sample_pos(0)),1, cv::Scalar(255,255,255),1);
        // cv::circle(cv_img, cv::Point2f(sample_neg(1),sample_neg(0)),1, cv::Scalar(255,255,255),1);

        // cv::imshow("test",cv_img);
        // cv::waitKey(20000);

      }
      
      return C_symmetry;

    }


template <typename T, typename C>
double TargetToJacobian(
    const Image<T>& image,
    Eigen::Matrix4d T_target_to_euc,
    C* camera,
    Eigen::Vector4d target_location,
    vector<vector<Eigen::VectorXd>> samples_targetframe,
    int num_samples_symmetry,
    cv::Mat cv_img,
    Eigen::Vector2d* Jacobian,
    Matrix<double, 2, 2>* H,
    Matrix<double, 2, 1>* b,
    double scaling_factor
    )
    {
      
      H->triangularView<Eigen::Upper>().setZero();
      b->setZero();
      double C_symmetry = 0;
      double residual = 0;
      double intensity_pos, intensity_neg;
      Eigen::Matrix<double, 1, 2> gradient_pos, gradient_neg, gradient_pos_, gradient_neg_;
      Eigen::VectorXd distorted_pixel_location;
      Eigen::MatrixXd Jacobian_cam_pos, Jacobian_cam_neg;
      Eigen::Matrix<double,2,3> Jacobian_cam_pos_T, Jacobian_cam_neg_T;
      Eigen::Matrix<double,3,2> Jacobian_transformation; 
      Eigen::Vector2d Jacobian_raw= Eigen::Vector2d::Zero();
      double learning_rate = 1.0;
      for (int i = 0; i < num_samples_symmetry ; i++)
      {
        vector<Eigen::VectorXd> sample_pair = samples_targetframe[i];
        Eigen::VectorXd target_location_pos, target_location_neg, pos_3D, neg_3D, image_pos, image_neg;

        // target_location_pos = sample_pair[0] + target_location;
        // target_location_neg = sample_pair[1] + target_location;
        target_location_pos = target_location;
        target_location_pos(0) += sample_pair[0](0);
        target_location_pos(1) += sample_pair[0](1);
        
        target_location_neg = target_location;
        target_location_neg(0) += sample_pair[1](0);
        target_location_neg(1) += sample_pair[1](1);

        pos_3D =  T_target_to_euc*target_location_pos;
        neg_3D =  T_target_to_euc*target_location_neg;

        Eigen::Vector3d pos_3D_input, neg_3D_input;

        pos_3D_input(0) = pos_3D(0);
        pos_3D_input(1) = pos_3D(1);
        pos_3D_input(2) = pos_3D(2);

        neg_3D_input(0) = neg_3D(0);
        neg_3D_input(1) = neg_3D(1);
        neg_3D_input(2) = neg_3D(2);
        // jacobian of camera model, image coordinates w.r.t. euclidean coordinates
        camera->vsEuclideanToKeypoint(pos_3D_input,image_pos,Jacobian_cam_pos);
        camera->vsEuclideanToKeypoint(neg_3D_input,image_neg,Jacobian_cam_neg);
        // getting the transpose
        Jacobian_cam_pos_T = Jacobian_cam_pos;
        Jacobian_cam_neg_T = Jacobian_cam_neg;
        Vec2f sample_pos = scaling_factor*Vec2f(image_pos(1),image_pos(0));   
        Vec2f sample_neg = scaling_factor*Vec2f(image_neg(1),image_neg(0)); 
        // image gradients
        image.InterpolateBilinearWithJacobian(sample_pos, &intensity_pos, &gradient_pos);
        image.InterpolateBilinearWithJacobian(sample_neg, &intensity_neg, &gradient_neg);
        gradient_pos_(0,0) = gradient_pos(0,1);
        gradient_pos_(0,1) = gradient_pos(0,0);
        gradient_neg_(0,0) = gradient_neg(0,1);
        gradient_neg_(0,1) = gradient_neg(0,0);

        // Jacobian of transformation between euclidean and target coordinates
        Jacobian_transformation = T_target_to_euc.block<3,2>(0,0);
        Eigen::Vector2d jacobian_pos = gradient_pos_*Jacobian_cam_pos_T*Jacobian_transformation;
        Eigen::Vector2d jacobian_neg = gradient_neg_*Jacobian_cam_neg_T*Jacobian_transformation;
        // SM_INFO_STREAM("Suggested changes: "<<(jacobian_pos-jacobian_neg));
        Jacobian_raw.setZero();
        Jacobian_raw = (jacobian_pos-jacobian_neg);
        *Jacobian += (intensity_pos-intensity_neg)*(jacobian_pos-jacobian_neg);
        C_symmetry += (intensity_pos-intensity_neg)*(intensity_pos-intensity_neg);
        residual = (intensity_pos-intensity_neg);
        // cv::circle(cv_img, cv::Point2f(sample_pos(1),sample_pos(0)),3, cv::Scalar(255,0,0),2);
        // cv::circle(cv_img, cv::Point2f(sample_neg(1),sample_neg(0)),3, cv::Scalar(255,0,0),2);
        // SM_INFO_STREAM("Positive location: "<<sample_pos);
        // SM_INFO_STREAM("Negative location: "<<sample_neg);

        // cv::imshow("test",cv_img);
        // cv::waitKey(20000);
        Eigen::MatrixXd Jacobian_mat = Jacobian_raw;
        // SM_INFO_STREAM("Jacobian RAW"<<Jacobian_mat);
        H->triangularView<Eigen::Upper>() += learning_rate*Jacobian_raw * Jacobian_raw.transpose();
        *b += residual * Jacobian_raw*learning_rate;
      }
      
      return C_symmetry;

    }
}