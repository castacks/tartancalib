#ifndef ASLAM_TARTAN_CALIB
#define ASLAM_TARTAN_CALIB
#include <Eigen/Dense>
#include <opencv2/imgproc/imgproc.hpp>
#include <boost/serialization/export.hpp>
#include <aslam/cameras/GridCalibrationTargetObservation.hpp>
#include <aslam/cameras/GridCalibrationTargetAprilgrid.hpp>
#include <vector>
#include <utility>
#include <Eigen/Core>
#include <opencv2/core/core.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <iostream>
#include <cmath>
#include <boost/math/special_functions/round.hpp>
#include <boost/serialization/nvp.hpp>
#include <sm/logging.hpp>
#include <boost/math/constants/constants.hpp>
#include <opencv2/viz.hpp>
#include <opencv2/highgui.hpp>
#include <aslam/cameras/GridDetector.hpp>
#include <opencv2/core/eigen.hpp>
#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <sm/python/Id.hpp>
#include <aslam/matplotlibcpp.h>
#include <ctime>
#include <fstream>
#include <sm/kinematics/Transformation.hpp>
#include "apriltags/TagDetector.h"
#include "apriltags/Tag36h11.h"
#include "implementation/Datasetwriterhelper.hpp"
#include "implementation/symmetry_refinement.h"
#include "implementation/tartan_refinement.h"
#include "implementation/symmetry_fit.h"
// #include "deltille/PolynomialSaddleDetectorContext.h"
#include <aslam/cameras/CameraGeometryBase.hpp>
#include <random>
// #include "deltille/utils.h"

namespace aslam
{
    namespace cameras
    {
      enum class DebugMode 
      {
        pointcloudprojection,
        pinholeprojection,
        originalprojection,
        targetpointcloud,
        individualprojections,
        reprojectionpoints,
        none
      };
      
      enum class ReprojectionMode 
      {
        pinhole,
        homography,
        cornerpredictor,
        none
      };


    template< typename C>
    class ReprojectionWrapper:
    public boost::enable_shared_from_this<ReprojectionWrapper<C>>
    {
      public:
        ReprojectionWrapper(
          std::vector<aslam::cameras::GridCalibrationTargetObservation>& obslist,
          const Eigen::MatrixXd & fov,
          const Eigen::MatrixXd & pose,
          const Eigen::MatrixXd & resolution,
          const ReprojectionMode reproj_type):
        obslist_(obslist),
        fov_(fov),
        pose_(pose),
        resolution_(resolution),
        reproj_type_(reproj_type)
        {

        };
      
        std::vector<aslam::cameras::GridCalibrationTargetObservation> obslist_;
        Eigen::MatrixXd fov_,pose_,resolution_;
        cv::Mat map_x_, map_y_;
        std::vector<cv::Mat> homography_;
        ReprojectionMode reproj_type_;
        Eigen::Matrix<float, 3, Eigen::Dynamic> xyz_;
    };


    template< typename C>
    class TartanCalibWorker:
      public boost::enable_shared_from_this<TartanCalibWorker<C>>
    {
        public:  
            TartanCalibWorker(
              const C* camera,
              aslam::cameras::GridDetector gd,
              const std::vector<aslam::cameras::GridCalibrationTargetObservation>& obslist,
              const Eigen::MatrixXd & fovs,
              const Eigen::MatrixXd & poses,
              const Eigen::MatrixXd & resolutions,
              const std::vector<std::string>& reproj_types,
              const std::vector<std::string>& debug_modes,
              const int min_init_corners_autocomplete,  
              const int min_tag_size_autocomplete,
              const float correction_threshold,
              const int min_resize_window_size,
              const int max_resize_window_size,
              const float refine_magnitude_reject,
              const bool symmetry_refinement,
              const float symmetry_edge_threshold,
              const std::string export_dataset_dir,
              const std::string debug_image_dir
            )
            : obslist_(obslist),
              gd_(gd),
              target_(gd.target()),
              target_april_(boost::dynamic_pointer_cast<GridCalibrationTargetAprilgrid>(target_)),
              tagDetector_(boost::make_shared<AprilTags::TagDetector>(AprilTags::tagCodes36h11, target_april_->options().blackTagBorder)),
              camera_(camera),
              fovs_(fovs),
              poses_(poses),
              resolutions_(resolutions),
              in_width_(obslist[0].imCols()),
              in_height_(obslist[0].imRows()), // assumption: this array only contains images of the same resolution, as we expect to create one TartanCalibWorker class per camera
              minInitCornersAutoComplete(min_init_corners_autocomplete),
              minTagSizeAutoComplete(min_tag_size_autocomplete),
              correction_threshold(correction_threshold),
              minResizewindowSize(min_resize_window_size),
              maxResizewindowSize(max_resize_window_size),
              refine_magnitude_reject(refine_magnitude_reject),
              symmetry_refinement(symmetry_refinement),
              symmetry_edge_threshold(symmetry_edge_threshold),
              export_dataset_dir(export_dataset_dir),
              debug_image_dir(debug_image_dir)
             { 
              camera_->getParameters(cam_params,true,true,true);
              num_frames_ = obslist.size();
              num_views_ = fovs.cols();
              new_obslist_ = obslist_; // eventually we're outputting this variable, but we initialize it with the right observations (i.e., images and time stamps)
              for (auto obs : new_obslist_)
              {
                obs.getCornersIdx(outCornerIdx_);
              }
                        
              // Export dataset to a binary usable in generic camera model referred to in "Why having 10000 parameters..." paper.
              // This export corresponds to BEFORE TartanCalib enhancements.
              // if (!export_dataset_fp.empty())
              // {
              //   export_dataset(export_dataset_fp);
              // }

              SM_ASSERT_TRUE(std::runtime_error,num_views_ == poses_.cols() && poses_.cols() == resolutions_.cols() && resolutions_.cols() == reproj_types.size(), "All inserted tartan matrices need the same number of columns." );
              
              // loading reprojectionwrapper classes
              for (int i = 0; i < num_views_; i++)
              {
                reprojection_wrappers_.push_back(ReprojectionWrapper<C>(obslist_,fovs.col(i),poses.col(i),resolutions.col(i),StringToReprojectionMode(reproj_types[i])));
              }

              for (auto debug_mode : debug_modes)
              {
                debug_modes_.push_back(StringToDebug(debug_mode));
              }

              for (auto obs: obslist_)
              {
                obs.getCornersIdx(outCornerIdx_);
                num_corners_start+=outCornerIdx_.size();
              }
              SM_INFO_STREAM("Started tartan calib with "<<num_corners_start<<" points.");
            
            };
            TartanCalibWorker(){};
            ~TartanCalibWorker(){}
            void compute_xyzs();
            void compute_xyz(aslam::cameras::ReprojectionWrapper<C>& reprojection);
            void compute_rotation(const Eigen::MatrixXd& pose );
            void compute_remap(aslam::cameras::ReprojectionWrapper<C>& reprojection);
            void project_to_original_image(void);
            void compute_remaps(); 
            void compute_reprojections();
            void compute_corners();
            cv::Mat get_mat(std::vector<aslam::cameras::GridCalibrationTargetObservation>,bool,float,std::vector<cv::Scalar> ,float);   
            void debug_show(void);
            std::vector<aslam::cameras::GridCalibrationTargetObservation> getNewObs(void);
            void homography_reprojection(aslam::cameras::ReprojectionWrapper<C>& reprojection );
            void match_quads(aslam::cameras::ReprojectionWrapper<C>& reprojection);
            bool export_dataset(std::string path);
            void write_observation(aslam::cameras::GridCalibrationTargetObservation obs, std::string filename);
           
            std::string ReprojectionModeToString(ReprojectionMode e)
            {
              switch (e)
              {
                case ReprojectionMode::pinhole: return "pinhole";
                case ReprojectionMode::homography: return "homography";
                case ReprojectionMode::cornerpredictor: return "cornerpredictor";
                case ReprojectionMode::none: return "none";
              }
            }
            ReprojectionMode StringToReprojectionMode(std::string e)
            {
              std::map<std::string, ReprojectionMode> reprojection_mode = boost::assign::map_list_of
              ("pinhole",ReprojectionMode::pinhole)
              ("homography",ReprojectionMode::homography)
              ("cornerpredictor",ReprojectionMode::cornerpredictor)
              ("none",ReprojectionMode::none);
              return reprojection_mode[e];
            }
            
            // DebugMode converters
            std::string DebugToString(DebugMode e)
            {
              switch (e)
              {
                case DebugMode::pointcloudprojection: return "pointcloudprojection";
                case DebugMode::pinholeprojection: return "pinholeprojection";
                case DebugMode::originalprojection: return "originalprojection";
                case DebugMode::targetpointcloud: return "targetpointcloud";
                case DebugMode::individualprojections: return "individualprojections";
                case DebugMode::reprojectionpoints: return "reprojectionpoints";
                case DebugMode::none: return "none";
              }
            }
            DebugMode StringToDebug(std::string e)
            {
              std::map<std::string, DebugMode> debug_mode = boost::assign::map_list_of
              ("pointcloudprojection",DebugMode::pointcloudprojection)
              ("pinholeprojection", DebugMode::pinholeprojection)
              ("originalprojection",DebugMode::originalprojection)
              ("targetpointcloud",DebugMode::targetpointcloud)
              ("individualprojections",DebugMode::individualprojections)
              ("reprojectionpoints",DebugMode::reprojectionpoints)
              ("none",DebugMode::none);
              return debug_mode[e];
            }
            // Converts a C++ vector to a python list
            template <class T>
            boost::python::list toPythonList(std::vector<T> vector) {
              typename std::vector<T>::iterator iter;
              boost::python::list list;
              for (iter = vector.begin(); iter != vector.end(); ++iter) {
                list.append(*iter);
              }
            return list;
            }

                
            /// \brief Serialization
            enum {CLASS_SERIALIZATION_VERSION = 4};
            BOOST_SERIALIZATION_SPLIT_MEMBER()

            /// \brief Serialization support
            template<class Archive>
            void save(Archive & ar, const unsigned int /*version*/) const
            {

            }
            template<class Archive>
            void load(Archive & ar, const unsigned int /*version*/)
            {
            }
          
        private:
            // key class variables
            std::vector<DebugMode> debug_modes_; //determines what the user will see
            std::vector<aslam::cameras::ReprojectionWrapper<C>> reprojection_wrappers_; 
            std::vector<std::string> reproj_types_;

            Eigen::Matrix<float, 3, Eigen::Dynamic> xyz_; // rows are xyz, each column is a point
            Eigen::Matrix<float,3,3> rot_mat_,rot_mat_x_,rot_mat_z_;
            Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> xx_,yy_; //vector with x and y values in the global
            Eigen::MatrixXd fovs_, poses_, resolutions_;
            std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> xyzs_;
            
            std::vector<unsigned int> outCornerIdx_;
            aslam::cameras::GridCalibrationTargetObservation obs_;
            std::vector<aslam::cameras::GridCalibrationTargetObservation> obslist_;
            int num_frames_,num_views_,num_points_,in_width_,in_height_,xyz_idx_,num_corners_;
            int num_corners_start = 0;
            int num_corners_end = 0;
            float rot_x_, rot_z_;
            bool verbose_,empty_pixels_;
            const C* camera_;

            // used to create remap
            Eigen::VectorXf xyz_ray_;
            Eigen::VectorXd distorted_pixel_location_; 

            // map parameters
            cv::Size resolution_out_;
            cv::Mat map_x_float_,map_y_float_;
            std::vector<cv::Mat> maps_x_, maps_y_;
            
            std::vector<std::vector<cv::Mat>> reprojections_;
            std::vector<aslam::cameras::GridCalibrationTargetObservation> new_obslist_;
            aslam::cameras::GridDetector gd_;
            std::vector<aslam::cameras::GridCalibrationTargetObservation> output_obslist_;
            std::string log_file = "log.txt";
            GridCalibrationTargetBase::Ptr target_;
            GridCalibrationTargetAprilgrid::Ptr target_april_;
            boost::shared_ptr<AprilTags::TagDetector> tagDetector_;
            Eigen::MatrixXd cam_params;

            int minInitCornersAutoComplete = 24; // we need at least this many corners to be able to do autocomplete, since the pose of the board is otherwise too uncertain.
            float minTagSizeAutoComplete = 0; // this is how many pixels a tag needs to be in size before we consider autocompleting it. This is just to make sure really small tags aren't detected and then detected poorly
            float correction_threshold = 2.0; // number of pixel offset between reprojection and detection we allow
            float minResizewindowSize = 2;
            float maxResizewindowSize = 8;
            double refine_magnitude_reject = 10.0;
            bool symmetry_refinement = false;
            double symmetry_edge_threshold = 10.0;
            std::string export_dataset_dir = "";
            std::string debug_image_dir="";
            bool dummy = false;
    };
    }
}


// namespace boost {
// namespace serialization {

// template<class Archive>
// void serialize(Archive & ar, const unsigned int /* version */) {
  
// }

// }  // namespace serialization
// }  // namespace boost

SM_BOOST_CLASS_VERSION_T1(aslam::cameras::TartanCalibWorker);
#include<aslam/implementation/TartanCalibWorker.hpp>
#endif /* ASLAM_TARTAN_CALIB */
