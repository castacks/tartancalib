#include <vector>
#include <numpy_eigen/boost_python_headers.hpp>
#include <aslam/cameras/CameraGeometryBase.hpp>
#include <sm/python/Id.hpp>
#include <boost/python/stl_iterator.hpp>
#include <aslam/cameras/GridCalibrationTargetObservation.hpp>
#include <aslam/TartanCalibWorker.hpp>
#include <sm/python/boost_serialization_pickle.hpp>

namespace detail {
// "Pass by reference" doesn't work with the Eigen type converters.
// So these functions must be wrapped. While we are at it, why 
// not make them nice.
template<typename C>
Eigen::VectorXd e2k(const C * camera, Eigen::Vector3d const & p) {
  Eigen::VectorXd k;
  camera->vsEuclideanToKeypoint(p, k);

  return k;
}

template<typename C>
boost::python::tuple e2kJp(const C * camera, Eigen::Vector3d const & p) {
  Eigen::MatrixXd Jp;
  Eigen::VectorXd k;
  bool isValid = camera->vsEuclideanToKeypoint(p, k, Jp);
  return boost::python::make_tuple(k, Jp, isValid);
}

template<typename C>
Eigen::VectorXd eh2k(const C * camera, Eigen::Vector4d const & p) {
  Eigen::VectorXd k;
  camera->vsHomogeneousToKeypoint(p, k);
  return k;
}

template<typename C>
boost::python::tuple eh2kJp(const C * camera, Eigen::Vector4d const & p) {
  Eigen::MatrixXd Jp;
  Eigen::VectorXd k;
  bool valid = camera->vsHomogeneousToKeypoint(p, k, Jp);
  return boost::python::make_tuple(k, Jp, valid);
}

template<typename C>
Eigen::Vector3d k2e(const C * camera, Eigen::VectorXd const & k) {
  Eigen::Vector3d p;
  camera->vsKeypointToEuclidean(k, p);
  return p;
}

template<typename C>
boost::python::tuple k2eJk(const C * camera, Eigen::VectorXd const & k) {
  Eigen::MatrixXd Jk;
  Eigen::VectorXd p;
  bool valid = camera->vsKeypointToEuclidean(k, p, Jk);
  return boost::python::make_tuple(p, Jk, valid);
}

template<typename C>
Eigen::Vector4d k2eh(const C * camera, Eigen::VectorXd const & k) {
  Eigen::VectorXd ph;
  camera->vsKeypointToHomogeneous(k, ph);
  return ph;
}

template<typename C>
boost::python::tuple k2ehJk(const C * camera, Eigen::VectorXd const & k) {
  Eigen::MatrixXd Jk;
  Eigen::VectorXd p;
  bool valid = camera->vsKeypointToHomogeneous(k, p, Jk);
  return boost::python::make_tuple(p, Jk, valid);
}

template<typename C>
boost::python::tuple estimateTransformation(const C * camera, aslam::cameras::GridCalibrationTargetObservation & obs)
{
  sm::kinematics::Transformation trafo;
  bool success = camera->estimateTransformation(obs, trafo);
  return boost::python::make_tuple(success, trafo);
}

// TODO: Condense all function arguments into a struct for neatness.
template<typename C>
boost::python::list getProjections(const C * camera,
  aslam::cameras::GridDetector gd,
  const boost::python::object& py_obslist,
  const Eigen::MatrixXd & fovs,
  const Eigen::MatrixXd & poses,
  const Eigen::MatrixXd & resolutions,
  const boost::python::list & reproj_pylist,
  const boost::python::list & debug_modes,
  const boost::python::list & TartanParameters)
{
  //convert obs python list to stl vector
  boost::python::stl_input_iterator<aslam::cameras::GridCalibrationTargetObservation> begin(py_obslist), end;
  std::vector<aslam::cameras::GridCalibrationTargetObservation> obslist(begin, end);

  //convert string python list to stl 
  std::vector<std::string> reproj_types;
  for (int i=0; i< len(reproj_pylist);i++)
  {
    reproj_types.push_back(boost::python::extract<std::string>(reproj_pylist[i]));
  }

  std::vector<std::string> debug_modes_;
  for (int i=0; i< len(debug_modes);i++)
  {
    debug_modes_.push_back(boost::python::extract<std::string>(debug_modes[i]));
  }
  
  const int min_init_corners_autocomplete = boost::python::extract<int>(TartanParameters[0]);
  const int min_tag_size_autocomplete = boost::python::extract<int>(TartanParameters[1]);
  const float correction_threshold = boost::python::extract<float>(TartanParameters[2]);
  const int min_resize_window_size = boost::python::extract<int>(TartanParameters[3]);
  const int max_resize_window_size = boost::python::extract<int>(TartanParameters[4]);
  const float refine_magnitude_reject = boost::python::extract<float>(TartanParameters[5]);
  const bool symmetry_refinement = boost::python::extract<bool>(TartanParameters[6]);
  const float symmetry_edge_threshold = boost::python::extract<float>(TartanParameters[7]);
  const std::string export_dataset_dir = boost::python::extract<std::string>(TartanParameters[8]);
  const std::string debug_image_dir = boost::python::extract<std::string>(TartanParameters[9]);

  auto tartan_ = aslam::cameras::TartanCalibWorker<C>(
    camera,
    gd,
    obslist,
    fovs,
    poses,
    resolutions,
    reproj_types,
    debug_modes_,
    min_init_corners_autocomplete,
    min_tag_size_autocomplete,
    correction_threshold,
    min_resize_window_size,
    max_resize_window_size,
    refine_magnitude_reject,
    symmetry_refinement,
    symmetry_edge_threshold,
    export_dataset_dir,
    debug_image_dir
  );

  tartan_.compute_xyzs();
  tartan_.compute_remaps();
  tartan_.compute_reprojections();
  tartan_.compute_corners();
  tartan_.project_to_original_image();
  tartan_.debug_show();
  std::vector<aslam::cameras::GridCalibrationTargetObservation> obslist_new = tartan_.getNewObs();
  
  boost::python::list ret = tartan_.toPythonList(obslist_new);

  return ret;
}

template<typename C>
bool initializeIntrinsics(C* camera, const boost::python::object& py_obslist)
{
  //convert python list to stl vector
  boost::python::stl_input_iterator<aslam::cameras::GridCalibrationTargetObservation> begin(py_obslist), end;
  std::vector<aslam::cameras::GridCalibrationTargetObservation> obslist(begin, end);

  bool success = camera->initializeIntrinsics(obslist);
  return success;
}

template <class T>
boost::python::list toPythonList(std::vector<T> vector) {
  typename std::vector<T>::iterator iter;
  boost::python::list list;
  for (iter = vector.begin(); iter != vector.end(); ++iter) {
    list.append(*iter);
  }
return list;
}

template<typename C>
boost::python::list getReprojectionErrors(C* camera, const boost::python::object& py_obslist)
{
  //convert python list to stl vector
  boost::python::stl_input_iterator<aslam::cameras::GridCalibrationTargetObservation> begin(py_obslist), end;
  std::vector<aslam::cameras::GridCalibrationTargetObservation> obslist(begin, end);
  sm::kinematics::Transformation T_target_camera;
  boost::python::list all_errs;
  Eigen::Vector2d y;
  Eigen::VectorXd yhat;
  std::vector<Eigen::VectorXd>  obs_err_list;
  for (auto obs: obslist)
  {
    camera->estimateTransformation(obs,T_target_camera);
    sm::kinematics::Transformation T_camera_target = T_target_camera.inverse();
    
    for (size_t i = 0; i < obs.target()->size(); ++i) {
      
      const Eigen::Vector3d& euclidean =  dynamic_cast< const Eigen::Vector3d &>(T_camera_target * obs.target()->point(i))  ;

      if (obs.imagePoint(i, y)
          && camera->vsEuclideanToKeypoint(euclidean,yhat )) {
        obs_err_list.push_back((y - yhat));
      }
    }

  }
  boost::python::list all_errs_python = toPythonList(obs_err_list);

  return all_errs_python; 
}

}  // namespace detail

template<typename T>
Eigen::MatrixXd getParameters(T * D, bool p, bool d, bool s) {

  Eigen::MatrixXd P;
  D->getParameters(P, p, d, s);
  return P;
}

void exportCameraGeometryBase() {
  sm::python::Id_python_converter<aslam::cameras::CameraId>::register_converter();

  using aslam::cameras::CameraGeometryBase;
  boost::python::class_<aslam::cameras::CameraGeometryBase,
      boost::shared_ptr<aslam::cameras::CameraGeometryBase>, boost::noncopyable>(
      "CameraGeometryBase", boost::python::no_init)
      .def("keypointDimension", &CameraGeometryBase::keypointDimension, "Get the dimension of the keypoint type")
      .def("isProjectionInvertible", &CameraGeometryBase::isProjectionInvertible, "Is the sensor model invertible? Ususally this is only true for a range/bearing sensor.")
      .def("temporalOffset", &CameraGeometryBase::vsTemporalOffset, "Given a keypoint, what is the offset from the start of integration for this image?\nDuration = temporalOffset(keypoint)")
      .def("createRandomKeypoint", &CameraGeometryBase::createRandomKeypoint, "Create a valid, random keypoint. This is useful for unit testing and experiments.")
      .def("createRandomVisiblePoint",&CameraGeometryBase::createRandomVisiblePoint, "Create a valid point in space visible by the camera.\np = createRandomVisiblePoint(depth).")
      .def("getProjections", &detail::getProjections<CameraGeometryBase>, "Get the reprojections used for tartancalib.")
      .def("getReprojectionErrors",&detail::getReprojectionErrors<CameraGeometryBase>,"Get all reprojection errors given an observation.")
      .def("euclideanToKeypoint", &detail::e2k<CameraGeometryBase>, "Map a 3x1 Euclidean point to a keypoint.\nk = euclideanToKeypoint(p)")
      .def("euclideanToKeypointJp", &detail::e2kJp<CameraGeometryBase>, "Map a 3x1 Euclidean point to a keypoint and get the Jacobian of the mapping with respect to small changes in the point.\n(k, Jp) = euclideanToKeypoint(p)")
      .def("homogeneousToKeypoint", &detail::eh2k<CameraGeometryBase>, "Map a 4x1 homogeneous Euclidean point to a keypoint.\nk = euclideanToKeypoint(p)")
      .def("homogeneousToKeypointJp", &detail::eh2kJp<CameraGeometryBase>, "Map a 4x1 homogeneous Euclidean point to a keypoint and get the Jacobian of the mapping with respect to small changes in the point.\n(k, Jp) = homogeneousToKeypointJp(p)")
      .def("keypointToHomogeneous", &detail::k2eh<CameraGeometryBase>,"Map a keypoint to a 4x1 homogeneous Euclidean point.\np = keypointToHomogeneous(k)")
      .def("keypointToHomogeneousJk", &detail::k2ehJk<CameraGeometryBase>, "Map a keypoint to a 4x1 homogeneous Euclidean point and get the Jacobian of the mapping with respect to small changes in the keypoint.\n(p, Jk) = keypointToHomogeneousJk(k)")
      .def("keypointToEuclidean", &detail::k2e<CameraGeometryBase>, "Map a keypoint to a 3x1 Euclidean point.\np = keypointToEuclidean(k)")
      .def("keypointToEuclideanJk", &detail::k2eJk<CameraGeometryBase>, "Map a keypoint to a 3x1 Euclidean point and get the Jacobian of the mapping with respect to small changes in the keypoint.\n(p, Jk) = keypointToEuclideanJk(k)")
      .def("isValid", &CameraGeometryBase::vsIsValid)
      .def("isEuclideanVisible", &CameraGeometryBase::vsIsEuclideanVisible)
      .def("isHomogeneousVisible", &CameraGeometryBase::vsIsHomogeneousVisible)
      .def("minimalDimensionsProjection", &CameraGeometryBase::minimalDimensionsProjection)
      .def("minimalDimensionsDistortion", &CameraGeometryBase::minimalDimensionsDistortion)
      .def("minimalDimensionsShutter", &CameraGeometryBase::minimalDimensionsShutter)
      .def("getParameters", &getParameters<CameraGeometryBase>)
      .def("setParameters", &CameraGeometryBase::setParameters)
      .def("estimateTransformation", &detail::estimateTransformation<CameraGeometryBase>, "estimate the transformation of the camera with respect to the calibration target, returns tuple (bool, sm.Transformation)")
      .def("initializeIntrinsics", &detail::initializeIntrinsics<CameraGeometryBase>, "intialize intrinsics on a list of observations")
      .def_pickle( sm::python::pickle_suite<CameraGeometryBase>())
      ;
}