#include <aslam/TartanCalibWorker.hpp>
#include <numpy_eigen/boost_python_headers.hpp>
#include <boost/python/operators.hpp>
#include <sm/python/boost_serialization_pickle.hpp>
#include <aslam/cameras/CameraGeometryBase.hpp>

using namespace boost::python; 
// using namespace aslam::cameras;
// class aslam::cameras::Tartan_Calib_Worker;
using aslam::cameras::CameraGeometryBase;

void exportTartan()
{       
        // using namespace aslam::cameras;
        class_<aslam::cameras::TartanCalibWorker<CameraGeometryBase>>("TartanCalibWorker")
        // .def(init<>())
        // .def(init<const std::vector<aslam::cameras::GridCalibrationTargetObservation>&, const Eigen::MatrixXd&, const Eigen::MatrixXd &, const Eigen::MatrixXd&, const bool>())
        .def("compute_xyzs",&aslam::cameras::TartanCalibWorker<CameraGeometryBase>::compute_xyzs)
        .def("compute_remaps",&aslam::cameras::TartanCalibWorker<CameraGeometryBase>::compute_remaps)
        .def("compute_reprojections",&aslam::cameras::TartanCalibWorker<CameraGeometryBase>::compute_reprojections)
        .def("compute_corners",&aslam::cameras::TartanCalibWorker<CameraGeometryBase>::compute_corners)
        .def("project_to_original_image",&aslam::cameras::TartanCalibWorker<CameraGeometryBase>::project_to_original_image)
        .def("debug_show",&aslam::cameras::TartanCalibWorker<CameraGeometryBase>::debug_show)
        .def("getNewObs",&aslam::cameras::TartanCalibWorker<CameraGeometryBase>::getNewObs)
        .def_pickle(sm::python::pickle_suite<aslam::cameras::TartanCalibWorker<CameraGeometryBase>>())
        ;
}
