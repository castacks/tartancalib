## ! DO NOT MANUALLY INVOKE THIS setup.py, USE CATKIN INSTEAD

from distutils.core import setup
from catkin_pkg.python_setup import generate_distutils_setup

# fetch values from package.xml
setup_args = generate_distutils_setup(
    packages=['kalibr_errorterms',
              'kalibr_common',
              'kalibr_camera_calibration',
              'kalibr_imu_camera_calibration',
              'tartan_logging'
              ],
    package_dir={'':'python'}
)

setup(**setup_args)
