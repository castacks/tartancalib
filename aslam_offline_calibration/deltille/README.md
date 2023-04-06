# Deltille detector
Deltille detector is a robust deltille and checkerboard detector. It comes with
detector library, example detector code, MATLAB bindings and code to replicate
the experiments in our
[paper](https://research.fb.com/publications/deltille-grids-for-geometric-camera-calibration).
Deltille detector uses parts of
[AprilTags](https://svn.csail.mit.edu/apriltags/) library and
[Kalibr](https://github.com/ethz-asl/kalibr) library.

## Dependencies

Deltille detector library is built of Linux or Windows but has no hard
dependencies on either. We tired to keep the external dependencies at minimum,
deltille detector library depends on **OpenCV 3.0+**, the example application
uses **Boost 1.63+**, optionally you may use **Matlab** to build Matlab
bindings and replicate some experiments from the paper.

## Building

The instructions provided here have been tested under Ubuntu 16.04.3 LTS x86_64
with Matlab version: R2016b (9.1.0.441665), but should work and compile on
Windows or MacOS and other versions of Matlab:

1. compile the main deltille detector library and binary
```
  $ DELTILLE_ROOT=<location of your clone of the deltille detector>
  $ cd $DELTILLE_ROOT && mkdir build && cd build
  $ cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_APPS=ON && make -j8
```

2. compile the main deltille detector library with Matlab support
```
  $ DELTILLE_ROOT=<location of your clone of the deltille detector>
  $ cd $DELTILLE_ROOT && mkdir build && cd build
  $ cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_MATLAB=ON && make -j8
```

## Rerun synthetic experiments
A set of Matlab scripts is provided to replicate the synthetic experiments in
the paper.

1. Build Matlab bindings

2. Edit your Matlab startup path to point to the compiled `*.mex*` files. When
   building on Linux with the default parameters the compiled binaries will be
   under the `matlab` directory. In Matlab, you can do the following:
   ` >> addpath('<DELTILLE_ROOT/build/matlab')`

3. Alternatively, you may start from the build directory, i.e.
```
  $ cd DELTILLE_ROOT/build
  $ matlab
```

4. To replicate the synthetic experiments, run 'runSyntheticExperiments'.

## Using the Matlab Wrapper
See `runSyntheticExperiments.m` for how to use the wrapper. The steps are:

1. Instatiate the `TargetDetector` class:
```
target_dectector = TargetDetector(dsc_file)
```
This will initialize the target loaded with the appropriate calibration target
definitions

2. Run the code on your images:
```
I = imread(...); % read th eimage
[c, ids] = target_detector.run(I);
hold off; imshow(I); hold on;
plot(1+c(1,:), 1+c(2,:), 'y+');
```
IMPORTANT: you must add 1 to the corners for visualization as Matlab is 1-based.

For additional details please consult our [paper](https://research.fb.com/publications/deltille-grids-for-geometric-camera-calibration)
 or contact the authors.

## Deltille detector dataset
The robot hand and fisheye images dataset used in the paper can be found in the related [repository](https://github.com/facebookincubator/deltille-dataset).

## License
Deltille detector code is licensed under the LGPL v2.1 license. For more
information, please see COPYING file.

## Citation
If you find this work useful, please cite the related paper:
```
@InProceedings{Ha_2017_ICCV,
author = {Ha, Hyowon and Perdoch, Michal and Alismail, Hatem and So Kweon, In and Sheikh, Yaser},
title = {Deltille Grids for Geometric Camera Calibration},
booktitle = {The IEEE International Conference on Computer Vision (ICCV)},
month = {Oct},
pages = {5344--5352},
year = {2017}}
```
