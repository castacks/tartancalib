
# This file is a modified version of Dockerfile_ros1_18_04.

ARG base_image=osrf/ros:melodic-desktop-full
FROM ${base_image}

# Allow using GUI apps.
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata \
 && ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime \
 && dpkg-reconfigure --frontend noninteractive tzdata \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Priliminary packages before installing others.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential \
        sudo \
        apt-utils \
        software-properties-common \
 && add-apt-repository ppa:deadsnakes/ppa -y \
 && apt-get update \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Refer to 
# https://github.com/ethz-asl/kalibr/wiki/installation
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
	git \
    wget \
    autoconf \
    automake \
    nano \
    vim \
    htop \
    tmux \
	python3-dev \
    python-pip \
    python-scipy \
    python-matplotlib \
	ipython \
    python-wxgtk4.0 \
    python-tk \
    python-igraph \
	libeigen3-dev \
    libboost-all-dev \
    libsuitesparse-dev \
	doxygen \
	libopencv-dev \
	libpoco-dev \
    libtbb-dev \
    libblas-dev \
    liblapack-dev \
    libv4l-dev \
    libglew-dev\
	python-catkin-tools \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    dill

# Create the workspace and build kalibr/tartancalib in it.
ENV WORKSPACE /catkin_ws

# Configure the build.
RUN mkdir -p $WORKSPACE/src \
 && cd $WORKSPACE \
 && catkin init \
 && catkin config --extend /opt/ros/melodic \
 && catkin config --cmake-args -DCMAKE_BUILD_TYPE=Release

# Copy the source code and build.
ADD . $WORKSPACE/src/kalibr
RUN	cd $WORKSPACE \
 &&	catkin build -j$(nproc) \
 && rm -rf build/

# When a user runs a command we will run this code before theirs
# This will allow for using the manual focal length if it fails to init
# https://github.com/ethz-asl/kalibr/pull/346
ENTRYPOINT export KALIBR_MANUAL_FOCAL_LENGTH_INIT=1 && \
	# /bin/bash -c "source \"$WORKSPACE/devel/setup.bash\"" && \ 
	cd $WORKSPACE && \
	/bin/bash
