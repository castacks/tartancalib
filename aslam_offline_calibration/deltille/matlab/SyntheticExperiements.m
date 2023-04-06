% Copyright (C) 2017-present, Facebook, Inc.
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

classdef SyntheticExperiements < handle

  properties (SetAccess = private)
    synth_gen_; % synthetic corner generator
    theta_;     % range of perspective angles
    phi_;       %
    blur_;
    noise_;
  end

  properties (Constant)
    kFixedNoise_ = 0.01;
    kFixedBlur_  = 1.0;
    kDefaultNumTrails_ = 1000;
    kNumRefinementIterations = 5;
  end


  methods
    function this = SyntheticExperiements(rows, cols, theta, phi, blur, noise)
      if nargin < 6, noise = eps:0.005:0.05+eps; end
      if nargin < 5, blur = 0.0:0.5:5.0; end
      if nargin < 4, phi = 0:359; end
      if nargin < 3, theta = eps:5:75+eps; end
      if nargin < 2, cols = 100; end
      if nargin < 1, rows = 100; end

      this.synth_gen_ = SyntheticCornerGenerator(rows, cols);
      this.theta_ = deg2rad(theta);
      this.phi_ = deg2rad(phi);
      this.blur_ = blur;
      this.noise_ = noise;
    end

    function [err_ours, err_harris, noise_pct] = runNoise(this, kernel_size, num_trials)
      % function [err_ours, err_harris] = runNoise(this)
      %
      % Run variable noise expierments
      %
      % INPUT:
      %   num_trials number of times to repeat the experiements with random poses
      %
      % OUTPUT
      %   err_ours    mean error from our algorithm. Cell array, the first are
      %               the results on the checkerboard target. The second are the
      %               results on the triangular 
      %   err_harris  ditto with harris
      %   noise_pct   percentage of noise (the x-axis in the plot)
      %

      if nargin < 3, num_trials = SyntheticExperiements.kDefaultNumTrails_; end
      if nargin < 2, kernel_size = 5; end

      x_init = this.synth_gen_.getInitialPoints();

      for i = 2 : -1 : 1
        err_ours{i}   = zeros(2, size(x_init, 2), length(this.noise_));
        err_harris{i} = zeros(size(err_ours{i}));
      end

      pattern_type = {'square', 'triangle'};

      for k = 1 : num_trials
        t = deg2rad(randrange(this.theta_(1), 60.0));
        p = deg2rad(randrange(this.phi_(1), this.phi_(end)));
        b = SyntheticExperiements.kFixedBlur_;
        for n = 1 : length(this.noise_)
          sigma_n = this.noise_(n);
          for i = 1 : 2
            [I, x_true] = this.synth_gen_.makeImage(i, t, p, b, sigma_n);

            x = computeHarrisCorner_sim(x_init, I, 5);
            err_harris{i}(:,:,k,n) = x - x_true;

            % NOTE the +-1 for conversion between Matlab and C
            x = 1.0 + saddle_subpix_refinement_mex(pattern_type{i}, I, x_init - 1.0, ...
              kernel_size, SyntheticExperiements.kNumRefinementIterations);
            err_ours{i}(:,:,k,n) = x - x_true;
          end
        end
      end

      noise_pct = 100.0 * this.noise_;
      for i = 1 : 2
        err_ours{i} = SyntheticExperiements.compute_mean_err(err_ours{i});
        err_harris{i} = SyntheticExperiements.compute_mean_err(err_harris{i});
      end
    end

    function [err_ours, err_harris, blur] = runBlur(this, kernel_size, num_trials)
      
      if nargin < 3, num_trials = SyntheticExperiements.kDefaultNumTrails_; end
      if nargin < 2, kernel_size = 5; end

      x_init = this.synth_gen_.getInitialPoints();

      for i = 2 : -1 : 1
        err_ours{i}   = zeros(2, size(x_init, 2), length(this.noise_));
        err_harris{i} = zeros(size(err_ours{i}));
      end

      pattern_type = {'square', 'triangle'};

      for k = 1 : num_trials
        t = deg2rad(randrange(this.theta_(1), 60.0));
        p = deg2rad(randrange(this.phi_(1), this.phi_(end)));
        sigma_n = SyntheticExperiements.kFixedNoise_;
        for n = 1 : length(this.blur_)
          b = this.blur_(n);
          for i = 1 : 2
            [I, x_true] = this.synth_gen_.makeImage(i, t, p, b, sigma_n);

            x = computeHarrisCorner_sim(x_init, I, 5);
            err_harris{i}(:,:,k,n) = x - x_true;

            % NOTE the +-1 for conversion between Matlab and C
            x = 1.0 + saddle_subpix_refinement_mex(pattern_type{i}, I, x_init - 1.0, ...
              kernel_size, SyntheticExperiements.kNumRefinementIterations);
            err_ours{i}(:,:,k,n) = x - x_true;
          end
        end
      end

      blur = this.blur_;
      for i = 1 : 2
        err_ours{i} = SyntheticExperiements.compute_mean_err(err_ours{i});
        err_harris{i} = SyntheticExperiements.compute_mean_err(err_harris{i});
      end
    end


    function [err_ours, err_harris, angles] = runPerspective(this, kernel_size, num_trials)
      
      if nargin < 3, num_trials = SyntheticExperiements.kDefaultNumTrails_; end
      if nargin < 2, kernel_size = 5; end

      x_init = this.synth_gen_.getInitialPoints();

      for i = 2 : -1 : 1
        err_ours{i}   = zeros(2, size(x_init, 2), length(this.noise_));
        err_harris{i} = zeros(size(err_ours{i}));
      end

      pattern_type = {'square', 'triangle'};

      for k = 1 : num_trials
        p = deg2rad(randrange(this.phi_(1), this.phi_(end)));
        b = this.kFixedBlur_;
        sigma_n = SyntheticExperiements.kFixedNoise_;
        for n = 1 : length(this.theta_)
          t = this.theta_(n);
          for i = 1 : 2
            [I, x_true] = this.synth_gen_.makeImage(i, t, p, b, sigma_n);

            x = computeHarrisCorner_sim(x_init, I, 5);
            err_harris{i}(:,:,k,n) = x - x_true;

            % NOTE the +-1 for conversion between Matlab and C
            x = 1.0 + saddle_subpix_refinement_mex(pattern_type{i}, I, x_init - 1.0, ...
              kernel_size, SyntheticExperiements.kNumRefinementIterations);
            err_ours{i}(:,:,k,n) = x - x_true;
          end
        end
      end

      angles = rad2deg(this.theta_);
      for i = 1 : 2
        err_ours{i} = SyntheticExperiements.compute_mean_err(err_ours{i});
        err_harris{i} = SyntheticExperiements.compute_mean_err(err_harris{i});
      end
    end

  end

  methods (Static, Access=private)
    function m = compute_mean_err(e)
      m = squeeze(mean(mean(sqrt(e(1,:,:,:).^2 + e(2,:,:,:).^2), 2), 3));
    end
  end

end
