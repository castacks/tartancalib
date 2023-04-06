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

blur = 1;
noise = eps:0.005:0.05+eps;
% TODO (fix this)
theta = eps:5:60+eps;    % perspective angle
phi = 0:359;

num_iters = 5;
num_trials = 1000;

harris_size = 5;
saddle_size = 5;
monkey_size = 5;

[kernel{1}, invA{1}] = initSaddleDetection_sim(saddle_size);
[kernel{2}, invA{2}] = initMonkeyDetection_sim(monkey_size);

do_approx = false;
if do_approx
  for i = 1 : length(kernel)
    kernel{i} = computeApproxSeparable(kernel{i});
  end
  funs = {@computeSaddlePoint_sim_approx, @computeMonkeyPoint_sim_approx};
else
  funs = {@computeSaddlePoint_sim, @computeMonkeyPoint_sim};
end


synth_gen = SyntheticCornerGenerator(100, 100);

x_init = synth_gen.getInitialPoints();
np = size(x_init, 2);


err_harris{1} = zeros(2, np, num_trials, length(noise));
err_harris{2} = zeros(2, np, num_trials, length(noise));

err_ours{1} = zeros(2, np, num_iters + 1, num_trials, length(noise));
err_ours{2} = zeros(2, np, num_iters + 1, num_trials, length(noise));

err_cpp{1} = zeros(2, np, num_trials, length(noise));
err_cpp{2} = zeros(2, np, num_trials, length(noise));

for k = 1 : num_trials
  fprintf('trial %d/%d\r', k, num_trials);

  t=(rand(1)*(theta(end)-theta(1))+theta(1))/180*pi;
  p=(rand(1)*(phi(end)-phi(1))+phi(1))/180*pi;
  b=rand(1)*(blur(end)-blur(1))+blur(1);

  for n = 1 : length(noise)
    sigma_n = noise(n);
    for i = 1 : 2
      [I, x_true] = synth_gen.makeImage(i, t, p, b, sigma_n);

      % harris
      x = computeHarrisCorner_sim(x_init, I, harris_size);
      err_harris{i}(:,:,k,n) = x - x_true;

      x = funs{i}(I, x_init, kernel{i}, invA{i}, num_iters);
      err_ours{i}(:,:,:,k,n) = x - x_true;

      if i == 1
        pattern_type = 'square';
      else 
        pattern_type = 'triangle';
      end

      x = saddle_subpix_refinement_mex(pattern_type, I, x_init - 1, 5, 5);
      err_cpp{i}(:,:,k,n) = (x+1.0) - x_true;
    end
  end
end

fprintf('\n');

figure(2); clf;
plotSimErrors(100*noise, err_harris, err_ours);
legend('Harris [square]', 'Harris [triangle]', 'Ours [square]', 'Ours [triangle]');

figure(3); clf;
plotSimErrors(100*noise, err_cpp, err_ours);
legend('Ours CPP [square]', 'OURS CPP [triangle]', 'Ours [square]', 'Ours [triangle]');
axis([0 5 0 0.1]);

