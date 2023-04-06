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

image_rows = 100;
image_cols = 100;
kernel_size = 5;
num_trials = 1000;

synth_exp = SyntheticExperiements(image_rows, image_cols);

legend_data = {'saddle', 'monkey saddle', 'opencv cb', 'opencv tri'};

[e1, e2, x] = synth_exp.runBlur(kernel_size, num_trials);
subplot(1,3,1); hold off;
plot(x, e1{1}, 'gs-', 'LineWidth', 1.5);  hold on;
plot(x, e1{2}, 'r^-', 'LineWidth', 1.5);
plot(x, e2{1}, 'ks--', 'LineWidth', 1.5);
plot(x, e2{2}, 'b^--', 'LineWidth', 1.5);
axis([0 x(end) 0.0 0.2]);
xlabel('Gaussian blur \sigma_n [px]');
ylabel('Mean Localization Error [px]');
set(gca, 'color', 'w'); set(gcf, 'color', 'w');
legend(legend_data{:}, 'Location', 'NorthEast');
set(gca, 'FontSize', 20);
grid on; grid minor;


[e1, e2, x] = synth_exp.runNoise(kernel_size, num_trials);
subplot(1,3,2); hold off;
plot(x, e1{1}, 'gs-', 'LineWidth', 1.5);  hold on;
plot(x, e1{2}, 'r^-', 'LineWidth', 1.5);
plot(x, e2{1}, 'ks--', 'LineWidth', 1.5);
plot(x, e2{2}, 'b^--', 'LineWidth', 1.5);
axis([0 x(end) 0.0 0.2]);
xlabel('Gaussian noise \sigma_n [%]');
%ylabel('Mean Localization Error [px]');
set(gca, 'color', 'w'); set(gcf, 'color', 'w');
%legend(legend_data{:}, 'Location', 'NorthWest');
set(gca, 'FontSize', 20);
grid on; grid minor;

[e1, e2, x] = synth_exp.runPerspective(kernel_size, num_trials);
subplot(1,3,3); hold off;
plot(x, e1{1}, 'gs-', 'LineWidth', 1.5);  hold on;
plot(x, e1{2}, 'r^-', 'LineWidth', 1.5);
plot(x, e2{1}, 'ks--', 'LineWidth', 1.5);
plot(x, e2{2}, 'b^--', 'LineWidth', 1.5);
axis([0 x(end) 0.0 0.2]);
xlabel('Perspective angle [deg]');
%ylabel('Mean Localization Error [px]');
set(gca, 'color', 'w'); set(gcf, 'color', 'w');
%legend(legend_data{:}, 'Location', 'NorthWest');
set(gca, 'FontSize', 20);
grid on; grid minor;
