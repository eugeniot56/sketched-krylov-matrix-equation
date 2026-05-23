% -------------------------------------------------------------------------
% Copyright (c) 2026, Eugenio Turchet
%
% This file is part of the project:
% "Sketching methods for generalized Lyapunov and Sylvester equations".
%
% This code is distributed under the BSD 3-Clause License.
% -------------------------------------------------------------------------



%% TEST_GENERALIZED_LYAPUNOV_KRYLOV
% Comparison between the full Krylov-Galerkin method and the sketched
% Krylov-Galerkin method for the generalized Lyapunov equation
%
%     A*X + X*A' + M*X*M' = b*b'.

clear;
clc;
close all;

this_file = mfilename('fullpath');
[this_folder, ~, ~] = fileparts(this_file);

repository_root = fullfile(this_folder, '..', '..');
src_folder = fullfile(repository_root, 'src');
external_folder = fullfile(repository_root, 'external');

addpath(genpath(src_folder));
addpath(genpath(external_folder));

%% Problem setup

grid_size = 900;
problem_dim = grid_size^2;

A = poisson2d(grid_size, grid_size) + 10 * speye(problem_dim);

e = ones(problem_dim, 1);
random_perturbation = rand(problem_dim, 1);

M = spdiags( ...
    [-e, 3 * e + 2 * random_perturbation, -e], ...
    -1:1, ...
    problem_dim, ...
    problem_dim ...
);

b = ones(problem_dim, 1);
%b = randn(problem_dim,1);
%% Solver parameters

max_iter = 300;
tolerance = 5e-6;

arnoldi_window_full = problem_dim;
arnoldi_window_sketch = 6;

residual_check = 1;
sketch_dim = 140;

%% Full Krylov-Galerkin method

fprintf('\nRunning full Krylov-Galerkin method...\n');

[V_full, Y_full, res_full, info_full] = solve_generalized_lyapunov_krylov( ...
    A, ...
    M, ...
    b, ...
    max_iter, ...
    tolerance, ...
    arnoldi_window_full, ...
    residual_check ...
);

%% Sketched Krylov-Galerkin method

fprintf('\nRunning sketched Krylov-Galerkin method...\n');

[V_sketch, Y_sketch, res_sketch, info_sketch] = solve_generalized_lyapunov_sketch_krylov( ...
    A, ...
    M, ...
    b, ...
    max_iter, ...
    tolerance, ...
    arnoldi_window_sketch, ...
    sketch_dim, ...
    residual_check ...
);

%% Plot residual histories

figure;
semilogy(res_full, '-*', 'LineWidth', 1.2);
hold on;
semilogy(res_sketch, '-*', 'LineWidth', 1.2);
grid on;

xlabel('Iteration');
ylabel('Relative residual');
title('Residual history');

legend( ...
    'Full Krylov-Galerkin', ...
    'Sketched Krylov-Galerkin', ...
    'Location', 'best' ...
);

%% Print summary

fprintf('\nSummary\n');
fprintf('--------------------------------------------------\n');

fprintf('Full Krylov-Galerkin:\n');
fprintf('  Iterations:      %d\n', info_full.iterations);
fprintf('  Subspace dim:    %d\n', info_full.subspace_dim);
fprintf('  Final residual:  %.3e\n', info_full.final_residual);
fprintf('  Total time:      %.4f s\n\n', info_full.total_time);

fprintf('Sketched Krylov-Galerkin:\n');
fprintf('  Iterations:      %d\n', info_sketch.iterations);
fprintf('  Subspace dim:    %d\n', info_sketch.subspace_dim);
fprintf('  Reduced dim:     %d\n', info_sketch.reduced_dim);
fprintf('  Final residual:  %.3e\n', info_sketch.final_residual);
fprintf('  Total time:      %.4f s\n', info_sketch.total_time);