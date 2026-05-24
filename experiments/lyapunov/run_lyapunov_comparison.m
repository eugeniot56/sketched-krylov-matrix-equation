% -------------------------------------------------------------------------
% Copyright (c) 2026, Eugenio Turchet
%
% This file is part of the project:
% "Sketching methods for generalized Lyapunov and Sylvester equations".
%
% This code is distributed under the BSD 3-Clause License.
% -------------------------------------------------------------------------

%% RUN_LYAPUNOV_COMPARISON
% Comparison between the full Krylov-Galerkin method and the sketched
% Krylov-Galerkin method for the generalized Lyapunov equation
%
%     A*X + X*A' + M*X*M' = b*b'.

clear;
clc;
close all;

format short e;
format compact;

%% Paths

this_file = mfilename('fullpath');
[this_folder, ~, ~] = fileparts(this_file);

repository_root = fullfile(this_folder, '..', '..');
src_folder = fullfile(repository_root, 'src');
external_folder = fullfile(repository_root, 'external');

results_folder = fullfile(repository_root, 'results');
figures_folder = fullfile(results_folder, 'figures', 'lyapunov');
tables_folder = fullfile(results_folder, 'tables');

if ~exist(figures_folder, 'dir')
    mkdir(figures_folder);
end

if ~exist(tables_folder, 'dir')
    mkdir(tables_folder);
end

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

fprintf('\n');
fprintf('Problem dimensions\n');
fprintf('--------------------------------------------------\n');
fprintf('Grid size        : %d x %d\n', grid_size, grid_size);
fprintf('Problem dimension: %d\n', problem_dim);
fprintf('--------------------------------------------------\n');

%% Solver parameters

max_iter = 300;
tolerance = 5e-6;

arnoldi_window_full = problem_dim;
arnoldi_window_sketch = 6;

residual_check = 1;
sketch_dim = 140;

fprintf('\n');
fprintf('Solver parameters\n');
fprintf('--------------------------------------------------\n');
fprintf('Maximum iterations       : %d\n', max_iter);
fprintf('Tolerance                : %.3e\n', tolerance);
fprintf('Full Arnoldi window      : %d\n', arnoldi_window_full);
fprintf('Sketched Arnoldi window  : %d\n', arnoldi_window_sketch);
fprintf('Sketch dimension         : %d\n', sketch_dim);
fprintf('Residual check frequency : %d\n', residual_check);
fprintf('--------------------------------------------------\n');

%% Full Krylov-Galerkin method

fprintf('\n');
fprintf('--------------------------------------------------\n');
fprintf('Running full Krylov-Galerkin method\n');
fprintf('--------------------------------------------------\n');

tic;

[V_full, Y_full, res_full, info_full] = solve_generalized_lyapunov_krylov( ...
    A, ...
    M, ...
    b, ...
    max_iter, ...
    tolerance, ...
    arnoldi_window_full, ...
    residual_check ...
);

time_full = toc;

%% Sketched Krylov-Galerkin method

fprintf('\n');
fprintf('--------------------------------------------------\n');
fprintf('Running sketched Krylov-Galerkin method\n');
fprintf('--------------------------------------------------\n');

tic;

[V_sketch, Y_sketch, res_sketch, info_sketch] = ...
    solve_generalized_lyapunov_sketch_krylov( ...
        A, ...
        M, ...
        b, ...
        max_iter, ...
        tolerance, ...
        arnoldi_window_sketch, ...
        sketch_dim, ...
        residual_check);

time_sketch = toc;

%% Memory usage and reduced dimensions

memory_V_full = whos('V_full');
memory_Y_full = whos('Y_full');

memory_full_MB = ...
    (memory_V_full.bytes + memory_Y_full.bytes) / 1024^2;

rank_full = rank(Y_full);
space_dim_full = info_full.subspace_dim;

memory_V_sketch = whos('V_sketch');
memory_Y_sketch = whos('Y_sketch');

memory_sketch_MB = ...
    (memory_V_sketch.bytes + memory_Y_sketch.bytes) / 1024^2;

rank_sketch = rank(Y_sketch);
space_dim_sketch = info_sketch.subspace_dim;

%% Residual evolution plot

figure;

semilogy( ...
    res_full, ...
    '-*', ...
    'LineWidth', 1.2);

hold on;

semilogy( ...
    res_sketch, ...
    '-*', ...
    'LineWidth', 1.2);

grid on;
grid minor;

xlabel('Iteration');
ylabel('Relative residual');
title('Residual evolution');

legend( ...
    'Full Krylov-Galerkin', ...
    'Sketched Krylov-Galerkin', ...
    'Location', 'best');

set(gca, ...
    'FontSize', 13, ...
    'LineWidth', 1);

exportgraphics( ...
    gcf, ...
    fullfile(figures_folder, 'residual_evolution.png'), ...
    'Resolution', 300);

%% Timing comparison plot

figure;

bar([time_full, time_sketch]);

set(gca, ...
    'XTickLabel', ...
    {'Full Krylov-Galerkin', 'Sketched Krylov-Galerkin'});

ylabel('Time (s)');
title('Timing comparison');

grid on;

set(gca, ...
    'FontSize', 13, ...
    'LineWidth', 1);

exportgraphics( ...
    gcf, ...
    fullfile(figures_folder, 'timing_comparison.png'), ...
    'Resolution', 300);

%% Save summary table

summary_table = table( ...
    ["Full Krylov-Galerkin"; "Sketched Krylov-Galerkin"], ...
    [info_full.iterations; info_sketch.iterations], ...
    [space_dim_full; space_dim_sketch], ...
    [rank_full; rank_sketch], ...
    [info_full.final_residual; info_sketch.final_residual], ...
    [time_full; time_sketch], ...
    [memory_full_MB; memory_sketch_MB], ...
    'VariableNames', ...
    {'Method', 'Iterations', 'SpaceDimension', 'Rank', ...
     'FinalResidual', 'TimeSeconds', 'MemoryMB'} ...
);

writetable( ...
    summary_table, ...
    fullfile(tables_folder, 'lyapunov_results.csv'));

%% Summary table

fprintf('\n');
fprintf('====================================================================================================\n');
fprintf('METHOD COMPARISON\n');
fprintf('====================================================================================================\n');

fprintf('\n');
fprintf('%-30s %-10s %-12s %-10s %-14s %-12s %-14s\n', ...
    'Method', ...
    'Iter.', ...
    'Space dim', ...
    'Rank', ...
    'Final Res.', ...
    'Time (s)', ...
    'Memory (MB)');

fprintf('----------------------------------------------------------------------------------------------------\n');

fprintf('%-30s %-10d %-12d %-10d %-14.3e %-12.4f %-14.4f\n', ...
    'Full Krylov-Galerkin', ...
    info_full.iterations, ...
    space_dim_full, ...
    rank_full, ...
    info_full.final_residual, ...
    time_full, ...
    memory_full_MB);

fprintf('%-30s %-10d %-12d %-10d %-14.3e %-12.4f %-14.4f\n', ...
    'Sketched Krylov-Galerkin', ...
    info_sketch.iterations, ...
    space_dim_sketch, ...
    rank_sketch, ...
    info_sketch.final_residual, ...
    time_sketch, ...
    memory_sketch_MB);

fprintf('====================================================================================================\n');

fprintf('\nNote: memory values are estimates based on the final main solution variables only.\n');
fprintf('They do not represent the peak memory usage during the whole computation.\n');

fprintf('\nResults saved in:\n');
fprintf('  Figures: %s\n', figures_folder);
fprintf('  Tables : %s\n', tables_folder);