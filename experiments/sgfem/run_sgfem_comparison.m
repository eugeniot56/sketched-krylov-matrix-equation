%% RUN_SGFEM_COMPARISON
% -------------------------------------------------------------------------
% Comparison between:
%
%   1) MultiRB method
%   2) Sketched Krylov-Galerkin method
%
% for generalized Sylvester equations arising from stochastic Galerkin FEM
% discretizations of elliptic PDEs.
%
% The generalized matrix equation has the form
%
%   K{1}*X*G{1} + K{2}*X*G{2} + ... + K{m}*X*G{m} = rhs1*rhs2'.
%
% -------------------------------------------------------------------------
% Reference:
%
% C.E. Powell, D. Silvester, V. Simoncini,
% "An Efficient Reduced Basis Solver for Stochastic Galerkin Matrix
% Equations",
% SIAM Journal on Scientific Computing,
% Vol. 39, No. 1, pp. A141--A163, 2017.
%
% -------------------------------------------------------------------------
% Original code:
% C.E. Powell, D. Silvester, V. Simoncini
%
% Modified by Eugenio Turchet, 2026.
% -------------------------------------------------------------------------

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
data_folder = fullfile(this_folder, 'data');

addpath(genpath(src_folder));
addpath(genpath(external_folder));

%% Select test problem

test_problem = input('\nTest problem 2, 3, 4, 5 or 6: ');

switch test_problem

    case 2
        data_file = 'TP_2.mat';

    case 3
        data_file = 'TP_3.mat';

    case 4
        data_file = 'TP_4.mat';

    case 5
        data_file = 'TP_5.mat';

  

    otherwise
        error('Invalid test problem. Choose 2, 3, 4, 5 or 6.');

end

data_path = fullfile(data_folder, data_file);

if ~isfile(data_path)
    error('Data file not found: %s', data_path);
end

load(data_path);

%% Problem dimensions

dim_K = size(K{1}, 1);
dim_G = size(G{1}, 1);
num_terms = length(K);

fprintf('\n');
fprintf('Problem dimensions\n');
fprintf('--------------------------------------------------\n');
fprintf('dim(K)          : %d\n', dim_K);
fprintf('dim(G)          : %d\n', dim_G);
fprintf('Number of terms : %d\n', num_terms);
fprintf('Total equations : %d\n', dim_K * dim_G);
fprintf('--------------------------------------------------\n');

%% Matrix preprocessing

permutation = symamd(Amean);

Amean_perm = Amean(permutation, permutation);

K_perm = cell(size(K));

for term = 1:num_terms
    K_perm{term} = K{term}(permutation, permutation);
end

rhs1 = fnew(1:dim_K);
rhs2 = eye(dim_G, 1);

rhs1_perm = rhs1(permutation);

%% Cholesky factorization of the permuted mean stiffness matrix

R = chol(Amean_perm, 'lower');

%% Shift selection

opts.tol = 1e-4;

K_perm{1} = (K_perm{1} + K_perm{1}') / 2;
K_perm{2} = (K_perm{2} + K_perm{2}') / 2;

emin2 = eigs(K_perm{2}, K_perm{1}, 1, 'SA', opts);
emax2 = eigs(K_perm{2}, K_perm{1}, 1, 'LA', opts);

mean_shift = mean([emin2, abs(emax2)]);

alphas = zeros(num_terms - 1, 1);
alphas(1) = 1 - mean_shift;

if num_terms > 2
    alphas(2:end) = 1;
end

%% Apply shifts

G_shifted = G;

for term = 2:num_terms

    K_perm{term} = (K_perm{term} + K_perm{term}') / 2;

    K_perm{term} = ...
        K_perm{term} + alphas(term - 1) * K_perm{1};

    G_shifted{1} = ...
        G_shifted{1} - alphas(term - 1) * G_shifted{term};

end

%% Parameter-free MultiRB setup

s_parameter = zeros(2, 1);
s_parameter(1) = 2 - alphas(1);
s_parameter(2) = 1;

%% MultiRB options

param.max_space_dim = 200;
param.period = 1;
param.rat_solve = '1';
param.res_method = '4';

fprintf('\n');
fprintf('--------------------------------------------------\n');
fprintf('Running MultiRB solver\n');
fprintf('--------------------------------------------------\n');

tic;

[X1, X2, dimV, final_err, avg_inner, error_vec, iv_vec] = ...
    MultiRB( ...
        K_perm, ...
        G_shifted, ...
        rhs1_perm, ...
        rhs2, ...
        Amean_perm, ...
        R, ...
        param, ...
        s_parameter);

time_multirb = toc;

X1(permutation, :) = X1;

%% Sketched Krylov-Galerkin setup

tolerance = 1e-6;
max_iter = 160;

arnoldi_window = 85;
sketch_dim = 300;

residual_check = 1;

fprintf('\n');
fprintf('--------------------------------------------------\n');
fprintf('Running sketched Krylov-Galerkin solver\n');
fprintf('--------------------------------------------------\n');

tic;

[V, Y, res, info_sketch] = ...
    solve_generalized_sylvester_sketch_krylov( ...
        K_perm, ...
        G_shifted, ...
        Amean_perm, ...
        R, ...
        rhs1_perm, ...
        rhs2, ...
        tolerance, ...
        max_iter, ...
        arnoldi_window, ...
        sketch_dim, ...
        residual_check);

time_sketch = toc;

%% Convergence history

figure;

semilogy(error_vec, '-o', 'LineWidth', 1.2);
hold on;

semilogy(res, '-o', 'LineWidth', 1.2);

grid on;

xlabel('Iteration');
ylabel('Relative error / residual');
title('Convergence history');

legend( ...
    'MultiRB', ...
    'Sketched Krylov-Galerkin', ...
    'Location', 'best');

%% Summary table

fprintf('\n');
fprintf('==============================================================\n');
fprintf('METHOD COMPARISON\n');
fprintf('==============================================================\n');

fprintf('\n');
fprintf('%-30s %-12s %-12s %-12s\n', ...
    'Method', ...
    'Iterations', ...
    'Final Res.', ...
    'Time (s)');

fprintf('--------------------------------------------------------------\n');

fprintf('%-30s %-12d %-12.3e %-12.4f\n', ...
    'MultiRB', ...
    length(error_vec), ...
    final_err, ...
    time_multirb);

fprintf('%-30s %-12d %-12.3e %-12.4f\n', ...
    'Sketched Krylov-Galerkin', ...
    info_sketch.iterations, ...
    info_sketch.final_residual, ...
    time_sketch);

fprintf('==============================================================\n');
