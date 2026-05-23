% -------------------------------------------------------------------------
% Copyright (c) 2026, Eugenio Turchet
%
% This file is part of the project:
% "Sketching methods for generalized Lyapunov and Sylvester equations".
%
% This code is distributed under the BSD 3-Clause License.
% -------------------------------------------------------------------------


function [V, Y, res, info] = solve_generalized_sylvester_sketch_krylov( ...
    A, K, M, M1, f1, f2, tolerance, max_iter, arnoldi_window, sketch_dim, residual_check)
%SOLVE_GENERALIZED_SYLVESTER_SKETCH_KRYLOV
%   Sketched Krylov-Galerkin method for generalized Sylvester-type equations
%
%       A{1}*X*K{1} + A{2}*X*K{2} + ... + A{m}*X*K{m} = f1*f2'.
%
%   The method builds a Krylov-type basis and solves a sketched projected
%   problem. The vector f1 is first preconditioned using M1.
%
%   REQUIRED EXTERNAL FUNCTIONS:
%       setup_sketching_operator.m
%       update_qr_factorization.m
%       truncated_Arnoldi_powerToll.m
%       somma.m
%       residual_norm_Powertoll.m
%       bicgstabl.m

    %% Input validation

    if nargin < 11
        error('Not enough input arguments.');
    end

    num_terms = length(A);

    if length(K) ~= num_terms
        error('A and K must have the same number of cells.');
    end

    [n, ~] = size(A{1});

    %% Initial vector

    f1 = M1' \ (M1 \ f1);
    norm_f1 = norm(f1);

    if norm_f1 == 0
        error('The preconditioned vector f1 must be nonzero.');
    end

    V = f1 / norm_f1;

    %% Sketching setup

    sketch = setup_sketching_handle(n, sketch_dim);

    Sf1_full = sketch(f1);
    Sf1 = Sf1_full / norm_f1;

    residual_norm0 = norm(Sf1_full * f2', 'fro');

    Q_sketch = Sf1 / norm(Sf1);
    R_sketch = norm(Sf1);

    Vhat = V / R_sketch;

    F_sketch = Sf1_full * f2';

    %% Initial reduced matrices

    S_A = cell(num_terms, 1);
    A_reduced = cell(num_terms, 1);

    S_A{1} = sketch(Vhat);
    A_reduced{1} = Q_sketch' * S_A{1};

    for term = 2:num_terms
        Z = M1 \ (A{term} * Vhat);
        AZ = M1' \ Z;

        S_A{term} = sketch(AZ);
        A_reduced{term} = Q_sketch' * S_A{term};
    end

    %% Iteration setup

    Y = [];
    Y_previous = [];

    res = NaN(max_iter + 1, 1);
    res(1) = 1;

    total_time = 0;

    last_solver_flag = [];
    last_solver_relres = [];
    last_solver_iter = [];

    fprintf('\n');
    fprintf('Sketched Krylov-Galerkin method for generalized Sylvester equation\n');
    fprintf('------------------------------------------------------------------\n');
    fprintf(' iter    space dim    rel. residual        time (s)\n');
    fprintf('------------------------------------------------------------------\n');

    %% Main iteration

    for iter = 1:max_iter

        iter_time = tic;

        % Generate new candidate directions
        W = [];

        for term = 2:num_terms
            Z = M1 \ (A{term} * V(:, iter));
            W = [W, M1' \ Z];
        end

        current_dim = size(V, 2);

        W = truncated_Arnoldi_powelltol( ...
            W, V, arnoldi_window, current_dim, tolerance);

        % Compress the new block using SVD
        [U, S, ~] = svd(W, 'econ');
        singular_values = diag(S);

        if isempty(singular_values) || singular_values(1) == 0
            warning('No new independent direction was generated. Stopping.');
            break;
        end

        cumulative_energy = cumsum(singular_values) / sum(singular_values);
        num_new_vectors = find(cumulative_energy >= 0.99, 1);

        V_new = U(:, 1:num_new_vectors);

        old_dim = size(V, 2);
        V = [V, V_new];

        %% Sketching update

        S_V_new = sketch(V_new);

        [Q_sketch, R_sketch] = update_qr_factorization( ...
            S_V_new, Q_sketch, R_sketch);

        Vhat = V / R_sketch;

        %% Update sketched operator blocks

        Vhat_new = Vhat(:, old_dim + 1:old_dim + num_new_vectors);

        S_A{1}(:, old_dim + 1:old_dim + num_new_vectors) = sketch(Vhat_new);
        A_reduced{1} = Q_sketch' * S_A{1};

        for term = 2:num_terms
            Z = M1 \ (A{term} * Vhat_new);
            AZ = M1' \ Z;

            S_A{term}(:, old_dim + 1:old_dim + num_new_vectors) = sketch(AZ);
            A_reduced{term} = Q_sketch' * S_A{term};
        end

        reduced_rows = size(A_reduced{1}, 1);
        reduced_cols = size(K{1}, 1);

        f_reduced = norm(Sf1_full) * eye(reduced_rows, 1);
        rhs_reduced = f_reduced * f2';

        %% Solve reduced problem periodically

        if mod(iter, residual_check) == 0

            reduced_operator = @(y) reshape( ...
                somma(A_reduced, K, reshape(y, reduced_rows, reduced_cols)), ...
                reduced_rows * reduced_cols, ...
                1 ...
            );

            preconditioner = @(x) reshape( ...
                pinv(A_reduced{1}) * reshape(x, reduced_rows, reduced_cols), ...
                reduced_rows * reduced_cols, ...
                1 ...
            );

            Y0 = zeros(reduced_rows, reduced_cols);

            if ~isempty(Y_previous)
                [old_rows, old_cols] = size(Y_previous);
                Y0(1:old_rows, 1:old_cols) = Y_previous;
            end

            [y, flag, relres, solver_iter] = bicgstabl( ...
                reduced_operator, ...
                rhs_reduced(:), ...
                tolerance / 100, ...
                100, ...
                preconditioner, ...
                [], ...
                Y0(:) ...
            );

            Y = reshape(y, reduced_rows, reduced_cols);
            Y_previous = Y;

            residual_value = residual_norm_powelltol(S_A, K, Y, F_sketch);
            rel_residual = residual_value / residual_norm0;

            elapsed_time = toc(iter_time);
            total_time = total_time + elapsed_time;

            res(iter + 1) = rel_residual;

            last_solver_flag = flag;
            last_solver_relres = relres;
            last_solver_iter = solver_iter;

            fprintf('%5d    %9d    %13.6e    %10.4f\n', ...
                iter, reduced_rows, rel_residual, elapsed_time);

            if rel_residual < tolerance
                break;
            end

        end

    end

    %% Output information

    res = res(~isnan(res));

    info.iterations = iter;
    info.subspace_dim = size(V, 2);
    info.reduced_rows = size(Y, 1);
    info.reduced_cols = size(Y, 2);
    info.final_residual = res(end);
    info.total_time = total_time;
    info.rank_basis = rank(V);

    info.last_solver_flag = last_solver_flag;
    info.last_solver_relres = last_solver_relres;
    info.last_solver_iterations = last_solver_iter;

end

function residual_norm = residual_norm_powelltol(A, K, Y, F)
%RESIDUAL_NORM_POWELLTOL
%   Computes the Frobenius norm of the residual associated with the
%   generalized matrix equation
%
%       A{1}*Y*K{1} + A{2}*Y*K{2} + ... + A{m}*Y*K{m} = F.
%
%   INPUT:
%       A   - cell array containing the left operators
%       K   - cell array containing the right operators
%       Y   - reduced solution matrix
%       F   - right-hand side matrix
%
%   OUTPUT:
%       residual_norm - Frobenius norm of the residual

residual_matrix = somma(A, K, Y) - F;

residual_norm = norm(residual_matrix, 'fro');

end

function Z = somma(K, G, Y)
%SOMMA
%   Computes the matrix sum
%
%       Z = K{1}*Y*G{1} + K{2}*Y*G{2} + ... + K{m}*Y*G{m}.
%
%   INPUT:
%       K   - cell array containing left matrices
%       G   - cell array containing right matrices
%       Y   - matrix to be multiplied
%
%   OUTPUT:
%       Z   - resulting matrix sum

    num_terms_K = length(K);
    num_terms_G = length(G);

    if num_terms_K ~= num_terms_G
        error('K and G must contain the same number of matrices.');
    end

    [num_rows, ~] = size(K{1});
    [~, num_cols] = size(G{1});

    Z = zeros(num_rows, num_cols);

    for term = 1:num_terms_K
        Z = Z + K{term} * Y * G{term};
    end

end