% -------------------------------------------------------------------------
% Copyright (c) 2026, Eugenio Turchet
%
% This file is part of the project:
% "Sketching methods for generalized Lyapunov and Sylvester equations".
%
% This code is distributed under the BSD 3-Clause License.
% -------------------------------------------------------------------------



function [V, Y, res, info] = solve_generalized_lyapunov_sketch_krylov( ...
    A, M, b, max_iter, tolerance, arnoldi_window, sketch_dim, residual_check)
%SOLVE_GENERALIZED_LYAPUNOV_SKETCH_KRYLOV
%   Sketched Krylov-Galerkin method for the generalized Lyapunov equation
%
%       A*X + X*A' + M*X*M' = b*b'.
%
%   The method builds a Krylov-type basis V and uses a sketching operator
%   to define and solve a reduced projected problem.
%
%   REQUIRED EXTERNAL FUNCTIONS:
%       truncated_Arnoldi.m
%       setup_sketching_handle.m
%       qr_aggiornata.m

    %% Input validation

    if nargin < 8
        error('Not enough input arguments.');
    end

    [n, nA] = size(A);

    if n ~= nA
        error('A must be square.');
    end

    if ~isequal(size(M), [n, n])
        error('M must have the same size as A.');
    end

    if size(b, 1) ~= n || size(b, 2) ~= 1
        error('b must be an n-by-1 column vector.');
    end

    if norm(b) == 0
        error('b must be nonzero.');
    end

    %% Initialization

    norm_b = norm(b);

    V = b / norm_b;

    sketch = setup_sketching_handle(n, sketch_dim);

    Sb = sketch(b);
    Sb_normalized = Sb / norm_b;

    rhs_norm_sketch = norm(Sb)^2;

    Q_sketch = Sb_normalized / norm(Sb_normalized);
    R_sketch = norm(Sb_normalized);

    Vhat = V / R_sketch;

    AVhat = A * Vhat;
    MVhat = M * Vhat;

    S_AVhat = sketch(AVhat);
    S_MVhat = sketch(MVhat);

    A_reduced = Q_sketch' * S_AVhat;
    M_reduced = Q_sketch' * S_MVhat;

    Y = [];
    res = NaN(max_iter, 1);

    total_time = 0;
    last_gmres_flag = [];
    last_gmres_relres = [];
    last_gmres_iter = [];

    fprintf('\n');
    fprintf('Sketched Krylov-Galerkin method\n');
    fprintf('---------------------------------------------------------------\n');
    fprintf(' iter    space dim    rel. residual     time (s)    GMRES it   rank\n');
    fprintf('---------------------------------------------------------------\n');

    %% Main iteration

    for iter = 1:max_iter

        iter_time = tic;

        % Generate new Krylov directions
        W = [A * V(:, iter), M * V(:, iter)];

        current_dim = size(V, 2);
        W = truncated_Arnoldi(W, V, arnoldi_window, current_dim);

        % Orthonormalize new directions in Euclidean norm
        [V_new, ~] = qr(W, 0);

        if isempty(V_new)
            warning('No new independent direction generated. Stopping.');
            break;
        end

        V = [V, V_new];

        % Sketch the new Euclidean basis vectors
        S_V_new = sketch(V_new);

        % Update QR factorization of S*V
        [Q_sketch, R_sketch] = update_qr_factorization(S_V_new, Q_sketch, R_sketch);

        % Construct basis orthonormal with respect to the sketched inner product
        Vhat = V / R_sketch;

        % Apply operators only to the last added sketched-orthonormal vectors
        num_new_vectors = size(V_new, 2);

        Vhat_new = Vhat(:, end - num_new_vectors + 1:end);

        AVhat_new = A * Vhat_new;
        MVhat_new = M * Vhat_new;

        S_AVhat = [S_AVhat, sketch(AVhat_new)];
        S_MVhat = [S_MVhat, sketch(MVhat_new)];

        % Reduced sketched matrices
        A_reduced = Q_sketch' * S_AVhat;
        M_reduced = Q_sketch' * S_MVhat;

        reduced_dim = size(A_reduced, 1);

        b_reduced = norm(Sb) * eye(reduced_dim, 1);
        rhs_reduced = b_reduced * b_reduced';

        if mod(iter, residual_check) == 0

            reduced_operator = @(y) apply_reduced_operator( ...
                y, A_reduced, M_reduced, reduced_dim);

            gmres_tolerance = tolerance / 100;
            gmres_max_iter = 500;

            [y, flag, relres, gmres_iter] = gmres( ...
                reduced_operator, ...
                rhs_reduced(:), ...
                [], ...
                gmres_tolerance, ...
                gmres_max_iter ...
            );

            Y = reshape(y, reduced_dim, reduced_dim);

            rel_residual = compute_sketched_relative_residual( ...
                S_AVhat, Q_sketch, S_MVhat, Sb, Y, rhs_norm_sketch);

            elapsed_time = toc(iter_time);
            total_time = total_time + elapsed_time;

            res(iter) = rel_residual;

            last_gmres_flag = flag;
            last_gmres_relres = relres;
            last_gmres_iter = gmres_iter;

            fprintf('%5d    %9d    %13.6e   %10.4f   %8d   %4d\n', ...
                iter, reduced_dim, rel_residual, elapsed_time, ...
                gmres_iter(end), rank(R_sketch));

            if rel_residual < tolerance
                break;
            end

        end

    end

    %% Post-processing

    res = res(~isnan(res));

    info.iterations = iter;
    info.subspace_dim = size(V, 2);
    info.reduced_dim = size(Y, 1);
    info.final_residual = res(end);
    info.total_time = total_time;
    info.gmres_flag = last_gmres_flag;
    info.gmres_relres = last_gmres_relres;
    info.gmres_iterations = last_gmres_iter;
    info.rank_R_sketch = rank(R_sketch);

end


function y = apply_reduced_operator(x, A_reduced, M_reduced, reduced_dim)
%APPLY_REDUCED_OPERATOR
%   Applies the vectorized reduced Lyapunov operator
%
%       L(Y) = A_reduced*Y + Y*A_reduced' + M_reduced*Y*M_reduced'.

    X = reshape(x, reduced_dim, reduced_dim);

    Y = A_reduced * X ...
        + X * A_reduced' ...
        + M_reduced * X * M_reduced';

    y = Y(:);

end


function rel_residual = compute_sketched_relative_residual( ...
    S_AVhat, Q_sketch, S_MVhat, Sb, Y, rhs_norm_sketch)
%COMPUTE_SKETCHED_RELATIVE_RESIDUAL
%   Computes the relative residual in the sketched space.

    reduced_dim = size(Y, 1);

    Z = [S_AVhat, Q_sketch, S_MVhat, Sb];

    O11 = sparse(reduced_dim, reduced_dim);
    O12 = sparse(reduced_dim, reduced_dim + 1);
    O13 = sparse(reduced_dim, 1);

    last_row = [sparse(1, 3 * reduced_dim), -1];

    small_residual_matrix = [ ...
        O11, Y,   O12; ...
        Y,   O11, O12; ...
        O11, O11, Y, O13; ...
        last_row ...
    ];

    sketched_residual = Z * small_residual_matrix * Z';

    rel_residual = norm(sketched_residual, 'fro') / rhs_norm_sketch;

end