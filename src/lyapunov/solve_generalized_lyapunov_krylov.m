% -------------------------------------------------------------------------
% Copyright (c) 2026, Eugenio Turchet
%
% This file is part of the project:
% "Sketching methods for generalized Lyapunov and Sylvester equations".
%
% This code is distributed under the BSD 3-Clause License.
% -------------------------------------------------------------------------



function [V, Y, res, info] = solve_generalized_lyapunov_krylov( ...
    A, M, b, max_iter, tolerance, arnoldi_window, residual_check)
%SOLVE_GENERALIZED_LYAPUNOV_KRYLOV
%   Krylov-Galerkin method for the generalized Lyapunov equation
%
%       A*X + X*A' + M*X*M' = b*b'.
%
%   The method computes an approximate solution of the form
%
%       X = V*Y*V',
%
%   where V is an orthonormal basis of a Krylov-type subspace and Y solves
%   the reduced projected equation.
%
%   INPUT:
%       A               - n-by-n matrix
%       M               - n-by-n matrix
%       b               - n-by-1 right-hand side vector
%       max_iter        - maximum number of outer iterations
%       tolerance       - stopping tolerance on the relative residual
%       arnoldi_window  - truncation parameter for truncated Arnoldi
%       residual_check  - frequency for reduced solve and residual check
%
%   OUTPUT:
%       V       - orthonormal basis of the projection subspace
%       Y       - reduced solution matrix
%       res     - relative residual history
%       info    - structure containing iteration information
%
%   REQUIRED EXTERNAL FUNCTIONS:
%       truncated_Arnoldi.m
%       pcg2.m

    %% Input validation

    if nargin < 7
        error('Not enough input arguments.');
    end

    [n, nA] = size(A);

    if n ~= nA
        error('A must be a square matrix.');
    end

    if ~isequal(size(M), [n, n])
        error('M must have the same size as A.');
    end

    if size(b, 1) ~= n || size(b, 2) ~= 1
        error('b must be an n-by-1 column vector.');
    end

    if max_iter <= 0 || floor(max_iter) ~= max_iter
        error('max_iter must be a positive integer.');
    end

    if tolerance <= 0
        error('tolerance must be positive.');
    end

    if arnoldi_window <= 0 || floor(arnoldi_window) ~= arnoldi_window
        error('arnoldi_window must be a positive integer.');
    end

    if residual_check <= 0 || floor(residual_check) ~= residual_check
        error('residual_check must be a positive integer.');
    end

    %% Initialization

    norm_b = norm(b);

    if norm_b == 0
        error('The right-hand side vector b must be nonzero.');
    end

    V = b / norm_b;

    rhs_norm = norm_b^2;

    Y = [];
    res = NaN(max_iter + 1, 1);
    res(1) = 1;

    total_time = 0;
    last_pcg_iter = [];

    fprintf('\n');
    fprintf(' iter    space dim    rel. residual        time (s)\n');
    fprintf('-----------------------------------------------------\n');

    %% Main iteration

    for iter = 1:max_iter

        iter_time = tic;

        % Generate candidate directions using A and M
        if iter > size(V, 2)
            warning('No more basis vectors available for expansion. Stopping.');
            break;
        end

        W = [A * V(:, iter), M * V(:, iter)];

        % Orthogonalize the new directions against the current basis
        current_dim = size(V, 2);
        W = truncated_Arnoldi(W, V, arnoldi_window, current_dim);

        % Compress the candidate directions using SVD
        [U, S, ~] = svd(W, 'econ');
        singular_values = diag(S);

        if isempty(singular_values) || singular_values(1) == 0
            warning('No new independent direction was generated. Stopping.');
            break;
        end

        numerical_rank = sum(singular_values / singular_values(1) > 1e-10);

        if numerical_rank == 0
            warning('The numerical rank of the new block is zero. Stopping.');
            break;
        end

        % Update the projection basis
        V = [V, U(:, 1:numerical_rank)];

        subspace_dim = size(V, 2);

        % Reduced matrices
        AV = A * V;
        MV = M * V;

        A_reduced = V' * AV;
        M_reduced = V' * MV;

        b_reduced = norm_b * eye(subspace_dim, 1);
        rhs_reduced = b_reduced * b_reduced';

        % Solve reduced problem and compute residual only every residual_check iterations
        if mod(iter, residual_check) == 0

            Y0 = zeros(subspace_dim, subspace_dim);

            pcg_tolerance = 1e-8;
            pcg_max_iter = subspace_dim^2;

            % Diagonal Kronecker-type preconditioner
            D = kron( ...
                spdiags(diag(A_reduced), 0, subspace_dim, subspace_dim), ...
                spdiags(diag(A_reduced), 0, subspace_dim, subspace_dim) ...
            );

            [y, last_pcg_iter] = pcg2( ...
                A_reduced, ...
                M_reduced, ...
                rhs_reduced, ...
                Y0, ...
                pcg_max_iter, ...
                pcg_tolerance, ...
                D ...
            );

            Y = reshape(y, subspace_dim, subspace_dim);

            rel_residual = compute_relative_residual( ...
                V, Y, AV, MV, b, rhs_norm ...
            );

            elapsed_time = toc(iter_time);
            total_time = total_time + elapsed_time;

            res(iter + 1) = rel_residual;

            fprintf('%5d    %9d    %13.6e    %10.4f\n', ...
                iter, subspace_dim, rel_residual, elapsed_time);

            if rel_residual < tolerance
                break;
            end

        end

    end

    %% Output information

    res = res(~isnan(res));

    info.iterations = iter;
    info.subspace_dim = size(V, 2);
    info.final_residual = res(end);
    info.total_time = total_time;
    info.last_pcg_iterations = last_pcg_iter;

end


function rel_residual = compute_relative_residual(V, Y, AV, MV, b, rhs_norm)
%COMPUTE_RELATIVE_RESIDUAL
%   Computes the relative Frobenius norm of the residual
%
%       R = A*V*Y*V' + V*Y*V'*A' + M*V*Y*V'*M' - b*b'
%
%   without forming the full approximate solution X = V*Y*V'.

    subspace_dim = size(V, 2);

    Z = [AV, V, MV, b];

    [~, S] = qr(Z, 0);

    Z11 = zeros(subspace_dim, subspace_dim);
    Z12 = zeros(subspace_dim, subspace_dim + 1);
    Z13 = zeros(subspace_dim, 1);

    last_row = [zeros(1, 3 * subspace_dim), -1];

    small_residual_matrix = [ ...
        Z11, Y,   Z12; ...
        Y,   Z11, Z12; ...
        Z11, Z11, Y, Z13; ...
        last_row ...
    ];

    projected_residual = S * small_residual_matrix * S';

    rel_residual = norm(projected_residual, 'fro') / rhs_norm;

end