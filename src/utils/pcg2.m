% -------------------------------------------------------------------------
% Copyright (c) 2026, Eugenio Turchet
%
% This file is part of the project:
% "Sketching methods for generalized Lyapunov and Sylvester equations".
%
% This code is distributed under the BSD 3-Clause License.
% -------------------------------------------------------------------------

function [x, iter] = pcg2( ...
    A_reduced, M_reduced, b, X0, max_iter, tolerance, D, verbose)

    if nargin < 8
        verbose = false;
    end

    [n_rows, n_cols] = size(X0);

    residual_matrix = ...
        b ...
        - A_reduced * X0 ...
        - X0 * A_reduced' ...
        - M_reduced * X0 * M_reduced';

    x = X0(:);
    r = residual_matrix(:);

    initial_residual = norm(r);

    if initial_residual == 0
        iter = 0;
        return;
    end

    residual = initial_residual;
    iter = 0;

    while (residual / initial_residual > tolerance) && (iter < max_iter)

        z = D \ r;

        iter = iter + 1;

        gamma = r' * z;

        if iter == 1
            p = z;
        else
            beta = gamma / gamma_old;
            p = z + beta * p;
        end

        P = reshape(p, n_rows, n_cols);

        W = ...
            A_reduced * P ...
            + P * A_reduced' ...
            + M_reduced * P * M_reduced';

        w = W(:);

        delta = p' * w;

        if abs(delta) < eps
            warning('PCG breakdown: p''*A*p is close to zero.');
            break;
        end

        alpha = gamma / delta;

        x = x + alpha * p;
        r = r - alpha * w;

        gamma_old = gamma;
        residual = norm(r);

    end

    if verbose
        fprintf('PCG iterations: %d   relative residual: %.3e\n', ...
            iter, residual / initial_residual);
    end

end