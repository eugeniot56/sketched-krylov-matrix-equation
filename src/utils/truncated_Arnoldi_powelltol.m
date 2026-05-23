function A = truncated_Arnoldi_powelltol(A, V, k, d, tolerance)
%TRUNCATED_ARNOLDI_POWELLTOL
%   Performs truncated Arnoldi orthogonalization using a Powell-type
%   tolerance parameter.
%
%   The columns of A are orthogonalized against the last k vectors of the
%   basis V, from index max(1, d-k+1) to d.
%
%   INPUT:
%       A           - block of vectors to orthogonalize
%       V           - current orthonormal basis
%       k           - truncation window size
%       d           - current basis dimension
%       tolerance   - tolerance parameter
%
%   OUTPUT:
%       A           - orthogonalized block of vectors

    start_index = max(1, d - k + 1);

    for i = start_index:d

        projection = V(:, i)' * A;
        A = A - V(:, i) * projection;

        projection = V(:, i)' * A;
        A = A - V(:, i) * projection;

    end

    column_norms = vecnorm(A);

    small_columns = column_norms <= tolerance;

    if any(small_columns)
        warning('Some columns have norm below the prescribed tolerance.');
    end

end
