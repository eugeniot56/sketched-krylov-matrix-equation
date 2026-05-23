function A = truncated_Arnoldi(A, V, k, d)
%TRUNCATED_ARNOLDI
%   Performs truncated Arnoldi orthogonalization of the columns of A
%   against the last k vectors of the orthonormal basis V.
%
%   INPUT:
%       A   - block of vectors to orthogonalize
%       V   - orthonormal basis matrix
%       k   - truncation window size
%       d   - current dimension of the basis V
%
%   OUTPUT:
%       A   - orthogonalized block of vectors
%
%   The orthogonalization is performed using two passes of modified
%   Gram-Schmidt for improved numerical stability.

    start_index = max(1, d - k + 1);

    for i = start_index:d

        projection = V(:, i)' * A;
        A = A - V(:, i) * projection;

        projection = V(:, i)' * A;
        A = A - V(:, i) * projection;

    end

end