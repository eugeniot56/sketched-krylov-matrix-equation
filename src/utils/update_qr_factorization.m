function [Q, R] = update_qr_factorization(V_new, Q, R)
%UPDATE_QR_FACTORIZATION
%   Updates the thin QR factorization after appending new columns.
%
%   Given an existing factorization
%
%       V_old = Q * R,
%
%   this function computes the QR factorization of
%
%       [V_old, V_new]
%
%   without recomputing it from scratch.
%
%   INPUT:
%       V_new   - new block of columns to append
%       Q       - current orthonormal factor
%       R       - current upper triangular factor
%
%   OUTPUT:
%       Q       - updated orthonormal factor
%       R       - updated upper triangular factor

    num_new_vectors = size(V_new, 2);

    projection = Q' * V_new;
    orthogonal_part = V_new - Q * projection;

    [Q_new, R_new] = qr(orthogonal_part, 0);

    Q = [Q, Q_new];

    R = [ ...
        R, projection; ...
        zeros(num_new_vectors, size(R, 2)), R_new ...
    ];

end