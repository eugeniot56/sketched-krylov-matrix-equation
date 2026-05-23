function sketch_operator = setup_sketching_handle(problem_dim, sketch_dim)
%SETUP_SKETCHING_OPERATOR
%   Constructs a randomized sketching operator based on a subsampled
%   randomized discrete cosine transform (SRDCT).
%
%   INPUT:
%       problem_dim   - dimension of the original space
%       sketch_dim    - dimension of the sketching space
%
%   OUTPUT:
%       sketch_operator - function handle implementing the sketch
%
%   The sketching operator has the form
%
%       S = D * F * E,
%
%   where:
%       E is a random Rademacher diagonal matrix,
%       F is the discrete cosine transform (DCT),
%       D is a random row sampling operator.
%
%   The returned handle can be applied to vectors or matrices:
%
%       Y = sketch_operator(X)

    %% Input validation

    if problem_dim <= 0 || floor(problem_dim) ~= problem_dim
        error('problem_dim must be a positive integer.');
    end

    if sketch_dim <= 0 || floor(sketch_dim) ~= sketch_dim
        error('sketch_dim must be a positive integer.');
    end

    if sketch_dim > problem_dim
        error('sketch_dim must not exceed problem_dim.');
    end

    %% Randomized diagonal Rademacher matrix

    rng('default');

    random_signs = 2 * round(rand(problem_dim, 1)) - 1;

    E = spdiags(random_signs, 0, problem_dim, problem_dim);

    %% Random row sampling operator

    sampled_indices = randperm(problem_dim, sketch_dim);

    D = speye(problem_dim);
    D = D(sampled_indices, :);

    %% Scaling factor

    scaling = sqrt(sketch_dim / problem_dim);

    %% Sketching operator

    sketch_operator = @(X) D * dct(E * X) / scaling;

end