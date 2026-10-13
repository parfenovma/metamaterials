module OperatorO3Statistics

using LinearAlgebra

export precision_dictionary,
       matched_output_covariance,
       matched_mean,
       covariance_factor,
       normalized_scores,
       wilson_upper_bound,
       crossing_amplitude

"""Apply `(I + B*B')^-1` to the nominal dictionary without forming it."""
function precision_dictionary(dictionary::AbstractMatrix, clutter::AbstractMatrix)
    size(clutter, 1) == size(dictionary, 1) ||
        throw(DimensionMismatch("clutter and dictionary rows must agree"))
    size(clutter, 2) == 0 && return Matrix{ComplexF64}(dictionary)
    small_system = Hermitian(
        Matrix{ComplexF64}(I, size(clutter, 2), size(clutter, 2)) +
        clutter' * clutter,
    )
    dictionary - clutter * (small_system \ (clutter' * dictionary))
end

"""Covariance of nominal matched-bank outputs under an unseen clutter basis."""
function matched_output_covariance(
    precision_applied_dictionary::AbstractMatrix,
    test_clutter::AbstractMatrix,
)
    size(test_clutter, 1) == size(precision_applied_dictionary, 1) ||
        throw(DimensionMismatch("test clutter and dictionary rows must agree"))
    base = precision_applied_dictionary' * precision_applied_dictionary
    if size(test_clutter, 2) > 0
        coupling = test_clutter' * precision_applied_dictionary
        base += coupling' * coupling
    end
    Matrix(Hermitian((base + base') / 2))
end

"""Mean matched-bank output for every true-signal column."""
function matched_mean(
    precision_applied_dictionary::AbstractMatrix,
    test_dictionary::AbstractMatrix,
)
    size(precision_applied_dictionary, 1) == size(test_dictionary, 1) ||
        throw(DimensionMismatch("nominal and test dictionary rows must agree"))
    precision_applied_dictionary' * test_dictionary
end

"""Stable lower factor for a Hermitian positive-semidefinite covariance."""
function covariance_factor(covariance::AbstractMatrix; relative_jitter=1.0e-10)
    size(covariance, 1) == size(covariance, 2) ||
        throw(DimensionMismatch("covariance must be square"))
    relative_jitter >= 0 || throw(ArgumentError("jitter must be nonnegative"))
    hermitian = Hermitian((covariance + covariance') / 2)
    scale = max(maximum(real, diag(hermitian)), eps(Float64))
    regularized = Hermitian(
        Matrix(hermitian) + relative_jitter * scale *
        Matrix{Float64}(I, size(covariance, 1), size(covariance, 1)),
    )
    cholesky(regularized).L
end

"""Maximum normalized GLRT score and template index for every trial column."""
function normalized_scores(outputs::AbstractMatrix, normalization::AbstractVector)
    size(outputs, 1) == length(normalization) ||
        throw(DimensionMismatch("one normalization value per template is required"))
    all(>(0), normalization) ||
        throw(ArgumentError("template normalizations must be positive"))
    maxima = Vector{Float64}(undef, size(outputs, 2))
    indices = Vector{Int}(undef, size(outputs, 2))
    for trial in axes(outputs, 2)
        best_score = -Inf
        best_index = 0
        for template in axes(outputs, 1)
            score = abs2(outputs[template, trial]) / normalization[template]
            if score > best_score
                best_score = score
                best_index = template
            end
        end
        maxima[trial] = best_score
        indices[trial] = best_index
    end
    maxima, indices
end

"""One-sided Wilson upper confidence bound for a binomial proportion."""
function wilson_upper_bound(successes::Integer, trials::Integer; z=1.6448536269514722)
    0 <= successes <= trials || throw(ArgumentError("invalid binomial counts"))
    trials > 0 || throw(ArgumentError("at least one trial is required"))
    z > 0 || throw(ArgumentError("z must be positive"))
    probability = successes / trials
    denominator = 1 + z^2 / trials
    center = probability + z^2 / (2trials)
    spread = z * sqrt(
        probability * (1 - probability) / trials + z^2 / (4trials^2),
    )
    (center + spread) / denominator
end

"""Log-linear first crossing of a monotone probability curve."""
function crossing_amplitude(amplitudes, probabilities, target)
    length(amplitudes) == length(probabilities) ||
        throw(DimensionMismatch("amplitudes and probabilities must agree"))
    !isempty(amplitudes) || throw(ArgumentError("curve must not be empty"))
    all(>(0), amplitudes) || throw(ArgumentError("amplitudes must be positive"))
    0 <= target <= 1 || throw(ArgumentError("target must be a probability"))
    monotone = accumulate(max, Float64.(probabilities))
    first_index = findfirst(>=(target), monotone)
    isnothing(first_index) && return Inf
    first_index == 1 && return Float64(first(amplitudes))
    lower_index = first_index - 1
    lower_probability = monotone[lower_index]
    upper_probability = monotone[first_index]
    upper_probability <= lower_probability + eps(Float64) &&
        return Float64(amplitudes[first_index])
    fraction = (target - lower_probability) / (upper_probability - lower_probability)
    exp(
        log(Float64(amplitudes[lower_index])) + fraction *
        (log(Float64(amplitudes[first_index])) - log(Float64(amplitudes[lower_index]))),
    )
end

end
