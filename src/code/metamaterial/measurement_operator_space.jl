module MeasurementOperatorSpace

using LinearAlgebra
using Statistics

export OperatorArchitecture,
       OperatorMetrics,
       architecture_space,
       validate,
       whiten_dictionary,
       normalized_gram,
       operator_metrics

"""Architecture-level point considered before any geometry is synthesized."""
Base.@kwdef struct OperatorArchitecture
    id::Symbol
    electronic_channels::Int
    acquisitions::Int
    mechanical_positions::Int = 1
    physical_states::Int = 1
    txrx_constraint::Symbol
    code_family::Symbol
    diversity::Tuple{Vararg{Symbol}} = ()
    lift_status::Symbol = :plane_strain_surrogate
    role::Symbol = :candidate
end

const TXRX_CONSTRAINTS = (:common_monostatic, :reciprocal_matched, :independent)
const LIFT_STATUSES = (
    :translation_kz,
    :axisymmetric_m,
    :group_sector,
    :plane_strain_surrogate,
    :full_3d_required,
)
const ROLES = (:baseline, :candidate, :bound)

function validate(architecture::OperatorArchitecture)
    architecture.electronic_channels > 0 ||
        throw(ArgumentError("electronic channel count must be positive"))
    architecture.acquisitions > 0 ||
        throw(ArgumentError("acquisition count must be positive"))
    architecture.mechanical_positions > 0 ||
        throw(ArgumentError("mechanical position count must be positive"))
    architecture.physical_states > 0 ||
        throw(ArgumentError("physical state count must be positive"))
    architecture.txrx_constraint in TXRX_CONSTRAINTS ||
        throw(ArgumentError("unsupported Tx/Rx constraint: $(architecture.txrx_constraint)"))
    architecture.lift_status in LIFT_STATUSES ||
        throw(ArgumentError("unsupported 2D-to-3D lift status: $(architecture.lift_status)"))
    architecture.role in ROLES ||
        throw(ArgumentError("unsupported architecture role: $(architecture.role)"))
    architecture
end

"""Frozen architecture classes for the first operator-space campaign."""
function architecture_space()
    architectures = OperatorArchitecture[
        OperatorArchitecture(
            id=:a_scan,
            electronic_channels=1,
            acquisitions=1,
            txrx_constraint=:common_monostatic,
            code_family=:uniform,
            diversity=(:time, :frequency),
            role=:baseline,
        ),
        OperatorArchitecture(
            id=:fixed_focus,
            electronic_channels=1,
            acquisitions=1,
            txrx_constraint=:reciprocal_matched,
            code_family=:fixed_focus,
            diversity=(:frequency,),
        ),
        OperatorArchitecture(
            id=:coded_broadband,
            electronic_channels=1,
            acquisitions=1,
            txrx_constraint=:common_monostatic,
            code_family=:calibrated_speckle,
            diversity=(:frequency, :spatial_fingerprint),
        ),
        OperatorArchitecture(
            id=:time_encoded,
            electronic_channels=1,
            acquisitions=1,
            txrx_constraint=:common_monostatic,
            code_family=:spatial_to_delay,
            diversity=(:time_slot, :frequency),
        ),
        OperatorArchitecture(
            id=:modal_encoded,
            electronic_channels=2,
            acquisitions=2,
            txrx_constraint=:reciprocal_matched,
            code_family=:mode_selective,
            diversity=(:pressure, :shear, :symmetry_sector),
        ),
        OperatorArchitecture(
            id=:few_port,
            electronic_channels=4,
            acquisitions=4,
            txrx_constraint=:reciprocal_matched,
            code_family=:fixed_passive_codebook,
            diversity=(:port, :frequency),
        ),
        OperatorArchitecture(
            id=:mechanical_index,
            electronic_channels=1,
            acquisitions=4,
            mechanical_positions=4,
            txrx_constraint=:common_monostatic,
            code_family=:translated_aperture,
            diversity=(:position, :frequency),
        ),
        OperatorArchitecture(
            id=:replaceable_mask,
            electronic_channels=1,
            acquisitions=4,
            physical_states=4,
            txrx_constraint=:common_monostatic,
            code_family=:physical_codebook,
            diversity=(:state, :frequency),
        ),
        OperatorArchitecture(
            id=:reconfigurable_state,
            electronic_channels=1,
            acquisitions=4,
            physical_states=4,
            txrx_constraint=:reciprocal_matched,
            code_family=:mechanically_reconfigured_surface,
            diversity=(:state, :surface_kinematics, :frequency),
        ),
        OperatorArchitecture(
            id=:ideal_k_code,
            electronic_channels=4,
            acquisitions=4,
            txrx_constraint=:independent,
            code_family=:unconstrained_complex,
            diversity=(:abstract_code,),
            role=:bound,
        ),
        OperatorArchitecture(
            id=:paut_oracle,
            electronic_channels=16,
            acquisitions=256,
            txrx_constraint=:independent,
            code_family=:full_matrix_capture,
            diversity=(:transmit_element, :receive_element, :frequency),
            role=:bound,
        ),
    ]
    validate.(architectures)
end

struct OperatorMetrics
    measurement_count::Int
    hypothesis_count::Int
    minimum_detectability::Float64
    median_detectability::Float64
    maximum_pair_coherence::Float64
    stable_rank::Float64
    numerical_rank::Int
    information_logdet_nats::Float64
end

function covariance_factor(noise_covariance::AbstractMatrix, measurement_count::Int)
    size(noise_covariance) == (measurement_count, measurement_count) ||
        throw(DimensionMismatch("noise covariance must match dictionary rows"))
    hermitian = Hermitian((noise_covariance + noise_covariance') / 2)
    isposdef(hermitian) || throw(ArgumentError("noise covariance must be positive definite"))
    cholesky(hermitian)
end

"""Return `C_n^(-1/2) A`, with hypotheses stored in columns."""
function whiten_dictionary(
    dictionary::AbstractMatrix;
    noise_covariance::AbstractMatrix=Matrix{Float64}(
        I, size(dictionary, 1), size(dictionary, 1),
    ),
)
    measurement_count, hypothesis_count = size(dictionary)
    measurement_count > 0 || throw(ArgumentError("dictionary must have at least one row"))
    hypothesis_count > 0 || throw(ArgumentError("dictionary must have at least one column"))
    factor = covariance_factor(noise_covariance, measurement_count)
    factor.L \ dictionary
end

"""Absolute normalized Gram matrix of a whitened defect dictionary."""
function normalized_gram(whitened_dictionary::AbstractMatrix)
    norms = vec(sqrt.(sum(abs2, whitened_dictionary; dims=1)))
    all(>(0), norms) || throw(ArgumentError("every defect hypothesis must have nonzero response"))
    normalized = whitened_dictionary ./ reshape(norms, 1, :)
    abs.(normalized' * normalized)
end

function information_logdet(whitened_dictionary::AbstractMatrix, signal_prior::AbstractMatrix)
    hypothesis_count = size(whitened_dictionary, 2)
    size(signal_prior) == (hypothesis_count, hypothesis_count) ||
        throw(DimensionMismatch("signal prior must match dictionary columns"))
    prior = Hermitian((signal_prior + signal_prior') / 2)
    isposdef(prior) || throw(ArgumentError("signal prior must be positive definite"))
    prior_factor = cholesky(prior)
    transformed = whitened_dictionary * prior_factor.L
    information = Hermitian(
        Matrix{eltype(whitened_dictionary)}(I, hypothesis_count, hypothesis_count) +
        transformed' * transformed,
    )
    isposdef(information) ||
        throw(ArgumentError("information matrix must be positive definite"))
    Float64(real(logdet(information)))
end

"""Geometry-independent metrics for a defect-minus-reference dictionary."""
function operator_metrics(
    dictionary::AbstractMatrix;
    noise_covariance::AbstractMatrix=Matrix{Float64}(
        I, size(dictionary, 1), size(dictionary, 1),
    ),
    signal_prior::AbstractMatrix=Matrix{Float64}(
        I, size(dictionary, 2), size(dictionary, 2),
    ),
    rank_tolerance::Union{Nothing,Real}=nothing,
)
    whitened = whiten_dictionary(dictionary; noise_covariance)
    singular_values = svdvals(whitened)
    largest = first(singular_values)
    largest > 0 || throw(ArgumentError("dictionary must contain a nonzero response"))
    tolerance = isnothing(rank_tolerance) ?
        max(size(whitened)...) * eps(Float64) * largest : Float64(rank_tolerance)
    tolerance >= 0 || throw(ArgumentError("rank tolerance must be nonnegative"))

    detectability = vec(sqrt.(sum(abs2, whitened; dims=1)))
    gram = normalized_gram(whitened)
    pair_coherence = size(gram, 1) == 1 ? 0.0 : maximum(
        gram[index, other]
        for index in axes(gram, 1), other in axes(gram, 2)
        if index != other
    )

    OperatorMetrics(
        size(whitened, 1),
        size(whitened, 2),
        Float64(minimum(detectability)),
        Float64(median(detectability)),
        Float64(pair_coherence),
        Float64(sum(abs2, singular_values) / largest^2),
        count(>(tolerance), singular_values),
        information_logdet(whitened, signal_prior),
    )
end

end
