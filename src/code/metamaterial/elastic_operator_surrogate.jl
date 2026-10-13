module ElasticOperatorSurrogate

using LinearAlgebra

export SurrogateConfig,
       DefectHypothesis,
       validate,
       element_positions,
       localization_hypotheses,
       mirror_hypothesis_index,
       elastic_green_tensor,
       response_matrix,
       response_library,
       uniform_port,
       fixed_focus_port,
       dispersive_port,
       cosine_ports,
       one_port_dictionary,
       multiport_dictionary,
       fmc_dictionary,
       ideal_compression_dictionary,
       maximum_reciprocity_error,
       maximum_mirror_error

Base.@kwdef struct SurrogateConfig
    pressure_speed::Float64 = 1.0
    shear_speed::Float64 = 0.5
    reference_frequency::Float64 = 1.0
    frequency_ratios::Vector{Float64} = collect(range(0.8, 1.2; length=7))
    aperture_wavelengths::Float64 = 8.0
    element_count::Int = 17
    pressure_attenuation_per_wavelength::Float64 = 0.01
    shear_attenuation_per_wavelength::Float64 = 0.025
    shear_amplitude::Float64 = 0.7
end

struct DefectHypothesis
    id::Symbol
    x_wavelengths::Float64
    y_wavelengths::Float64
    scattering_tensor::Matrix{Float64}
end

function validate(config::SurrogateConfig)
    config.pressure_speed > 0 || throw(ArgumentError("P speed must be positive"))
    0 < config.shear_speed < config.pressure_speed ||
        throw(ArgumentError("S speed must lie between zero and P speed"))
    config.reference_frequency > 0 || throw(ArgumentError("reference frequency must be positive"))
    !isempty(config.frequency_ratios) || throw(ArgumentError("frequency grid must not be empty"))
    all(>(0), config.frequency_ratios) || throw(ArgumentError("frequencies must be positive"))
    config.aperture_wavelengths > 0 || throw(ArgumentError("aperture must be positive"))
    config.element_count >= 3 || throw(ArgumentError("at least three virtual elements are required"))
    isodd(config.element_count) || throw(ArgumentError("element count must be odd for exact C_s symmetry"))
    config.pressure_attenuation_per_wavelength >= 0 ||
        throw(ArgumentError("P attenuation must be nonnegative"))
    config.shear_attenuation_per_wavelength >= 0 ||
        throw(ArgumentError("S attenuation must be nonnegative"))
    config.shear_amplitude >= 0 || throw(ArgumentError("S amplitude must be nonnegative"))
    config
end

function validate(hypothesis::DefectHypothesis)
    hypothesis.x_wavelengths > 0 || throw(ArgumentError("defect must lie beyond the aperture"))
    size(hypothesis.scattering_tensor) == (2, 2) ||
        throw(DimensionMismatch("scattering tensor must be 2 x 2"))
    isapprox(hypothesis.scattering_tensor, hypothesis.scattering_tensor'; atol=1.0e-12) ||
        throw(ArgumentError("scattering tensor must be symmetric"))
    hypothesis
end

function element_positions(config::SurrogateConfig=SurrogateConfig())
    validate(config)
    collect(range(
        -config.aperture_wavelengths / 2,
        config.aperture_wavelengths / 2;
        length=config.element_count,
    ))
end

function localization_hypotheses(
    ;
    depths=collect(range(2.5, 6.0; length=8)),
    lateral_positions=collect(range(-2.0, 2.0; length=17)),
)
    hypotheses = DefectHypothesis[]
    for (depth_index, depth) in enumerate(depths)
        for (lateral_index, lateral) in enumerate(lateral_positions)
            suffix = "x$(depth_index)_y$(lateral_index)"
            push!(hypotheses, DefectHypothesis(
                Symbol("iso_$(suffix)"),
                Float64(depth),
                Float64(lateral),
                Matrix{Float64}(I, 2, 2),
            ))
        end
    end
    validate.(hypotheses)
end

function mirror_hypothesis_index(hypotheses::AbstractVector{DefectHypothesis})
    indices = Int[]
    for hypothesis in hypotheses
        matches = findall(other ->
            isapprox(other.x_wavelengths, hypothesis.x_wavelengths; atol=1.0e-12) &&
            isapprox(other.y_wavelengths, -hypothesis.y_wavelengths; atol=1.0e-12) &&
            isapprox(other.scattering_tensor, [1.0 0.0; 0.0 -1.0] *
                hypothesis.scattering_tensor * [1.0 0.0; 0.0 -1.0]; atol=1.0e-12),
            hypotheses,
        )
        length(matches) == 1 || error("hypothesis set is not closed under reflection")
        push!(indices, only(matches))
    end
    indices
end

function modal_green_components(
    source_y::Real,
    target_x::Real,
    target_y::Real,
    frequency_ratio::Real,
    config::SurrogateConfig,
)
    dx = Float64(target_x)
    dy = Float64(target_y - source_y)
    distance = hypot(dx, dy)
    distance > 0 || throw(ArgumentError("source and target must not coincide"))
    ep = [dx, dy] / distance
    es = [-ep[2], ep[1]]
    frequency = config.reference_frequency * Float64(frequency_ratio)
    kp = 2pi * frequency / config.pressure_speed
    ks = 2pi * frequency / config.shear_speed
    pressure_phase = exp(
        im * kp * distance -
        config.pressure_attenuation_per_wavelength * Float64(frequency_ratio) * distance,
    ) / sqrt(distance)
    shear_phase = config.shear_amplitude * exp(
        im * ks * distance -
        config.shear_attenuation_per_wavelength * Float64(frequency_ratio) * distance,
    ) / sqrt(distance)
    pressure_phase .* (ep * ep'), shear_phase .* (es * es')
end

function elastic_green_tensor(
    source_y::Real,
    target_x::Real,
    target_y::Real,
    frequency_ratio::Real,
    config::SurrogateConfig=SurrogateConfig(),
)
    validate(config)
    gp, gs = modal_green_components(source_y, target_x, target_y, frequency_ratio, config)
    gp + gs
end

function response_matrix(
    hypothesis::DefectHypothesis,
    frequency_ratio::Real,
    config::SurrogateConfig=SurrogateConfig(),
)
    validate(config)
    validate(hypothesis)
    positions = element_positions(config)
    element_count = length(positions)
    source_polarization = ComplexF64[1, 0]
    green_vectors = Matrix{ComplexF64}(undef, 2, element_count)
    for (index, position) in enumerate(positions)
        green_vectors[:, index] = elastic_green_tensor(
            position,
            hypothesis.x_wavelengths,
            hypothesis.y_wavelengths,
            frequency_ratio,
            config,
        ) * source_polarization
    end
    transpose(green_vectors) * hypothesis.scattering_tensor * green_vectors
end

function response_library(
    hypotheses::AbstractVector{DefectHypothesis},
    config::SurrogateConfig=SurrogateConfig(),
)
    validate(config)
    validate.(hypotheses)
    matrices = Array{ComplexF64}(undef,
        config.element_count,
        config.element_count,
        length(config.frequency_ratios),
        length(hypotheses),
    )
    for (hypothesis_index, hypothesis) in enumerate(hypotheses)
        for (frequency_index, frequency_ratio) in enumerate(config.frequency_ratios)
            matrices[:, :, frequency_index, hypothesis_index] = response_matrix(
                hypothesis, frequency_ratio, config,
            )
        end
    end
    matrices
end

uniform_port(config::SurrogateConfig=SurrogateConfig()) =
    fill(inv(sqrt(config.element_count)), config.element_count)

function fixed_focus_port(
    frequency_ratio::Real,
    focus_x::Real,
    focus_y::Real,
    config::SurrogateConfig=SurrogateConfig(),
)
    source_polarization = ComplexF64[1, 0]
    steering = ComplexF64[
        dot(source_polarization, elastic_green_tensor(
            position, focus_x, focus_y, frequency_ratio, config,
        ) * source_polarization)
        for position in element_positions(config)
    ]
    phases = map(value -> abs(value) > 0 ? conj(value) / abs(value) : one(value), steering)
    phases / norm(phases)
end

function dispersive_port(
    frequency_ratio::Real,
    config::SurrogateConfig=SurrogateConfig();
    symmetry::Symbol=:even,
)
    positions = element_positions(config)
    normalized = positions ./ maximum(abs, positions)
    delays = if symmetry == :even
        0.15 .+ 0.75 .* abs.(normalized) .^ 1.7 .+
        0.10 .* cos.(2pi .* normalized)
    elseif symmetry == :asymmetric
        raw = @. 0.52 + 0.23 * sin(2pi * (0.37 * normalized + 0.11)) +
            0.25 * normalized + 0.08 * sin(5pi * normalized)
        raw .- minimum(raw)
    else
        throw(ArgumentError("symmetry must be :even or :asymmetric"))
    end
    weights = exp.(-2pi * im * Float64(frequency_ratio) .* delays)
    weights / norm(weights)
end

function cosine_ports(port_count::Integer, config::SurrogateConfig=SurrogateConfig())
    1 <= port_count <= config.element_count ||
        throw(ArgumentError("port count must lie between one and element count"))
    count = config.element_count
    ports = Matrix{Float64}(undef, count, port_count)
    sample_indices = collect(0:(count - 1))
    for mode in 0:(port_count - 1)
        vector = cos.(pi .* (sample_indices .+ 0.5) .* mode ./ count)
        ports[:, mode + 1] = vector / norm(vector)
    end
    isapprox(ports' * ports, Matrix{Float64}(I, port_count, port_count); atol=1.0e-12) ||
        error("cosine port basis lost orthonormality")
    ports
end

function one_port_dictionary(
    library::AbstractArray{<:Complex,4},
    port_function::Function,
    config::SurrogateConfig=SurrogateConfig(),
)
    hypothesis_count = size(library, 4)
    dictionary = Matrix{ComplexF64}(undef, length(config.frequency_ratios), hypothesis_count)
    for (frequency_index, frequency_ratio) in enumerate(config.frequency_ratios)
        port = port_function(frequency_ratio)
        length(port) == config.element_count || throw(DimensionMismatch("one weight per element is required"))
        isapprox(sum(abs2, port), 1.0; atol=1.0e-12) || error("port weights must have unit energy")
        for hypothesis_index in 1:hypothesis_count
            response = view(library, :, :, frequency_index, hypothesis_index)
            dictionary[frequency_index, hypothesis_index] = transpose(port) * response * port
        end
    end
    dictionary
end

function multiport_dictionary(
    library::AbstractArray{<:Complex,4},
    ports::AbstractMatrix,
    config::SurrogateConfig=SurrogateConfig(),
)
    size(ports, 1) == config.element_count || throw(DimensionMismatch("port rows must match elements"))
    port_count = size(ports, 2)
    isapprox(ports' * ports, Matrix{eltype(ports)}(I, port_count, port_count); atol=1.0e-10) ||
        error("port columns must be orthonormal")
    hypothesis_count = size(library, 4)
    dictionary = Matrix{ComplexF64}(undef,
        length(config.frequency_ratios) * port_count^2,
        hypothesis_count,
    )
    energy_scale = inv(sqrt(port_count))
    for hypothesis_index in 1:hypothesis_count
        row = 1
        for frequency_index in eachindex(config.frequency_ratios)
            response = view(library, :, :, frequency_index, hypothesis_index)
            compressed = transpose(ports) * response * ports .* energy_scale
            for transmit_port in 1:port_count, receive_port in 1:port_count
                dictionary[row, hypothesis_index] = compressed[receive_port, transmit_port]
                row += 1
            end
        end
    end
    dictionary
end

function fmc_dictionary(
    library::AbstractArray{<:Complex,4},
    config::SurrogateConfig=SurrogateConfig(),
)
    hypothesis_count = size(library, 4)
    rows_per_frequency = config.element_count^2
    dictionary = Matrix{ComplexF64}(undef,
        length(config.frequency_ratios) * rows_per_frequency,
        hypothesis_count,
    )
    energy_scale = inv(sqrt(config.element_count))
    for hypothesis_index in 1:hypothesis_count
        row = 1
        for frequency_index in eachindex(config.frequency_ratios)
            response = view(library, :, :, frequency_index, hypothesis_index)
            for transmit_element in 1:config.element_count, receive_element in 1:config.element_count
                dictionary[row, hypothesis_index] =
                    energy_scale * response[receive_element, transmit_element]
                row += 1
            end
        end
    end
    dictionary
end

function ideal_compression_dictionary(dictionary::AbstractMatrix, measurement_count::Integer)
    1 <= measurement_count <= min(size(dictionary)...) ||
        throw(ArgumentError("ideal measurement count exceeds dictionary rank dimensions"))
    decomposition = svd(dictionary; full=false)
    compressed = Diagonal(decomposition.S[1:measurement_count]) *
        decomposition.Vt[1:measurement_count, :]
    retained_energy = sum(abs2, decomposition.S[1:measurement_count]) /
        sum(abs2, decomposition.S)
    compressed, retained_energy, decomposition.S
end

function maximum_reciprocity_error(library::AbstractArray{<:Complex,4})
    maximum(
        norm(matrix - transpose(matrix)) / max(norm(matrix), eps())
        for matrix in eachslice(library; dims=(3, 4))
    )
end

function maximum_mirror_error(
    library::AbstractArray{<:Complex,4},
    hypotheses::AbstractVector{DefectHypothesis},
)
    mirrors = mirror_hypothesis_index(hypotheses)
    reversal = reverse(1:size(library, 1))
    maximum(
        norm(
            view(library, :, :, frequency_index, mirrors[hypothesis_index]) -
            view(library, reversal, reversal, frequency_index, hypothesis_index),
        ) / max(norm(view(library, :, :, frequency_index, hypothesis_index)), eps())
        for frequency_index in axes(library, 3), hypothesis_index in axes(library, 4)
    )
end

end
