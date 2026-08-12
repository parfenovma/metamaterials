module ElasticPortModes

using LinearAlgebra

export PortModeConfig,
       PortMode,
       solve_port_modes,
       propagating_modes,
       evanescent_modes,
       fundamental_quasi_longitudinal,
       mode_overlap,
       track_fundamental_branch

"""Material, frequency and 1D cross-section discretization for a straight strip port."""
Base.@kwdef struct PortModeConfig
    height_m::Float64 = 4.2e-3
    density::Float64 = 1210.0
    pressure_wave_speed::Float64 = 2340.0
    shear_wave_speed::Float64 = 1170.0
    frequency_hz::Float64 = 242.0e3
    element_count::Int = 80
    propagation_tolerance::Float64 = 1.0e-7
    residual_tolerance::Float64 = 1.0e-7
end

"""
One guided mode with spatial dependence `u(y) exp(gamma*x) exp(i*omega*t)`.

For a right-going propagating mode `gamma = -i*k`. Propagating modes are
normalized to unit absolute time-averaged power per metre out of plane.
Evanescent modes are normalized to unit cross-sectional displacement norm.
"""
struct PortMode
    gamma_per_m::ComplexF64
    wavenumber_per_m::ComplexF64
    kind::Symbol
    direction::Symbol
    power_w_per_m::Float64
    relative_residual::Float64
    parity_score::Float64
    parity::Symbol
    p_fraction::Float64
    s_fraction::Float64
    axial_displacement_fraction::Float64
    y_m::Vector{Float64}
    displacement_x::Vector{ComplexF64}
    displacement_y::Vector{ComplexF64}
    traction_xx_pa::Vector{ComplexF64}
    traction_xy_pa::Vector{ComplexF64}
end

function validate_config(config::PortModeConfig)
    config.height_m > 0 || throw(ArgumentError("port height must be positive"))
    config.density > 0 || throw(ArgumentError("density must be positive"))
    config.pressure_wave_speed > 0 || throw(ArgumentError("P-wave speed must be positive"))
    config.shear_wave_speed > 0 || throw(ArgumentError("S-wave speed must be positive"))
    config.pressure_wave_speed > sqrt(2) * config.shear_wave_speed ||
        throw(ArgumentError("material parameters give a non-positive first Lame constant"))
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.element_count >= 4 || throw(ArgumentError("at least four cross-section elements are required"))
    config.propagation_tolerance > 0 || throw(ArgumentError("propagation tolerance must be positive"))
    config.residual_tolerance > 0 || throw(ArgumentError("residual tolerance must be positive"))
    nothing
end

function material_constants(config)
    mu = config.density * config.shear_wave_speed^2
    lambda = config.density * config.pressure_wave_speed^2 - 2mu
    lambda, mu
end

"""
Assemble `Q(gamma) = A0 + gamma*A1 + gamma^2*A2` for a free strip.

Linear Lagrange elements are used only across the port. The free conditions
`sigma_xy = sigma_yy = 0` at the upper and lower faces are natural boundary
conditions of this weak form.
"""
function assemble_quadratic_pencil(config::PortModeConfig)
    validate_config(config)
    lambda, mu = material_constants(config)
    rho = config.density
    omega = 2pi * config.frequency_hz
    node_count = config.element_count + 1
    dof_count = 2node_count
    y = collect(range(-config.height_m / 2, config.height_m / 2; length=node_count))
    a0 = zeros(ComplexF64, dof_count, dof_count)
    a1 = zeros(ComplexF64, dof_count, dof_count)
    a2 = zeros(ComplexF64, dof_count, dof_count)
    constitutive = ComplexF64[
        lambda + 2mu lambda 0;
        lambda lambda + 2mu 0;
        0 0 mu
    ]
    gauss_points = (-inv(sqrt(3.0)), inv(sqrt(3.0)))

    for element in 1:config.element_count
        h = y[element + 1] - y[element]
        jacobian = h / 2
        element_dofs = [2element - 1, 2element, 2element + 1, 2element + 2]
        local_a0 = zeros(ComplexF64, 4, 4)
        local_a1 = zeros(ComplexF64, 4, 4)
        local_a2 = zeros(ComplexF64, 4, 4)
        for xi in gauss_points
            shape = ((1 - xi) / 2, (1 + xi) / 2)
            derivative = (-1 / h, 1 / h)
            b0 = ComplexF64[
                0 0 0 0;
                0 derivative[1] 0 derivative[2];
                derivative[1] 0 derivative[2] 0
            ]
            b1 = ComplexF64[
                shape[1] 0 shape[2] 0;
                0 0 0 0;
                0 shape[1] 0 shape[2]
            ]
            shape_matrix = ComplexF64[
                shape[1] 0 shape[2] 0;
                0 shape[1] 0 shape[2]
            ]
            mass = rho * (shape_matrix' * shape_matrix)
            k00 = b0' * constitutive * b0
            k01 = b0' * constitutive * b1
            k10 = b1' * constitutive * b0
            k11 = b1' * constitutive * b1
            local_a0 .+= (omega^2 * mass - k00) * jacobian
            local_a1 .+= (k10 - k01) * jacobian
            local_a2 .+= k11 * jacobian
        end
        a0[element_dofs, element_dofs] .+= local_a0
        a1[element_dofs, element_dofs] .+= local_a1
        a2[element_dofs, element_dofs] .+= local_a2
    end
    (; y, a0, a1, a2)
end

function linearized_eigenpairs(a0, a1, a2)
    dof_count = size(a0, 1)
    identity = Matrix{ComplexF64}(I, dof_count, dof_count)
    zero_block = zeros(ComplexF64, dof_count, dof_count)
    # A2 is positive definite for an elastic strip. Reducing to a standard
    # companion problem is substantially better conditioned than passing the
    # SI-scaled generalized pencil directly to QZ.
    companion = [zero_block identity; -(a2 \ a0) -(a2 \ a1)]
    decomposition = eigen(companion)
    decomposition.values, decomposition.vectors[1:dof_count, :]
end

function cross_section_norm(y, ux, uy)
    result = 0.0
    for element in 1:(length(y) - 1)
        h = y[element + 1] - y[element]
        result += h / 2 * (
            abs2(ux[element]) + abs2(uy[element]) +
            abs2(ux[element + 1]) + abs2(uy[element + 1])
        )
    end
    result
end

function component_norm(y, values)
    result = 0.0
    for element in 1:(length(y) - 1)
        h = y[element + 1] - y[element]
        result += h / 2 * (abs2(values[element]) + abs2(values[element + 1]))
    end
    result
end

function modal_diagnostics(y, ux, uy, gamma, config)
    lambda, mu = material_constants(config)
    omega = 2pi * config.frequency_hz
    power = 0.0
    p_measure = 0.0
    s_measure = 0.0
    gauss_points = (-inv(sqrt(3.0)), inv(sqrt(3.0)))
    for element in 1:(length(y) - 1)
        h = y[element + 1] - y[element]
        jacobian = h / 2
        derivative_x = (ux[element + 1] - ux[element]) / h
        derivative_y = (uy[element + 1] - uy[element]) / h
        for xi in gauss_points
            n1, n2 = (1 - xi) / 2, (1 + xi) / 2
            value_x = n1 * ux[element] + n2 * ux[element + 1]
            value_y = n1 * uy[element] + n2 * uy[element + 1]
            sigma_xx = (lambda + 2mu) * gamma * value_x + lambda * derivative_y
            sigma_xy = mu * (derivative_x + gamma * value_y)
            velocity_x = im * omega * value_x
            velocity_y = im * omega * value_y
            power += -0.5 * real(
                sigma_xx * conj(velocity_x) + sigma_xy * conj(velocity_y),
            ) * jacobian
            divergence = gamma * value_x + derivative_y
            rotation = gamma * value_y - derivative_x
            p_measure += abs2(divergence) * jacobian
            s_measure += abs2(rotation) * jacobian
        end
    end
    total = p_measure + s_measure
    p_fraction = total > 0 ? p_measure / total : NaN
    power, p_fraction, 1 - p_fraction
end

function nodal_tractions(y, ux, uy, gamma, config)
    lambda, mu = material_constants(config)
    count = length(y)
    derivative_x = zeros(ComplexF64, count)
    derivative_y = zeros(ComplexF64, count)
    weights = zeros(Int, count)
    for element in 1:(count - 1)
        h = y[element + 1] - y[element]
        dx = (ux[element + 1] - ux[element]) / h
        dy = (uy[element + 1] - uy[element]) / h
        for node in (element, element + 1)
            derivative_x[node] += dx
            derivative_y[node] += dy
            weights[node] += 1
        end
    end
    derivative_x ./= weights
    derivative_y ./= weights
    traction_xx = (lambda + 2mu) .* gamma .* ux .+ lambda .* derivative_y
    traction_xy = mu .* (derivative_x .+ gamma .* uy)
    traction_xx, traction_xy
end

function vector_parity(ux, uy)
    mirrored_x = reverse(ux)
    mirrored_y = -reverse(uy)
    denominator = sum(abs2, ux) + sum(abs2, uy)
    denominator > 0 || return NaN
    numerator = real(sum(conj.(ux) .* mirrored_x) + sum(conj.(uy) .* mirrored_y))
    clamp(numerator / denominator, -1.0, 1.0)
end

function classify_mode(gamma, power, config)
    scale = max(abs(imag(gamma)), 2pi * config.frequency_hz / config.shear_wave_speed, 1.0)
    if abs(real(gamma)) <= config.propagation_tolerance * scale
        direction = power > 0 ? :right : power < 0 ? :left : :neutral
        return :propagating, direction
    end
    :evanescent, real(gamma) < 0 ? :right : :left
end

function make_mode(gamma, vector, pencil, config)
    node_count = length(pencil.y)
    ux = ComplexF64[vector[2index - 1] for index in 1:node_count]
    uy = ComplexF64[vector[2index] for index in 1:node_count]
    displacement_norm = sqrt(cross_section_norm(pencil.y, ux, uy))
    displacement_norm > 0 || return nothing
    ux ./= displacement_norm
    uy ./= displacement_norm
    power, p_fraction, s_fraction = modal_diagnostics(pencil.y, ux, uy, gamma, config)
    kind, direction = classify_mode(gamma, power, config)
    if kind == :propagating && abs(power) > 100eps(Float64)
        normalization = sqrt(abs(power))
        ux ./= normalization
        uy ./= normalization
        power, p_fraction, s_fraction = modal_diagnostics(pencil.y, ux, uy, gamma, config)
    end
    traction_xx, traction_xy = nodal_tractions(pencil.y, ux, uy, gamma, config)
    qphi = pencil.a0 * vector + gamma * (pencil.a1 * vector) + gamma^2 * (pencil.a2 * vector)
    denominator = (
        norm(pencil.a0 * vector) +
        abs(gamma) * norm(pencil.a1 * vector) +
        abs2(gamma) * norm(pencil.a2 * vector)
    )
    residual = denominator > 0 ? norm(qphi) / denominator : Inf
    parity_score = vector_parity(ux, uy)
    parity = parity_score >= 0.8 ? :symmetric : parity_score <= -0.8 ? :antisymmetric : :mixed
    axial_norm = component_norm(pencil.y, ux)
    transverse_norm = component_norm(pencil.y, uy)
    axial_fraction = axial_norm / (axial_norm + transverse_norm)
    PortMode(
        ComplexF64(gamma),
        ComplexF64(im * gamma),
        kind,
        direction,
        Float64(power),
        Float64(residual),
        Float64(parity_score),
        parity,
        Float64(p_fraction),
        Float64(s_fraction),
        Float64(axial_fraction),
        copy(pencil.y),
        ux,
        uy,
        traction_xx,
        traction_xy,
    )
end

"""Solve all finite discrete port modes at the configured real frequency."""
function solve_port_modes(config::PortModeConfig=PortModeConfig())
    pencil = assemble_quadratic_pencil(config)
    # Solve for the dimensionless exponent zeta = gamma*height. Without this
    # scaling the companion pencil mixes entries separated by many powers of
    # the SI length unit and loses several digits in the propagating roots.
    values, vectors = linearized_eigenpairs(
        pencil.a0,
        pencil.a1 / config.height_m,
        pencil.a2 / config.height_m^2,
    )
    modes = PortMode[]
    for index in eachindex(values)
        gamma = values[index] / config.height_m
        isfinite(real(gamma)) && isfinite(imag(gamma)) || continue
        mode = make_mode(gamma, vectors[:, index], pencil, config)
        mode === nothing && continue
        mode.relative_residual <= config.residual_tolerance || continue
        push!(modes, mode)
    end
    sort!(modes; by=mode -> (
        mode.kind == :propagating ? 0 : 1,
        mode.direction == :right ? 0 : 1,
        abs(real(mode.gamma_per_m)),
        abs(imag(mode.gamma_per_m)),
    ))
end

propagating_modes(modes; direction=nothing) = filter(
    mode -> mode.kind == :propagating && (direction === nothing || mode.direction == direction),
    modes,
)

evanescent_modes(modes; direction=nothing) = filter(
    mode -> mode.kind == :evanescent && (direction === nothing || mode.direction == direction),
    modes,
)

"""
Return the right-going fundamental symmetric branch.

At low frequency this is the extensional S0-like mode. At frequencies near a
branch interaction its instantaneous divergence/curl content can change, so
symmetry is a safer single-frequency identifier than simply maximizing the P
diagnostic. Frequency continuation should be used when several symmetric
propagating branches coexist.
"""
function fundamental_quasi_longitudinal(modes)
    candidates = propagating_modes(modes; direction=:right)
    isempty(candidates) && error("no right-going propagating mode was found")
    symmetric = filter(mode -> mode.parity == :symmetric, candidates)
    pool = isempty(symmetric) ? candidates : symmetric
    first(sort(pool; by=mode -> (-mode.p_fraction, abs(mode.wavenumber_per_m))))
end

function cross_section_inner_product(first_mode, second_mode)
    first_mode.y_m == second_mode.y_m ||
        throw(ArgumentError("mode overlap requires the same cross-section mesh"))
    y = first_mode.y_m
    result = 0.0 + 0im
    for element in 1:(length(y) - 1)
        h = y[element + 1] - y[element]
        first_left = (
            conj(first_mode.displacement_x[element]) * second_mode.displacement_x[element] +
            conj(first_mode.displacement_y[element]) * second_mode.displacement_y[element]
        )
        first_right = (
            conj(first_mode.displacement_x[element + 1]) * second_mode.displacement_x[element + 1] +
            conj(first_mode.displacement_y[element + 1]) * second_mode.displacement_y[element + 1]
        )
        result += h / 2 * (first_left + first_right)
    end
    result
end

"""Normalized displacement overlap used for frequency continuation."""
function mode_overlap(first_mode::PortMode, second_mode::PortMode)
    numerator = abs(cross_section_inner_product(first_mode, second_mode))
    first_norm = real(cross_section_inner_product(first_mode, first_mode))
    second_norm = real(cross_section_inner_product(second_mode, second_mode))
    denominator = sqrt(max(first_norm * second_norm, 0.0))
    denominator > 0 ? clamp(numerator / denominator, 0.0, 1.0) : NaN
end

"""
Track the right-going fundamental symmetric branch from low to high frequency.

The first frequency should lie in the low-frequency regime where the S0-like
mode is unambiguous. Subsequent modes are selected by maximum normalized
cross-sectional overlap within the same vector-parity sector.
"""
function track_fundamental_branch(
    config::PortModeConfig,
    frequencies_hz::AbstractVector{<:Real},
)
    issorted(frequencies_hz) || throw(ArgumentError("continuation frequencies must be sorted"))
    isempty(frequencies_hz) && return PortMode[]
    branch = PortMode[]
    overlaps = Float64[]
    for frequency_hz in frequencies_hz
        frequency_config = PortModeConfig(
            height_m=config.height_m,
            density=config.density,
            pressure_wave_speed=config.pressure_wave_speed,
            shear_wave_speed=config.shear_wave_speed,
            frequency_hz=Float64(frequency_hz),
            element_count=config.element_count,
            propagation_tolerance=config.propagation_tolerance,
            residual_tolerance=config.residual_tolerance,
        )
        modes = solve_port_modes(frequency_config)
        candidates = filter(
            mode -> mode.kind == :propagating &&
                    mode.direction == :right &&
                    mode.parity == :symmetric,
            modes,
        )
        isempty(candidates) && error("no right-going symmetric mode at $frequency_hz Hz")
        if isempty(branch)
            selected = fundamental_quasi_longitudinal(modes)
            push!(overlaps, 1.0)
        else
            candidate_overlaps = mode_overlap.(Ref(branch[end]), candidates)
            selected_index = argmax(candidate_overlaps)
            selected = candidates[selected_index]
            push!(overlaps, candidate_overlaps[selected_index])
        end
        push!(branch, selected)
    end
    (; modes=branch, overlaps)
end

end
