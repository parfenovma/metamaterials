const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PILOT_ROOT = joinpath(PROJECT_ROOT, "tmp", "pilot_220khz")
const PILOT_FREQUENCY_HZ = 220.0e3
const SWEEP_FREQUENCIES_HZ = collect(160.0e3:20.0e3:280.0e3)
const PILOT_FINAL_TIME_S = 120.0e-6
const PILOT_AMPLITUDE_MM = 2.5

const PILOT_MESH_DIR = joinpath(PILOT_ROOT, "meshes")
const PILOT_MODEL_DIR = joinpath(PILOT_ROOT, "models")
const PILOT_SIGNAL_DIR = joinpath(PILOT_ROOT, "signals")
const PILOT_CHARACTERISTIC_DIR = joinpath(PILOT_ROOT, "characteristics")

const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

# Load stage dependencies before defining the functions below. Julia 1.12 uses
# strict world-age semantics for bindings created by a dynamic include.
if !isdefined(@__MODULE__, :MetamaterialProfiles)
    include(joinpath(@__DIR__, "profiles.jl"))
end
using .MetamaterialProfiles

if REQUESTED_STAGE == "mesh"
    include(joinpath(@__DIR__, "step1_mesher.jl"))
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "solve-one-")
    include(joinpath(@__DIR__, "step2_solver.jl"))
elseif REQUESTED_STAGE == "analyze"
    include(joinpath(@__DIR__, "step3_analyzer.jl"))
end

"""Small, deliberately coarse experiment used to verify the complete pipeline."""
function pilot_profiles()
    profiles = Main.MetamaterialProfiles
    profiles.WallProfile[
        profiles.SinusoidalProfile(0.0; periods=2),
        profiles.SinusoidalProfile(PILOT_AMPLITUDE_MM; periods=2),
        profiles.ExponentialProfile(PILOT_AMPLITUDE_MM; periods=2, sharpness=0.5),
        profiles.ExponentialProfile(PILOT_AMPLITUDE_MM; periods=1, sharpness=0.5),
        profiles.ExponentialProfile(PILOT_AMPLITUDE_MM; periods=2, sharpness=1.0),
        profiles.ExponentialProfile(PILOT_AMPLITUDE_MM; periods=2, sharpness=3.0),
        profiles.PowerProfile(PILOT_AMPLITUDE_MM; periods=2, power=2.0),
        profiles.PowerProfile(PILOT_AMPLITUDE_MM; periods=2, power=4.0),
        profiles.RoundedNotchProfile(
            PILOT_AMPLITUDE_MM;
            notch_count=4,
            notch_width_mm=1.4,
            end_margin_mm=1.2,
        ),
        profiles.LegacySinusoidalProfile(PILOT_AMPLITUDE_MM),
    ]
end

function pilot_signal_path(profile, frequency_hz::Real=PILOT_FREQUENCY_HZ)
    slug = Main.MetamaterialProfiles.profile_slug(profile)
    frequency_khz = frequency_hz / 1000.0
    joinpath(PILOT_SIGNAL_DIR, "data_$(slug)_F_$(frequency_khz).jld2")
end

function run_mesh_stage()
    mesh_config = MeshConfig(output_dir=PILOT_MESH_DIR)
    missing = filter(pilot_profiles()) do profile
        slug = Main.MetamaterialProfiles.profile_slug(profile)
        !isfile(joinpath(PILOT_MESH_DIR, "mesh_$(slug).msh"))
    end
    if isempty(missing)
        println("[=] All pilot meshes already exist")
    else
        generate_meshes(missing; mesh=mesh_config)
    end
end

function run_conversion_stage()
    mkpath(PILOT_MODEL_DIR)
    for profile in pilot_profiles()
        slug = Main.MetamaterialProfiles.profile_slug(profile)
        mesh_path = joinpath(PILOT_MESH_DIR, "mesh_$(slug).msh")
        model_path = joinpath(PILOT_MODEL_DIR, "model_$(slug).json")
        if isfile(model_path)
            println("  [=] Existing model: $model_path")
        else
            convert_mesh(mesh_path; output_dir=PILOT_MODEL_DIR)
        end
    end
end

function pilot_simulation(frequency_hz::Real=PILOT_FREQUENCY_HZ)
    simulation = SimulationConfig(
        frequencies_hz=[Float64(frequency_hz)],
        final_time_s=PILOT_FINAL_TIME_S,
        samples_per_period=30,
        pulse_cycles=4.0,
        save_vtk=false,
        skip_existing=true,
        model_dir=PILOT_MODEL_DIR,
        vtk_dir=joinpath(PILOT_ROOT, "vtk"),
        signal_dir=PILOT_SIGNAL_DIR,
    )
end

function run_single_solver_stage(profile_index::Integer, frequency_index::Integer)
    profiles = pilot_profiles()
    1 <= profile_index <= length(profiles) || error("invalid pilot profile index: $profile_index")
    1 <= frequency_index <= length(SWEEP_FREQUENCIES_HZ) ||
        error("invalid sweep frequency index: $frequency_index")
    frequency_hz = SWEEP_FREQUENCIES_HZ[frequency_index]
    run_acoustic_simulation(
        profiles[profile_index],
        frequency_hz;
        simulation=pilot_simulation(frequency_hz),
    )
end

function pilot_process_count()
    requested = parse(Int, get(ENV, "METAMATERIALS_PILOT_JOBS", "4"))
    requested > 0 || error("METAMATERIALS_PILOT_JOBS must be positive")
    min(requested, length(pilot_profiles()), Sys.CPU_THREADS)
end

function child_command(stage; threads::Integer=1)
    julia = Base.julia_cmd()
    script = @__FILE__
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=$threads $script --stage=$stage`
end

function run_solver_stage()
    profiles = pilot_profiles()
    mkpath(PILOT_SIGNAL_DIR)
    pending = [
        (profile_index, frequency_index)
        for frequency_index in eachindex(SWEEP_FREQUENCIES_HZ)
        for profile_index in eachindex(profiles)
        if !isfile(pilot_signal_path(
            profiles[profile_index],
            SWEEP_FREQUENCIES_HZ[frequency_index],
        ))
    ]
    if isempty(pending)
        println("[=] All pilot FEM results already exist")
        return
    end

    process_count = min(pilot_process_count(), length(pending))
    println("=== Pilot FEM: $(length(pending)) pending profiles, $process_count parallel processes ===")
    semaphore = Base.Semaphore(process_count)
    @sync for (profile_index, frequency_index) in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command(
                    "solve-one-$(profile_index)-$(frequency_index)";
                    threads=1,
                ))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function profile_row(
    profile,
    requested_frequency_hz,
    point,
    excess_envelope_delay_s,
)
    profiles = Main.MetamaterialProfiles
    shape_parameter_name = "none"
    shape_parameter_value = NaN
    periods = hasproperty(profile, :periods) ? getproperty(profile, :periods) : 2
    if profile isa profiles.ExponentialProfile
        shape_parameter_name = "sharpness"
        shape_parameter_value = profile.sharpness
    elseif profile isa profiles.PowerProfile
        shape_parameter_name = "power"
        shape_parameter_value = profile.power
    elseif profile isa profiles.LegacySinusoidalProfile
        shape_parameter_name = "correction_percent"
        shape_parameter_value = profile.correction_percent
    elseif profile isa profiles.RoundedNotchProfile
        shape_parameter_name = "notch_width_mm"
        shape_parameter_value = profile.notch_width_mm
        periods = profile.notch_count
    end

    (
        slug=profiles.profile_slug(profile),
        family=string(nameof(typeof(profile))),
        amplitude_mm=profiles.amplitude_mm(profile),
        periods=periods,
        shape_parameter_name,
        shape_parameter_value,
        minimum_gap_mm=profiles.minimum_gap_mm(profile),
        requested_frequency_hz,
        sampled_frequency_hz=point.frequency_hz,
        H_real=real(point.transfer),
        H_imag=imag(point.transfer),
        amplitude=point.amplitude,
        magnitude_squared=point.magnitude_squared,
        phase_rad=angle(point.transfer),
        group_delay_s=point.group_delay_s,
        excess_envelope_delay_s,
        valid=point.valid,
    )
end

"""Unwrap the independently sampled carrier phases for each geometry."""
function unwrap_frequency_rows(rows)
    result = copy(rows)
    for slug in unique(row.slug for row in rows)
        indices = sort(
            findall(row -> row.slug == slug, rows);
            by=index -> rows[index].requested_frequency_hz,
        )
        isempty(indices) && continue
        previous_raw = rows[first(indices)].phase_rad
        unwrapped = previous_raw
        result[first(indices)] = merge(rows[first(indices)], (phase_rad=unwrapped,))
        for index in Iterators.drop(indices, 1)
            raw = rows[index].phase_rad
            unwrapped += mod(raw - previous_raw + pi, 2.0 * pi) - pi
            result[index] = merge(rows[index], (phase_rad=unwrapped,))
            previous_raw = raw
        end
    end
    result
end

function write_library_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function run_analysis_stage()
    profiles = pilot_profiles()
    mkpath(PILOT_CHARACTERISTIC_DIR)

    rows = NamedTuple[]
    for frequency_hz in SWEEP_FREQUENCIES_HZ
        reference_path = pilot_signal_path(first(profiles), frequency_hz)
        isfile(reference_path) || error("reference result is missing: $reference_path")
        for profile in profiles
            sample_path = pilot_signal_path(profile, frequency_hz)
            isfile(sample_path) || error("FEM result is missing: $sample_path")
            output_path = save_analysis(
                sample_path;
                reference_path,
                output_dir=PILOT_CHARACTERISTIC_DIR,
            )
            saved = JLD2.load(output_path)
            point = Main.SpectralAnalysis.value_at_frequency(
                saved["calibrated_transmission"],
                frequency_hz,
            )
            push!(rows, profile_row(
                profile,
                frequency_hz,
                point,
                saved["excess_envelope_delay_s"],
            ))
        end
    end
    rows = unwrap_frequency_rows(rows)

    csv_path = joinpath(PILOT_ROOT, "library_frequency_sweep.csv")
    jld_path = joinpath(PILOT_ROOT, "library_frequency_sweep.jld2")
    write_library_csv(csv_path, rows)
    JLD2.jldsave(
        jld_path;
        format_version=1,
        requested_frequencies_hz=SWEEP_FREQUENCIES_HZ,
        final_time_s=PILOT_FINAL_TIME_S,
        reference_slug=Main.MetamaterialProfiles.profile_slug(first(profiles)),
        rows,
    )
    println("[+] Pilot library: $csv_path")
    println("[+] Machine-readable library: $jld_path")
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_conversion_stage()
    elseif stage == "solve"
        run_solver_stage()
    elseif startswith(stage, "solve-one-")
        parts = split(stage, '-')
        length(parts) == 4 || error("invalid single-solver stage: $stage")
        profile_index = parse(Int, parts[3])
        frequency_index = parse(Int, parts[4])
        run_single_solver_stage(profile_index, frequency_index)
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown pilot stage: $stage")
    end
end

function list_pilot()
    println("Pilot root: $PILOT_ROOT")
    println("frequencies_hz: $(join(SWEEP_FREQUENCIES_HZ, ", "))")
    for profile in pilot_profiles()
        slug = Main.MetamaterialProfiles.profile_slug(profile)
        gap = Main.MetamaterialProfiles.minimum_gap_mm(profile)
        println("  $slug, minimum_gap_mm=$gap")
    end
end

function run_pipeline()
    for stage in ("mesh", "convert", "solve", "analyze")
        println("\n=== Stage: $stage ===")
        run(child_command(stage; threads=1))
    end
end

function main(args=ARGS)
    if isempty(args)
        run_pipeline()
    elseif args == ["--list"]
        list_pilot()
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_pilot_sweep.jl [--list | --stage=mesh|convert|solve|analyze]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
