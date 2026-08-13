module DesignAluminiumHornTTDPassiveFeed

using JLD2

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = joinpath(
    PROJECT_ROOT, "tmp", "aluminium_horn_ttd_full_diffuser_aperture_242khz",
)
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_feed")
const WEIGHT_PATH = let provisional = joinpath(
        DESIGN_ROOT, "15_selected_pulse_drive.jld2",
    )
    isfile(provisional) ? provisional : joinpath(
        PROJECT_ROOT, "results", "aluminium_horn_ttd_p6_242khz",
        "fifteen_channel_diffuser", "15_selected_pulse_drive.jld2",
    )
end
const BRANCH_WIDTH_MM = 7.0
const LOWER_FREQUENCY_HZ = 226.04e3
const ALUMINIUM_CP_M_S = 6122.102437409232
const ADIABATIC_LIMIT = 0.05
const VERIFIED_SYSTEM_GAIN = 2.1260609847797616
const TARGET_GAIN = 2.0

mutable struct FeedNode
    id::Int
    power::Float64
    leaves::Vector{Int}
    left::Union{Nothing, FeedNode}
    right::Union{Nothing, FeedNode}
end

FeedNode(id, power, leaf) = FeedNode(id, power, [leaf], nothing, nothing)
isleaf(node) = isnothing(node.left)

function selected_full_weights()
    half = Float64.(JLD2.load(WEIGHT_PATH)["pressure_weights"])
    vcat(reverse(half[2:end]), half)
end

function balanced_tree(power)
    count = length(power)
    full_mask = (1 << count) - 1
    mask_power = zeros(Float64, full_mask + 1)
    for mask in 1:full_mask
        bit = trailing_zeros(mask)
        mask_power[mask + 1] = mask_power[(mask & (mask - 1)) + 1] + power[bit + 1]
    end
    memo = Dict{Int, Tuple{Float64, Int}}()
    function solve(mask)
        haskey(memo, mask) && return memo[mask]
        leaf_count = count_ones(mask)
        leaf_count == 1 && return (memo[mask] = (1.0, 0))
        left_count = leaf_count ÷ 2
        first_bit = trailing_zeros(mask)
        best_score = -Inf
        best_split = 0
        submask = mask
        while submask > 0
            if count_ones(submask) == left_count &&
               (isodd(leaf_count) || !iszero(submask & (1 << first_bit)))
                other = xor(mask, submask)
                left_score = first(solve(submask))
                right_score = first(solve(other))
                fraction = mask_power[submask + 1] / mask_power[mask + 1]
                score = min(fraction, 1 - fraction, left_score, right_score)
                if score > best_score + 1e-12 ||
                   (isapprox(score, best_score; atol=1e-12) && submask < best_split)
                    best_score = score
                    best_split = submask
                end
            end
            submask = (submask - 1) & mask
        end
        best_split > 0 || error("balanced feed partition was not found")
        memo[mask] = (best_score, best_split)
    end
    solve(full_mask)
    next_id = Ref(count + 1)
    internal = FeedNode[]
    function build(mask)
        if count_ones(mask) == 1
            leaf = trailing_zeros(mask) + 1
            return FeedNode(leaf, power[leaf], leaf)
        end
        split = last(solve(mask))
        left = build(split)
        right = build(xor(mask, split))
        parent = FeedNode(
            next_id[], mask_power[mask + 1], sort(vcat(left.leaves, right.leaves)),
            left, right,
        )
        next_id[] += 1
        push!(internal, parent)
        parent
    end
    build(full_mask), internal
end

"""
    efficiency_optimal_tree(power, eta_equal, eta_other; contiguous=false)

Find the binary tree that maximizes the retained useful power in the
optimistic cascade model. Every node divides power according to the target
leaf powers. A 50/50 node retains `eta_equal`; every other node retains
`eta_other`.

With `contiguous=true`, every subtree is restricted to a contiguous interval
of aperture channels. The unrestricted result is only a mathematical ceiling:
it may require branch crossings and neither result proves that all requested
split ratios share the two measured efficiencies.
"""
function efficiency_optimal_tree(power, eta_equal, eta_other; contiguous=false)
    powers = Float64.(power)
    !isempty(powers) || throw(ArgumentError("at least one leaf is required"))
    all(>(0), powers) || throw(ArgumentError("leaf powers must be positive"))
    0 <= eta_equal <= 1 || throw(ArgumentError("eta_equal must be in [0, 1]"))
    0 <= eta_other <= 1 || throw(ArgumentError("eta_other must be in [0, 1]"))
    count = length(powers)
    count <= 20 || throw(ArgumentError(
        "unrestricted subset search is limited to 20 leaves",
    ))

    split_efficiency(fraction) = isapprox(fraction, 0.5; atol=1e-10) ?
                                 Float64(eta_equal) : Float64(eta_other)
    next_id = Ref(count + 1)
    internal = FeedNode[]

    if contiguous
        interval_power(first_index, last_index) =
            sum(@view powers[first_index:last_index])
        memo = Dict{Tuple{Int, Int}, Tuple{Float64, Int}}()
        function solve_interval(first_index, last_index)
            first_index == last_index && return (1.0, 0)
            key = (first_index, last_index)
            haskey(memo, key) && return memo[key]
            total = interval_power(first_index, last_index)
            best_score = -Inf
            best_split = 0
            for split in first_index:(last_index - 1)
                left_power = interval_power(first_index, split)
                fraction = left_power / total
                score = split_efficiency(fraction) * (
                    fraction * first(solve_interval(first_index, split)) +
                    (1 - fraction) *
                    first(solve_interval(split + 1, last_index))
                )
                if score > best_score + 1e-14
                    best_score = score
                    best_split = split
                end
            end
            memo[key] = (best_score, best_split)
        end
        function build_interval(first_index, last_index)
            first_index == last_index &&
                return FeedNode(first_index, powers[first_index], first_index)
            split = last(solve_interval(first_index, last_index))
            left = build_interval(first_index, split)
            right = build_interval(split + 1, last_index)
            parent = FeedNode(
                next_id[], left.power + right.power,
                collect(first_index:last_index), left, right,
            )
            next_id[] += 1
            push!(internal, parent)
            parent
        end
        score = first(solve_interval(1, count))
        return build_interval(1, count), internal, score
    end

    full_mask = (1 << count) - 1
    mask_power = zeros(Float64, full_mask + 1)
    for mask in 1:full_mask
        bit = trailing_zeros(mask)
        mask_power[mask + 1] =
            mask_power[(mask & (mask - 1)) + 1] + powers[bit + 1]
    end
    memo = Dict{Int, Tuple{Float64, Int}}()
    function solve_subset(mask)
        haskey(memo, mask) && return memo[mask]
        count_ones(mask) == 1 && return (memo[mask] = (1.0, 0))
        first_bit = 1 << trailing_zeros(mask)
        best_score = -Inf
        best_split = 0
        submask = (mask - 1) & mask
        while submask > 0
            if !iszero(submask & first_bit)
                other = xor(mask, submask)
                fraction = mask_power[submask + 1] / mask_power[mask + 1]
                score = split_efficiency(fraction) * (
                    fraction * first(solve_subset(submask)) +
                    (1 - fraction) * first(solve_subset(other))
                )
                if score > best_score + 1e-14
                    best_score = score
                    best_split = submask
                end
            end
            submask = (submask - 1) & mask
        end
        memo[mask] = (best_score, best_split)
    end
    function build_subset(mask)
        if count_ones(mask) == 1
            leaf = trailing_zeros(mask) + 1
            return FeedNode(leaf, powers[leaf], leaf)
        end
        split = last(solve_subset(mask))
        left = build_subset(split)
        right = build_subset(xor(mask, split))
        parent = FeedNode(
            next_id[], left.power + right.power,
            sort(vcat(left.leaves, right.leaves)), left, right,
        )
        next_id[] += 1
        push!(internal, parent)
        parent
    end
    score = first(solve_subset(full_mask))
    build_subset(full_mask), internal, score
end

function collect_depths!(depths, node, depth=0)
    if isleaf(node)
        depths[only(node.leaves)] = depth
        return
    end
    collect_depths!(depths, node.left, depth + 1)
    collect_depths!(depths, node.right, depth + 1)
end

function collect_node_depths!(depths, node, depth=0)
    depths[node.id] = depth
    isleaf(node) && return
    collect_node_depths!(depths, node.left, depth + 1)
    collect_node_depths!(depths, node.right, depth + 1)
end

function transformer_length_mm(minimum_fraction)
    ratio = inv(minimum_fraction)
    ALUMINIUM_CP_M_S * log(ratio) /
    (8 * LOWER_FREQUENCY_HZ * ADIABATIC_LIMIT) * 1e3
end

function main()
    weights = selected_full_weights()
    power = abs2.(weights)
    normalized_power = power ./ sum(power)
    root, internal = balanced_tree(power)
    node_depth = Dict{Int, Int}()
    collect_node_depths!(node_depth, root)
    node_rows = NamedTuple[]
    minimum_fraction = 1.0
    for node in internal
        left_fraction = node.left.power / node.power
        right_fraction = node.right.power / node.power
        minimum_fraction = min(minimum_fraction, left_fraction, right_fraction)
        push!(node_rows, (;
            node_id=node.id,
            depth=node_depth[node.id],
            left_id=node.left.id,
            right_id=node.right.id,
            leaves=join(node.leaves, ';'),
            left_leaves=join(node.left.leaves, ';'),
            right_leaves=join(node.right.leaves, ';'),
            node_power_fraction=node.power / root.power,
            left_power_fraction=left_fraction,
            right_power_fraction=right_fraction,
            left_junction_width_mm=BRANCH_WIDTH_MM * left_fraction,
            right_junction_width_mm=BRANCH_WIDTH_MM * right_fraction,
            small_to_large_output_amplitude_ratio=
                sqrt(min(left_fraction, right_fraction) /
                     max(left_fraction, right_fraction)),
        ))
    end
    depths = zeros(Int, length(weights))
    collect_depths!(depths, root)
    max_depth = maximum(depths)
    minimum_width_mm = BRANCH_WIDTH_MM * minimum_fraction
    module_length_mm = transformer_length_mm(minimum_fraction)
    level_minimum_fraction = [minimum(
        min(row.left_power_fraction, row.right_power_fraction)
        for row in node_rows if row.depth == depth
    ) for depth in 0:(max_depth - 1)]
    level_transformer_length_mm = transformer_length_mm.(level_minimum_fraction)
    maximum_feed_axial_length_mm = sum(level_transformer_length_mm)
    minimum_feed_efficiency = (TARGET_GAIN / VERIFIED_SYSTEM_GAIN)^2
    worst_index = argmin([
        min(row.left_power_fraction, row.right_power_fraction) for row in node_rows
    ])
    worst = node_rows[worst_index]

    mkpath(OUTPUT_ROOT)
    nodes_path = joinpath(OUTPUT_ROOT, "passive_feed_tree_nodes.csv")
    open(nodes_path, "w") do io
        println(io, "node_id,depth,left_id,right_id,leaves,left_leaves,right_leaves,node_power_fraction,left_power_fraction,right_power_fraction,left_junction_width_mm,right_junction_width_mm,small_to_large_output_amplitude_ratio")
        for row in sort(node_rows; by=row -> row.node_id)
            println(io, join((
                row.node_id, row.depth, row.left_id, row.right_id, row.leaves,
                row.left_leaves, row.right_leaves, row.node_power_fraction,
                row.left_power_fraction, row.right_power_fraction,
                row.left_junction_width_mm, row.right_junction_width_mm,
                row.small_to_large_output_amplitude_ratio,
            ), ','))
        end
    end
    leaves_path = joinpath(OUTPUT_ROOT, "passive_feed_leaf_targets.csv")
    centers_mm = collect(-57.4:8.2:57.4)
    open(leaves_path, "w") do io
        println(io, "leaf,y_mm,target_pressure_weight,target_power_fraction,tree_depth,padding_stages")
        for index in eachindex(weights)
            println(io, join((
                index, centers_mm[index], weights[index], normalized_power[index],
                depths[index], max_depth - depths[index],
            ), ','))
        end
    end
    summary_path = joinpath(OUTPUT_ROOT, "passive_feed_summary.csv")
    open(summary_path, "w") do io
        println(io, "leaf_count,splitter_count,maximum_tree_depth,minimum_child_power_fraction,minimum_junction_width_mm,common_transformer_length_mm,level_minimum_power_fractions,level_transformer_lengths_mm,maximum_feed_axial_length_mm,worst_node_id,worst_small_power_fraction,worst_large_power_fraction,worst_target_amplitude_ratio,verified_ideal_feed_gain,minimum_total_feed_power_efficiency_for_gain_2")
        println(io, join((
            length(weights), length(internal), max_depth, minimum_fraction,
            minimum_width_mm, module_length_mm,
            join(level_minimum_fraction, ';'),
            join(level_transformer_length_mm, ';'),
            maximum_feed_axial_length_mm,
            worst.node_id,
            min(worst.left_power_fraction, worst.right_power_fraction),
            max(worst.left_power_fraction, worst.right_power_fraction),
            worst.small_to_large_output_amplitude_ratio,
            VERIFIED_SYSTEM_GAIN, minimum_feed_efficiency,
        ), ','))
    end
    design_path = joinpath(OUTPUT_ROOT, "passive_feed_design.jld2")
    JLD2.jldsave(
        design_path;
        format_version=1,
        full_pressure_weights=weights,
        normalized_leaf_power=normalized_power,
        leaf_depth=depths,
        maximum_tree_depth=max_depth,
        minimum_child_power_fraction=minimum_fraction,
        minimum_junction_width_mm=minimum_width_mm,
        common_transformer_length_mm=module_length_mm,
        level_minimum_power_fraction=level_minimum_fraction,
        level_transformer_length_mm,
        maximum_feed_axial_length_mm,
        minimum_total_feed_power_efficiency_for_gain_2=minimum_feed_efficiency,
        worst_node_id=worst.node_id,
        worst_small_power_fraction=
            min(worst.left_power_fraction, worst.right_power_fraction),
        worst_large_power_fraction=
            max(worst.left_power_fraction, worst.right_power_fraction),
        worst_target_amplitude_ratio=worst.small_to_large_output_amplitude_ratio,
    )
    println("[+] passive feed: $(length(weights)) leaves, $(length(internal)) splitters")
    println("[+] minimum child fraction=$minimum_fraction, width=$minimum_width_mm mm")
    println("[+] worst target amplitude ratio=$(worst.small_to_large_output_amplitude_ratio)")
    println("[+] max depth=$max_depth, level lengths=$level_transformer_length_mm mm")
    println("[+] maximum feed axial length=$maximum_feed_axial_length_mm mm")
    println("[+] feed power efficiency required for G>=2: $minimum_feed_efficiency")
    println("[+] $summary_path")
    println("[+] $design_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
