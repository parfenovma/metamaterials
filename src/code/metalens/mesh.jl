using Gmsh: gmsh
using LinearAlgebra

const C_SOUND = 2340.0
const FREQ = 100_000.0
const N_CHANNELS = 5
const H_CHANNEL = 7.0
const W_CHANNEL = 17.0
const DOMAIN_LENGTH = 100.0
const FOCAL_LENGTH = 65.0

const DELAY_DATA = [
    (2.5, 11.3), (2.6, 12.0), (2.7, 12.3), (2.8, 13.0), (3.0, 14.3)
]

function get_best_amplitude(target_delay_us)
    best_A = DELAY_DATA[1][1]
    min_diff = Inf
    for (A, delay) in DELAY_DATA
        diff = abs(delay - target_delay_us)
        if diff < min_diff
            min_diff = diff
            best_A = A
        end
    end
    return best_A
end

function generate_crystal_points(A, y_offset)
    B = 0.4; CC = 4.0; yy = H_CHANNEL; xx = W_CHANNEL - 2.0*B; N_pts = 100
    C = A * CC / 100.0
    pts_bottom = Tuple{Float64, Float64}[]
    push!(pts_bottom, (0.0, 0.0))
    for i in 0:N_pts
        x = B + (4.0 * i * xx) / (5.0 * N_pts)
        angle_deg = (4.0 * 180.0 * i) / (N_pts + 1.0)
        if i <= N_pts/4 || i >= 3*N_pts/4
            y = (A / 2.0) * (1.0 - cosd(angle_deg))
        else
            y = (A / 2.0 - C / 2.0) * (1.0 - cosd(angle_deg)) + C
        end
        push!(pts_bottom, (x, y))
    end
    push!(pts_bottom, (xx + 2.0*B, 0.0)) 
    
    pts_top = Tuple{Float64, Float64}[]
    for p in pts_bottom
        push!(pts_top, (W_CHANNEL - p[1], yy - p[2]))
    end
    all_pts_mm = vcat(pts_bottom, pts_top)
    return [(x * 1e-3, (y + y_offset) * 1e-3) for (x, y) in all_pts_mm]
end


function build_metalens()
    gmsh.clear()
    gmsh.model.add("Acoustic_Metalens")
    
    println("=== Evaluating lens profile ...===")
    delay_center_us = DELAY_DATA[end][2]
    
    channel_surfaces = Int[]
    
    for i in 1:N_CHANNELS
        y_center_mm = (i - (N_CHANNELS + 1) / 2) * H_CHANNEL
        
        delta_dist_mm = sqrt(y_center_mm^2 + FOCAL_LENGTH^2) - FOCAL_LENGTH
        delta_time_us = (delta_dist_mm / C_SOUND) * 1000.0
        
        target_delay = delay_center_us - delta_time_us
        
        A = get_best_amplitude(target_delay)
        println("Channel $i | Y = $(round(y_center_mm, digits=1)) mm | Required: $(round(target_delay, digits=2)) microseconds | Chosen A = $A")
        
        y_offset = (i - 1) * H_CHANNEL
        pts = generate_crystal_points(A, y_offset)
        
        point_tags = [gmsh.model.occ.addPoint(p[1], p[2], 0.0) for p in pts]
        line_tags = Int[]
        for j in eachindex(point_tags)
            p1 = point_tags[j]
            p2 = point_tags[j == length(point_tags) ? 1 : j+1] 
            push!(line_tags, gmsh.model.occ.addLine(p1, p2))
        end
        
        loop_tag = gmsh.model.occ.addCurveLoop(line_tags)
        surf_tag = gmsh.model.occ.addPlaneSurface([loop_tag])
        push!(channel_surfaces, surf_tag)
    end
    
    total_lens_height_m = N_CHANNELS * H_CHANNEL * 1e-3
    free_space_surf = gmsh.model.occ.addRectangle(
        W_CHANNEL * 1e-3, 0.0, 0.0, 
        DOMAIN_LENGTH * 1e-3, total_lens_height_m
    )

    gmsh.model.occ.fragment([(2, free_space_surf)], [(2, s) for s in channel_surfaces])
    gmsh.model.occ.synchronize()
    
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    lines_dimtags = gmsh.model.getEntities(1)
    
    sources = Int[]; interfaces = Int[]
    walls = Int[]; free_space_boundaries = Int[]
    
    for (dim, tag) in lines_dimtags
        c_mass = gmsh.model.occ.getCenterOfMass(1, tag)
        x_c, y_c = c_mass[1], c_mass[2]
        
        if x_c < 1e-5
            push!(sources, tag)
        elseif abs(x_c - W_CHANNEL*1e-3) < 1e-5
            push!(interfaces, tag)
        elseif x_c > W_CHANNEL*1e-3 + 1e-5
            push!(free_space_boundaries, tag)
        else
            push!(walls, tag)
        end
    end
    
    gmsh.model.addPhysicalGroup(1, sources, 101); gmsh.model.setPhysicalName(1, 101, "Source")
    gmsh.model.addPhysicalGroup(1, walls, 102); gmsh.model.setPhysicalName(1, 102, "RigidWalls")
    gmsh.model.addPhysicalGroup(1, interfaces, 103); gmsh.model.setPhysicalName(1, 103, "LensInterface")
    gmsh.model.addPhysicalGroup(1, free_space_boundaries, 104); gmsh.model.setPhysicalName(1, 104, "RadiationBoundary")
    
    surfaces = [tag for (dim, tag) in gmsh.model.getEntities(2)]
    gmsh.model.addPhysicalGroup(2, surfaces, 201); gmsh.model.setPhysicalName(2, 201, "AcousticDomain")
    
    println("\n=== Mesh ===")
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", walls)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", 0.5 * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", 1.5 * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0 * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 20.0 * 1e-3)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    
    gmsh.model.mesh.generate(2)
    
    filename = "Acoustic_Metalens_F$(Int(FOCAL_LENGTH)).msh"
    gmsh.write(filename)
    println("  [+] Сетка линзы сохранена: $filename")
end

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
build_metalens()
# gmsh.fltk.run()
gmsh.finalize()
