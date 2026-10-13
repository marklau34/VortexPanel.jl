using SparseArrays
using VortexPanel.DOFMap

"""Rows of C sum to one and rows of K to zero"""
function check_row_sums(dm)
    @test maximum(abs.(sum(dm.C, dims=2) .- 1)) < 1e-12
    @test isempty(dm.K) || maximum(abs.(sum(dm.K, dims=2))) < 1e-12
end

linear_field(x) = 0.3 + 1.7x[1] - 0.9x[2] + 2.3x[3]

const WING_G = build_geometry(reference_wing(1))
const WING_DM = build_dofmap(WING_G)

@testset "clean wing" begin
    g, dm = WING_G, WING_DM
    @test dm.nslot == 702 && dm.nq == 702 && count(>(0), dm.lower_slot) == 18
    @test dm.C == sparse(I, 702, 702)
    check_row_sums(dm)
    @test size(dm.K) == (18, 702)
    # each Kutta row is +1 on the trailing edge node's slot and -1 on its lower copy
    for k in 1:18
        v = dm.wake_te.node[k]
        cols, vals = findnz(dm.K[k, :])
        @test Dict(zip(cols, vals)) == Dict(dm.q_of_slot[v] => 1.0, dm.q_of_slot[dm.lower_slot[v]] => -1.0)
    end
    # lower surface corners at trailing edge nodes take the copies, every other corner its node
    for t in eachindex(g.tris), k in 1:3
        v, s = g.tris[t][k], dm.corner_slot[t][k]
        @test s == (dm.lower_slot[v] > 0 && g.kutta.corner_side[t][k] == TELower ? dm.lower_slot[v] : v)
    end
    # one trailing edge line along +y with a downstream bisector
    w = dm.wake_te
    @test w.lines == [collect(1:18)] && w.closed == [false]
    @test issorted([p[2] for p in w.x])
    @test all(b -> b ≈ SVector(1.0, 0.0, 0.0), w.bisector)
    @test all(t -> t ≈ SVector(0.0, 1.0, 0.0), w.tangent)
    # a uniform doublet: uniform on every slot, nothing shed
    @test dm.C * ones(dm.nq) ≈ ones(dm.nslot)
    @test norm(dm.K * ones(dm.nq)) < 1e-12
    @test sum(dm.slot_area) ≈ sum(g.quality.area)
    rows = bordered_rows(dm)
    @test length(rows) == 1 && rows[1].e == ones(dm.nq) && sum(rows[1].c) ≈ sum(g.quality.area)
    @test reproduction_error(dm, g, linear_field) < 1e-12
    files = MeshLibrary.write_vtk(dm, g, joinpath(mktempdir(), "wing"))
    @test length(files) == 2 && all(isfile, files)
end

@testset "cap trailing edge copy" begin
    g = WING_G
    dl = build_dofmap(g; opts=DOFMapOptions(cap_te_copy=:lower))
    for t in findall(==(CapTri), g.tri_kind), k in 1:3
        v = g.tris[t][k]
        @test WING_DM.corner_slot[t][k] == v
        @test dl.corner_slot[t][k] == (WING_DM.lower_slot[v] > 0 ? WING_DM.lower_slot[v] : v)
    end
    @test_throws ArgumentError build_dofmap(g; opts=DOFMapOptions(cap_te_copy=:both))
end

@testset "Galerkin reduction" begin
    dm = WING_DM
    te_of_wake_node = [1, 1, 2, 18]
    W = wake_rows(dm, te_of_wake_node)
    @test W == dm.K[te_of_wake_node, :]
    A = sprand(dm.nslot, dm.nslot + 4, 0.01)
    b = rand(dm.nslot)
    Ar, br = reduce_system(dm, A, b; W=W)
    @test size(Ar) == (dm.nq, dm.nq) && length(br) == dm.nq
    @test Ar ≈ dm.C' * A * [dm.C; W]
    @test_throws DimensionMismatch reduce_system(dm, A, b)
end

@testset "T-junctions" begin
    X, faces = quad_cube(split_face=true)
    g = build_geometry(X, faces)
    dm = build_dofmap(g)
    @test dm.nslot == length(g.x) && count(==(:tjunction), dm.constraint) == 2
    @test dm.nq == dm.nslot - 2
    @test Set(dm.slot_node[s] for s in 1:dm.nslot if dm.constraint[s] == :tjunction) == Set(tj.hanging for tj in g.tjunctions)
    check_row_sums(dm)
    @test reproduction_error(dm, g, linear_field) < 1e-12
    @test stitch_offset(dm, g) < 1e-12
    @test isempty(dm.K) && length(bordered_rows(dm)) == 1
    @test length(MeshLibrary.write_vtk(dm, g, joinpath(mktempdir(), "cube"))) == 1 # no trailing edge
end

@testset "stitched non-matching seams" begin
    X, faces = seam_cube()
    g = build_geometry(X, faces)
    @test count(p -> p.kind == :gap, g.gaps) == 2 && isempty(g.tjunctions) && isempty(g.openings)
    dm = build_dofmap(g)
    @test count(==(:gap), dm.constraint) == 4 # the two interior nodes of each thirds side follow the halves side
    @test all(s -> g.x[dm.slot_node[s]][1] in (1 / 3, 2 / 3), findall(==(:gap), dm.constraint))
    check_row_sums(dm)
    @test reproduction_error(dm, g, linear_field) < 1e-12 # exact through coincident stitches
    @test stitch_offset(dm, g) < 1e-12
    @test g.bodies[1].closed && g.bodies[1].volume ≈ 1
end

@testset "stitched offset gap" begin
    mesh = cube_sphere(4)
    top = mesh.elemset[findfirst(s -> s.name == "face_3+", mesh.elemset)].elems
    ε = 0.005
    for v in unique(reduce(vcat, [collect(mesh.elements[e].con) for e in top]))
        mesh.nodes[v] = MeshLibrary.Point((1 + ε) * mesh.nodes[v].coord)
    end
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    @test count(==(:gap), dm.constraint) == 16 # one closed side of 16 nodes follows the other
    check_row_sums(dm)
    grad = norm(SVector(1.7, -0.9, 2.3))
    @test reproduction_error(dm, g, linear_field) <= grad * g.gaps[1].width * (1 + 1e-9)
    @test stitch_offset(dm, g) ≈ g.gaps[1].width rtol = 0.05
end

@testset "non-matching trailing edge" begin
    mesh = reference_wing(1; split_te=true)
    copies = 685:702 # PAGE's lower trailing edge copies, shifted along the span except at the tips
    dy = minimum(diff(sort(unique(round.(getindex.(mesh.nodes.coord[copies], 2), digits=12)))))
    for v in copies
        p = mesh.nodes[v].coord
        abs(abs(p[2]) - WING_SPAN / 2) < 1e-9 && continue
        mesh.nodes[v] = MeshLibrary.Point(p + SVector(0.0, 0.3dy, 0.0))
    end
    g = build_geometry(mesh)
    dm = build_dofmap(g)
    @test length(g.kutta.pair_chains) == 1
    @test count(>(0), dm.lower_slot) == 2 # only the shared tip nodes are duplicated
    @test size(dm.K, 1) == 34 && dm.wake_te.lines == [collect(1:34)]
    check_row_sums(dm)
    # different linear doublet fields on the two surfaces: K gives their difference at every union node
    f_up(x) = 1.0 + 0.4x[2] + 0.2x[1]
    f_low(x) = -0.5 + 1.1x[2] - 0.3x[1]
    x = g.x[dm.slot_node]
    q = [dm.slot_side[s] == TELower ? f_low(x[s]) : f_up(x[s]) for s in dm.free_slot]
    μw = dm.K * q
    @test maximum(abs(μw[k] - (f_up(p) - f_low(p))) for (k, p) in enumerate(dm.wake_te.x)) < 1e-12
end

@testset "zipped tips" begin
    g = build_geometry(reference_wing(1); opts=GeometryOptions(openings=:zip))
    dm = build_dofmap(g)
    @test count(==(:zip), dm.constraint) == 2 * 18 # interior nodes of one side of each tip
    @test count(>(0), dm.lower_slot) == 18 && size(dm.K, 1) == 18 # the trailing edge keeps its copies
    check_row_sums(dm)
    @test reproduction_error(dm, g, x -> 0.3 + 1.7x[1] - 0.9x[2]) < 1e-12 # the sides pair at equal x
    @test stitch_offset(dm, g) ≈ maximum(p.width for p in g.gaps if p.kind == :zip) rtol = 1e-9
    @test !dm.body_closed[1] # a zipped body is not closed, but a uniform μ on it still induces nothing
    @test length(bordered_rows(dm)) == 1
end

@testset "abused wing" begin
    wing = reference_wing(1; split_te=false)
    project = naca_projector(chord=WING_CHORD, thickness=WING_THICKNESS)
    m, log = abuse(wing, (jitter=0.2, gaps=(2, 0.01), tjunctions=3, project=project); seed=11)
    g = build_geometry(m)
    dm = build_dofmap(g)
    @test count(==(:tjunction), dm.constraint) == length(g.tjunctions) == 6
    @test count(==(:gap), dm.constraint) == 4 # the two interior nodes of one side of each 3-edge slit
    @test dm.nq == dm.nslot - 10
    @test size(dm.K, 1) == 18
    check_row_sums(dm)
    grad = norm(SVector(1.7, -0.9, 2.3))
    @test reproduction_error(dm, g, linear_field) <= grad * maximum(p.width for p in g.gaps) * (1 + 1e-9)
    @test isempty(dm.log)
end
