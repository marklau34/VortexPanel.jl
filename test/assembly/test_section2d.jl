# The two-dimensional referee of validation stage 6 (validation/common/panel2d.jl): a linear-vortex panel method,
# checked against the exact Joukowski airfoil, that gives the NACA 0015 section the wings' midspan is judged by.

@testset "Joukowski airfoil" begin
    err = map((50, 100, 200)) do n
        pts, slope = joukowski_section(n)
        abs(panel2d(pts, 5.0).cl_gamma / (slope * sind(5.0)) - 1)
    end
    @test err[3] < 3e-5
    @test err[1] / err[2] > 3.5 && err[2] / err[3] > 3.5 # second order
end

@testset "NACA 0015 section" begin
    r = panel2d(naca00_section(200), 5.0)
    @test r.cl_gamma ≈ 0.617012 rtol = 1e-4 # the converged value (1600 panels per surface)
    @test abs(r.cl / r.cl_gamma - 1) < 1e-4 # pressure and circulation agree at a sharp trailing edge
    @test abs(r.cd) < 1e-4
    @test r.cm ≈ -0.009165 atol = 2e-5
    r0 = panel2d(naca00_section(100), 0.0)
    @test abs(r0.cl_gamma) < 1e-12 && abs(r0.cm) < 1e-12
end
