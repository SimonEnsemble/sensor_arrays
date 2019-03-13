using SVD

@testset "Sensor Tests" begin
    H = make_h_matrix(3, 5, "CO2", "C2H6")
    G = [(8.91e-5) (4.84e-5); (1.87e-5) (4.68e-5)] .- [5.84e-6; 9.56e-6]
    @test isapprox(H, G)
end
