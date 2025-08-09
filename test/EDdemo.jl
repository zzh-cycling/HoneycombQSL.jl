using Test
using HoneycombQSL

@testset "flux_path" begin
    m=3;n=4
    @test flux_path(1,1,m,n) ==([1, 2, 3, 8, 9, 10] .+1)
    @test flux_path(1,4,m,n) ==([7, 0, 1, 14, 15, 8] .+1)
    @test flux_path(2,4,m,n) ==([15, 8, 9, 22, 23, 16] .+1)
    @test flux_path(3,4,m,n) ==([23, 16, 17, 6, 7, 0] .+1)
    @test flux_path(3,1,m,n) == ([17, 18 ,19, 0, 1, 2] .+1)
    @test flux_path(3,2,m,n) == ([19, 20, 21, 2, 3, 4] .+1)
    @test flux_path(3,3,m,n) == ([21, 22, 23, 4, 5, 6] .+1)
end