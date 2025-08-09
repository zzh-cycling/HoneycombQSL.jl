using Test

@testset "tenary" begin
    @test tern_to_base10(dit"120120011221112020;3") == 223302561
    @test ternary_pad(30, 12) == dit"000000001010;3"
end

@testset "neighbors" begin
    @test Set(get_neighbors(3, 3, :x)) == Set([(0, 3), (2, 5), (1, 4), (6, 9), (7, 10), (8, 11), (12, 15), (13, 16), (14, 17)])
    @test Set(get_neighbors(3, 3, :y)) == Set([(3, 6), (4, 7), (5, 8), (9, 12), (10, 13), (11, 14), (0, 15), (1, 16), (2, 17)])
    @test Set(get_neighbors(3, 3, :z)) == Set([(1, 3), (2, 4), (0, 5), (7, 9), (8, 10), (6, 11), (13, 15), (14, 16), (12, 17)])
end

@testset "lattice_translations" begin
    translation_order_list = lattice_translations(3, 2)
    correct_translations = [
        [8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7],
        [4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3],
        [1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10],
        [9, 8, 11, 10, 1, 0, 3, 2, 5, 4, 7, 6],
        [5, 4, 7, 6, 9, 8, 11, 10, 1, 0, 3, 2]
    ]
    @test Set(translation_order_list) == Set(correct_translations)
end