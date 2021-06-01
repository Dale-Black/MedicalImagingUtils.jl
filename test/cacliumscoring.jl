include("./imports.jl")
@testset ExtendedTestSet "calciumscoring" begin
    @testset ExtendedTestSet "calciumscoring" begin
        ii1 = [1 2 0 1
               1 5 0 1]
        ii2 = [1 2 0 1
               1 5 0 1]
        img = cat(ii1, ii2, dims=3)

        spc = [0.625, 0.625, 0.5]

        mm1 = [1 0 0 1
               1 1 0 1]
        mm2 = [1 0 0 1
               1 1 0 1]
        msk = cat(mm1, mm2, dims=3)
        answer = (0.6510416666666667, 1.953125)
        test = compute_calcium_scores(img, spc, msk)
        @test (test[1] ≈ answer[1]) && (test[2] ≈ answer[2])
    end
end