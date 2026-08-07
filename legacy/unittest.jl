# Archived from the original repository. These tests target NMFfirst.jl and are
# not part of the BilinearProjections.jl package test suite.

using Test

function run_tests()
    @testset "Data Encoding" begin
        setup()
        opendata("test_data.txt")
        initparam()
        readdata(1)
        encode(1)
        @test zp1[1, 1] ≈ -0.47043023366006696
        @test zp1[1, 2] ≈ -0.4730650187757215
        @test zp1[1, 3] ≈ -0.45970744255343265
    end

    @testset "Weight Matrix Update" begin
        setup()
        opendata("test_data.txt")
        initparam()
        readdata(1)
        encode(1)
        updateweight()
        @test wp2[1, 1] ≈ -0.02002024210989662
        @test wp2[2, 1] ≈ 0.013304136375350156
        @test wp2[3, 1] ≈ -0.006245227883797618
    end

    @testset "Data Decoding" begin
        setup()
        opendata("test_data.txt")
        initparam()
        readdata(1)
        encode(1)
        decode(1)
        @test zp2[1, 1] ≈ -0.005725636871695334
        @test zp2[1, 2] ≈ 0.004642767498101088
        @test zp2[1, 3] ≈ -0.00631154012968767
    end

    @testset "Error Calculation" begin
        setup()
        opendata("test_data.txt")
        initparam()
        readdata(1)
        encode(1)
        decode(1)
        @test calcerr(1) ≈ 0.006177029847716439
    end

    closedata()
    println("All tests passed!")
end

run_tests()
