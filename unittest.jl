using Test

# Define a function to run the tests
function run_tests()
    @testset "Data Encoding" begin
        # Initialize the variables
        setup()
        opendata("test_data.txt")
        initparam()
        
        # Read the test data
        readdata(1)
        
        # Encode the data
        encode(1)
        
        # Check the encoded values
        @test zp1[1, 1] ≈ -0.47043023366006696
        @test zp1[1, 2] ≈ -0.4730650187757215
        @test zp1[1, 3] ≈ -0.45970744255343265
    end
    
    @testset "Weight Matrix Update" begin
        # Initialize the variables
        setup()
        opendata("test_data.txt")
        initparam()
        
        # Read the test data
        readdata(1)
        
        # Encode the data
        encode(1)
        
        # Update the weight matrix
        updateweight()
        
        # Check the updated weight matrix
        @test wp2[1, 1] ≈ -0.02002024210989662
        @test wp2[2, 1] ≈ 0.013304136375350156
        @test wp2[3, 1] ≈ -0.006245227883797618
    end
    
    @testset "Data Decoding" begin
        # Initialize the variables
        setup()
        opendata("test_data.txt")
        initparam()
        
        # Read the test data
        readdata(1)
        
        # Encode the data
        encode(1)
        
        # Decode the data
        decode(1)
        
        # Check the decoded values
        @test zp2[1, 1] ≈ -0.005725636871695334
        @test zp2[1, 2] ≈ 0.004642767498101088
        @test zp2[1, 3] ≈ -0.00631154012968767
    end
    
    @testset "Error Calculation" begin
        # Initialize the variables
        setup()
        opendata("test_data.txt")
        initparam()
        
        # Read the test data
        readdata(1)
        
        # Encode the data
        encode(1)
        
        # Decode the data
        decode(1)
        
        # Calculate the error
        err = calcerr(1)
        
        # Check the error value
        @test err ≈ 0.006177029847716439
    end
    
    # Add more tests if needed
    
    # Close the data file
    closedata()
    
    # Print completion message
    println("All tests passed!")
end

# Run the tests
run_tests()
