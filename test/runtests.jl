using SafeTestsets 

@safetestset "setup" begin include("setup.jl") end
@safetestset "temperature correction" begin include("tempcorrection.jl") end
@safetestset "math" begin include("math.jl") end
@safetestset "environment" begin include("environment.jl") end
@safetestset "balance" begin include("balance.jl") end
@safetestset "diffeq" begin include("diffeq.jl") end
