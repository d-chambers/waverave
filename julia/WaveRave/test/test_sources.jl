using WaveRave
using Test

using Debugger

@testset "Test get source time function" begin
    # test for 3 point stencil
    rick = ExplosiveRickerSource(location=[0,0])
    time = collect(0:0.1:100)
    stf = rick(time)
    @test length(stf) == length(time)
end
