using Test, DynamicEnergyBudgets, Unitful
using DynamicEnergyBudgets: Records, define_organs, split_state, sum_flux!

# @testset "build records" begin
    # time=0u"hr":1u"hr":10u"hr"
    # a = build_record(build_flux(1.0u"mol/hr", [:one,:two], [:A,:B]), time)
    # @test axisnames(a) == (:time,)
    # @test axisvalues(a) == (0u"hr":1u"hr":10u"hr",)
    # @test axisnames(a[1]) == (:state, :transformations)
    # @test axisvalues(a[1]) == (([:one, :two]), ([:A, :B]))
    # @test a[1][:one,:A] == 0.0u"mol/hr"
    # @test a[2][:one,:A] == 0.0u"mol/hr"
    # a[1][:one,:A] = 10.0u"mol/hr"
    # @test a[0.0u"hr"][:one,:A] == 10.0u"mol/hr"
    # @test a[1.0u"hr"][:one,:A] == 0.0u"mol/hr"
    # o = Organ()
    # Records(o, time)
# end

@testset "sum_flux" begin
    o = Plant();
    or1, or2 = define_organs(o, 1);

    or1.J .= oneunit(or1.J[1,1])
    or2.J .= oneunit(or2.J[1,1]) * 2
    du = fill(0.0u"mol/hr", 6)
    sum_flux!(du, (or1, or2))
    @test du == [5.0, 5.0, 5.0, 10.0, 10.0, 10.0] .* u"mol/hr" 
end

@testset "split_state" begin
    o = Plant();
    os = define_organs(o, 1);
    u = [1.0:6.0...]u"mol"
    us = split_state(os, u)
    @test us[1][1] == 1.0u"mol"
    @test us[2][1] == 4.0u"mol"
    @test sum(us[1]) == sum(1:3) .* u"mol"
    @test sum(us[2]) == sum(4:6) .* u"mol"
end
