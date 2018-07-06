using Revise
using DynamicEnergyBudgets
using DynamicEnergyBudgets: Records, split_state, build_record, build_flux, keep_records!, sum_flux!, offset_apply!, set_state!
using AxisArrays
using Unitful

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "build records" begin
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
end

# @testset "setflux" begin
#     o = Organism();
#     apply(keep_records!, o.organs, o.records, 1)
#     rec = o.records[1];
#     or = o.organs[1];
#     @test or.J[1,1] == zero(or.J[1,1])
#     @test rec.J[1][1,1] == zero(or.J[1,1])

#     or.J[1,1] = 10oneunit(or.J[1,1])
#     @test or.J[1,1] == 10oneunit(or.J[1,1])
#     @test rec.J[1][1,1] == 10oneunit(or.J[1,1])

#     apply(keep_records!, o.organs, o.records, 5u"hr")
#     or.J[1,1] = 20oneunit(or.J[1,1])
#     @test or.J[1,1] == 20oneunit(or.J[1,1])
#     @test rec.J[5u"hr"][1,1] == 20oneunit(or.J[1,1])
#     or.J1[1,1] = 40oneunit(or.J1[1,1])
#     @test or.J1[1,1] == 40oneunit(or.J[1,1])
#     @test rec.J1[5u"hr"][1,1] == 40oneunit(or.J1[1,1])
# end

@testset "sum_flux" begin
    du = fill(0.0u"mol/d", 12)
    o = Organism();
    or1 = o.organs[1];
    or2 = o.organs[2];
    or1.J .= oneunit(or1.J[1,1]) .* [1,2,3,4,5,6]
    or2.J .= oneunit(or2.J[1,1]) .* [1,2,3,4,5,6] / 2
    sum_flux!(du, o)
    du
    @test du == [6,12,18,24,30,36,3,6,9,12,15,18] .* u"mol/d" 
end

@testset "set_state" begin
    u = fill(2.0u"mol", 12)
    o = Organism();
    or1 = o.organs[1];
    or2 = o.organs[2];
    @test or1.state[2] == 0.0001u"mol"
    offset_apply!(set_state!, o.organs, u, 0)
    @test or1.state[1] == 2.0u"mol"
    @test or2.state[2] == 2.0u"mol"
    @test sum(or1.state) == 12.0u"mol"
    @test sum(or2.state) == 12.0u"mol"
end


@testset "split_state" begin
    o = Organism();
    u = fill(2.0u"mol", 12)
    us = split_state(o, u)
    @test us[1][1] == 2.0u"mol"
    @test us[2][2] == 2.0u"mol"
    @test sum(us[1]) == 12.0u"mol"
    @test sum(us[1]) == 12.0u"mol"
end
