using DynamicEnergyBudgets
using DynamicEnergyBudgets: build_axis, setflux!, sumflux!, offset_apply!, setstate!
using AxisArrays
using Unitful

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "build axis arrays" begin
    time=0u"hr":1u"hr":10u"hr"
    a = build_axis([:one,:two], [:A,:B], time)
    @test axisnames(a) == (:time,)
    @test axisvalues(a) == (0u"hr":1u"hr":10u"hr",)
    @test axisnames(a[1]) == (:state, :transformations)
    @test axisvalues(a[1]) == (([:one, :two]), ([:A, :B]))
    @test a[1][:one,:A] == 0.0u"mol/hr"
    @test a[2][:one,:A] == 0.0u"mol/hr"
    a[1][:one,:A] = 10.0u"mol/hr"
    @test a[0.0u"hr"][:one,:A] == 10.0u"mol/hr"
    @test a[1.0u"hr"][:one,:A] == 0.0u"mol/hr"
end

@testset "setflux" begin
    o = Organism();
    or = o.nodes[1];
    @test or.J[1,1] == zero(or.J[1,1])
    @test or.Jrecord[1][1,1] == zero(or.J[1,1])

    or.J[1,1] = 10oneunit(or.J[1,1])
    @test or.J[1,1] == 10oneunit(or.J[1,1])
    @test or.Jrecord[1][1,1] == 10oneunit(or.J[1,1])

    apply(setflux!, o.nodes, 5u"hr")
    or.J[1,1] = 20oneunit(or.J[1,1])
    @test or.J[1,1] == 20oneunit(or.J[1,1])
    @test or.Jrecord[5u"hr"][1,1] == 20oneunit(or.J[1,1])
    or.J1[1,1] = 40oneunit(or.J1[1,1])
    @test or.J1[1,1] == 40oneunit(or.J[1,1])
    @test or.J1record[5u"hr"][1,1] == 40oneunit(or.J1[1,1])
end

@testset "sumflux" begin
    o = Organism();
    or1 = o.nodes[1];
    or2 = o.nodes[2];
    or1.J .= oneunit(or1.J[1,1]) .* [1,2,3,4,5,6]
    or2.J .= oneunit(or2.J[1,1]) .* [1,2,3,4,5,6] / 2
    du = fill(0.0u"mol/hr", 12)
    offset_apply!(sumflux!, du, o.nodes, 0)
    du
    @test du == [6,12,18,24,30,36,3,6,9,12,15,18] .* 1.0u"mol/hr"
end

@testset "setstate" begin
    u = fill(2.0u"mol", 12)
    o = Organism();
    or1 = o.nodes[1];
    or2 = o.nodes[2];
    @test or1.state[2] == 0.0001u"mol"
    offset_apply!(setstate!, o.nodes, u, 0)
    @test or1.state[1] == 2.0u"mol"
    @test or2.state[2] == 2.0u"mol"
    @test sum(or1.state) == 12.0u"mol"
    @test sum(or2.state) == 12.0u"mol"
end
