using DynamicEnergyBudgets, Unitful, Test

using DynamicEnergyBudgets: tempcorr

reftemp = 310.0u"K"
arrtemp = 2000.0u"K"
lowerbound = 280.0u"K"
arrlower = 20000.0u"K"
upperbound = 315.0u"K"
arrupper = 70000.0u"K"

tc = TempCorrLowerUpper(
    reftemp=reftemp,
    arrtemp=arrtemp,
    tbelow=lowerbound,
    arrlower=arrlower,
    tabove=upperbound,
    arrupper=arrupper,
)

@test_broken tempcorr(tc, 100.0u"°C" |> u"K") ≈ 0.0 atol=1e-10
@test tempcorr(tc, -100.0u"°C" |> u"K") ≈ 0.0 atol=1e-10
@test tempcorr(tc, reftemp) == 1.0

tc = TempCorrLower(reftemp, arrtemp, lowerbound, arrlower)
@test tempcorr(tc, 100.0u"°C") > 1.0
@test tempcorr(tc, -100.0u"°C") ≈ 0.0 atol=1e-10
@test tempcorr(tc, reftemp) == 1.0

tc = TempCorr(reftemp, arrtemp)
@test tempcorr(tc, 100.0u"°C") > 1.0
@test tempcorr(tc, -100.0u"°C") > 0.0
@test tempcorr(tc, -100.0u"°C") < 1.0
@test tempcorr(tc, reftemp) == 1.0
