using DynamicEnergyBudgets, Unitful, Test

reftemp = 310.0u"K"
arrtemp = 2000.0u"K"
lowerbound = 280.0u"K"
arrlower = 20000.0u"K"
upperbound = 315.0u"K"
arrupper = 70000.0u"K"

tc = TempCorrLowerUpper(reftemp, arrtemp, lowerbound, arrlower, upperbound, arrupper)
@test tempcorr(100.0u"°C", tc) ≈ 0.0 atol=1e-10
@test tempcorr(-100.0u"°C", tc) ≈ 0.0 atol=1e-10
@test tempcorr(reftemp, tc) == 1.0

tc = TempCorrLower(reftemp, arrtemp, lowerbound, arrlower)
@test tempcorr(100.0u"°C", tc) > 1.0
@test tempcorr(-100.0u"°C", tc) ≈ 0.0 atol=1e-10
@test tempcorr(reftemp, tc) == 1.0

tc = TempCorr(reftemp, arrtemp)
@test tempcorr(100.0u"°C", tc) > 1.0
@test tempcorr(-100.0u"°C", tc) > 0.0
@test tempcorr(-100.0u"°C", tc) < 1.0
@test tempcorr(reftemp, tc) == 1.0
