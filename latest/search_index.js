var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#DynamicEnergyBudgets",
    "page": "Home",
    "title": "DynamicEnergyBudgets",
    "category": "module",
    "text": "This is a generalised DEB model. It was developed for plant modelling, but can potentially  model any organisms and symbioses.\n\nThis model can also be run in microclimates provided by the NicheMapr R package, and can use wide a range of photosynthesis and stomatal conductance formulations from Photosynthesis.jl.\n\nIt is also an in-progress attempt at using Julia\'s multiple-dispatch methods to abstract and generalise DEB theory and maintain a short, maintainable codebase for multiple models - potentially any organism.\n\nCode is adapted from the original DEBtool plant model by Bas Kooijman.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets-1",
    "page": "Home",
    "title": "DynamicEnergyBudgets",
    "category": "section",
    "text": "Modules = [DynamicEnergyBudgets]\nOrder   = [:module]"
},

{
    "location": "index.html#Index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Types-1",
    "page": "Home",
    "title": "Types",
    "category": "section",
    "text": "Modules = [DynamicEnergyBudgets]\nOrder   = [:type]"
},

{
    "location": "index.html#Functions-1",
    "page": "Home",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [DynamicEnergyBudgets]\nOrder   = [:function]"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractAssimilation",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractAssimilation",
    "category": "type",
    "text": "abstract AbstractAssimilation\n\nAssimilation \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractCarbonAssimilation",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractCarbonAssimilation",
    "category": "type",
    "text": "abstract AbstractCarbonAssimilation <: DynamicEnergyBudgets.AbstractAssimilation\n\nParaent of all Carbon assimilation types\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractNitrogenAssimilation",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractNitrogenAssimilation",
    "category": "type",
    "text": "abstract AbstractNitrogenAssimilation <: DynamicEnergyBudgets.AbstractAssimilation\n\nParaent of all Nitrogen assimilation types\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractScaling",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractScaling",
    "category": "type",
    "text": "abstract AbstractScaling\n\nSurface area scaling rules \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractStateFeedback",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractStateFeedback",
    "category": "type",
    "text": "abstract AbstractStateFeedback\n\nState feedback parameters. These modfy state based on state. \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractTempCorr",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractTempCorr",
    "category": "type",
    "text": "abstract AbstractTempCorr\n\nTemperature correction parameters\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Autophagy",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Autophagy",
    "category": "type",
    "text": "Autophagy. Parameters for self reabsorbtion when metabolic rates fall \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.C3Photosynthesis",
    "page": "Home",
    "title": "DynamicEnergyBudgets.C3Photosynthesis",
    "category": "type",
    "text": "Uses  FvCB photosynthesis model from Photosynthesis.jl \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.CarbonVars",
    "page": "Home",
    "title": "DynamicEnergyBudgets.CarbonVars",
    "category": "type",
    "text": "Variables for carbon assimilation \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.KooijmanArea",
    "page": "Home",
    "title": "DynamicEnergyBudgets.KooijmanArea",
    "category": "type",
    "text": "Surface areai scaling curve. Simulates growth and shade crowding later in life. \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.KooijmanSLAPhotosynthesis",
    "page": "Home",
    "title": "DynamicEnergyBudgets.KooijmanSLAPhotosynthesis",
    "category": "type",
    "text": "Parameters for simple photosynthesis module. With specific leaf area to convert area to mass \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Maturity",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Maturity",
    "category": "type",
    "text": "Maturity parameters. Seperated to make maturity modeling optional, reducing complexity \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.N_Assimilation",
    "page": "Home",
    "title": "DynamicEnergyBudgets.N_Assimilation",
    "category": "type",
    "text": "Parameters for lumped Nitrogen assimilation \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.NitrogenVars",
    "page": "Home",
    "title": "DynamicEnergyBudgets.NitrogenVars",
    "category": "type",
    "text": "Variables for nitgroen assimilation \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Organ",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Organ",
    "category": "type",
    "text": "type Organ{S, P, SH, V, F, F1}\n\nBasic model components. For a plants, organs might be roots, shoots and leaves \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Organism",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Organism",
    "category": "type",
    "text": "struct Organism{R, N, P}\n\nA single organism. Models can be run directly for an organism without an environment. \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Params",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Params",
    "category": "type",
    "text": "Model parameters that vary between organs \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Scenario",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Scenario",
    "category": "type",
    "text": "A scenario is an organism (or potentially organisms) growing in response to  real or simulated environmental data.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.SharedParams",
    "page": "Home",
    "title": "DynamicEnergyBudgets.SharedParams",
    "category": "type",
    "text": "Model parameters shared between organs \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.TempCorr",
    "page": "Home",
    "title": "DynamicEnergyBudgets.TempCorr",
    "category": "type",
    "text": "Simple temperature correction parameters \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.TempCorrLower",
    "page": "Home",
    "title": "DynamicEnergyBudgets.TempCorrLower",
    "category": "type",
    "text": "Temperature correction with lower boudn parameters\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.TempCorrLowerUpper",
    "page": "Home",
    "title": "DynamicEnergyBudgets.TempCorrLowerUpper",
    "category": "type",
    "text": "Temperature correction with lower and upper bound parameters\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Vars",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Vars",
    "category": "type",
    "text": "Model variables \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractAllometry",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractAllometry",
    "category": "type",
    "text": "abstract AbstractAllometry\n\nAllometry. Scaling rules to relate size to mass. \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.AbstractNH4_NO3Assimilation",
    "page": "Home",
    "title": "DynamicEnergyBudgets.AbstractNH4_NO3Assimilation",
    "category": "type",
    "text": "abstract AbstractNH4_NO3Assimilation <: DynamicEnergyBudgets.AbstractNitrogenAssimilation\n\nParent of Ammonia/Nitrate assimilation types\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Kooijman_NH4_NO3Assimilation",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Kooijman_NH4_NO3Assimilation",
    "category": "type",
    "text": "Parameters for Ammonia/Nitrate assimilation \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.Records",
    "page": "Home",
    "title": "DynamicEnergyBudgets.Records",
    "category": "type",
    "text": "struct Records{V, F, F1}\n\nRecords of variables and flux for ploting and analysis \n\n\n\n"
},

{
    "location": "index.html#Parameters-and-Types-1",
    "page": "Home",
    "title": "Parameters and Types",
    "category": "section",
    "text": "DynamicEnergyBudgets uses the julia type system to organise  both the parameters and the logic of the model. This means the model is completely modular and customisable.Modules = [DynamicEnergyBudgets]\nOrder   = [:type]\nPages   = [\"types.jl\"]"
},

{
    "location": "index.html#DynamicEnergyBudgets.debmodel!-Tuple{Any,Number}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.debmodel!",
    "category": "method",
    "text": "A generalised multi-reserve, multi-organ Dynamic Energy Budget model.\n\nApplies metabolism, translocation and assimilation mehtods to N organs.\n\nsettings is a struct with required model data, DEBSettings or similar. t is the timestep\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.runmodel!-Tuple{Any,Any,DynamicEnergyBudgets.Organism,Number}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.runmodel!",
    "category": "method",
    "text": "Method to run a DEB organism model, with a static environment.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.runmodel!-Tuple{Any,Any,DynamicEnergyBudgets.Scenario,Number}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.runmodel!",
    "category": "method",
    "text": "Method to run a DEB scenario, with environment applied to organism(s).\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.allometric_height-Tuple{DynamicEnergyBudgets.SqrtAllometry,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.allometric_height",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.catabolism!-Tuple{Any,DynamicEnergyBudgets.AbstractStateCNE,Number}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.catabolism!",
    "category": "method",
    "text": "catabolism!(o, u::AbstractStateCNE, t::Number)\n\nCatabolism for E, C and N, or C, N and E reserves. Does not finalise flux in J - operates only on J1 (intermediate storage)\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.dissipation!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.dissipation!",
    "category": "method",
    "text": "Dissipation for any reserve. Growth, maturity and maintenence are grouped as dissipative processes.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.feedback!-Tuple{Any,DynamicEnergyBudgets.Autophagy,DynamicEnergyBudgets.AbstractState}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.feedback!",
    "category": "method",
    "text": "Function to apply feedback on growth the process, such as autopagy in resource shortage.\n\nWithout a function like this you will likely be occasionally breaking the  laws of thermodynamics by introducing negative rates.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.find_rate-Union{Tuple{Any,Tuple{Tuple{Vararg{T,N}} where T,Tuple{Vararg{T,N}} where T,Vararg{Any,N} where N}}, Tuple{N}} where N",
    "page": "Home",
    "title": "DynamicEnergyBudgets.find_rate",
    "category": "method",
    "text": "Calculate rate formula. TODO: use Roots.jl for this\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.germinated-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.germinated",
    "category": "method",
    "text": "Check if germination has happened. Independent for each organ, although this may not make sense.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.growth!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.growth!",
    "category": "method",
    "text": "Allocates reserves to growth.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.maintenence!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.maintenence!",
    "category": "method",
    "text": "Allocates reserve drain due to maintenance.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.maturity!-Tuple{Any,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.maturity!",
    "category": "method",
    "text": "maturity!(f, o, u)\n\nAllocates reserve drain due to maturity maintenance. Stores in M state variable if it exists.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.metabolism!-Tuple{Any,Number}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.metabolism!",
    "category": "method",
    "text": "Metabolism is an identical process for all organs, with potentially different parameters or area and rate functions.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.product!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.product!",
    "category": "method",
    "text": "Allocates waste products from growth and maintenance.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.reserve_drain!-NTuple{4,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.reserve_drain!",
    "category": "method",
    "text": "Generalised reserve drain for any flux column col (ie :gro) and any combination of reserves.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.reserve_loss!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.reserve_loss!",
    "category": "method",
    "text": "Generalised reserve loss to track carbon. \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.reuse_rejected!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.reuse_rejected!",
    "category": "method",
    "text": "Reallocate state rejected from synthesizing units.\n\nTODO add a 1-organs method Also how does this interact with assimilation?\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.scaling-Tuple{DynamicEnergyBudgets.KooijmanArea,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.scaling",
    "category": "method",
    "text": "\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.translocate!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.translocate!",
    "category": "method",
    "text": "Versions for E, CN and CNE reserves.\n\nTranslocation is occurs between adjacent organs.  This function is identical both directiono, so on represents whichever is not the current organs. Will not run with less than 2 organs.\n\nFIXME this will be broken for organs > 2\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.translocation!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.translocation!",
    "category": "method",
    "text": "Versions for E, CN and CNE reserves.\n\nTranslocation is occurs between adjacent organs.  This function is identical both directiono, so on represents whichever is not the current organs. Will not run with less than 2 organs.\n\nFIXME this will be broken for organs > 2\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.κtra-Tuple{Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.κtra",
    "category": "method",
    "text": "κtra is the difference paramsbetween κsoma and κrep\n\n\n\n"
},

{
    "location": "index.html#Model-1",
    "page": "Home",
    "title": "Model",
    "category": "section",
    "text": "Core functions of the model.Modules = [DynamicEnergyBudgets]\nOrder   = [:function]\nPages   = [\"model.jl\"]"
},

{
    "location": "index.html#DynamicEnergyBudgets.correct_temps!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.correct_temps!",
    "category": "method",
    "text": "Scale variables by temperature\n\n\n\n"
},

{
    "location": "index.html#Environment-1",
    "page": "Home",
    "title": "Environment",
    "category": "section",
    "text": "Modules = [DynamicEnergyBudgets]\nOrder   = [:function]\nPages   = [\"environment.jl\"]"
},

{
    "location": "index.html#DynamicEnergyBudgets.assimilation!-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.assimilation!",
    "category": "method",
    "text": "assimilation!(o1, o2)\n\nRuns assimilation methods, depending on formulation and state.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.assimilation!-Tuple{DynamicEnergyBudgets.AbstractCarbonAssimilation,Any,Any,DynamicEnergyBudgets.AbstractStateCNE}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.assimilation!",
    "category": "method",
    "text": "assimilation!(f::AbstractCarbonAssimilation, o1, o2, u::AbstractStateCNE)\n\nRuns nitrogen uptake, and combines N with translocated C.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.assimilation!-Tuple{DynamicEnergyBudgets.AbstractNH4_NO3Assimilation,Any,Any,DynamicEnergyBudgets.AbstractStateCNE}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.assimilation!",
    "category": "method",
    "text": "assimilation!(f::AbstractNH4_NO3Assimilation, o1, o2, u::AbstractStateCNE)\n\nRuns nitrogen uptake for nitrate and ammonia, and combines N with translocated C. Unused ammonia is discarded.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.assimilation!-Tuple{DynamicEnergyBudgets.N_Assimilation,Any,Any,DynamicEnergyBudgets.AbstractStateCNE}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.assimilation!",
    "category": "method",
    "text": "assimilation!(f::N_Assimilation, o1, o2, u::AbstractStateCNE)\n\nRuns nitrogen uptake, and combines N with translocated C.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.photosynthesis-Tuple{DynamicEnergyBudgets.C3Photosynthesis,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.photosynthesis",
    "category": "method",
    "text": "photosynthesis(f::C3Photosynthesis, o1, o2)\n\nReturns carbon assimilated in mols per time.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.photosynthesis-Tuple{DynamicEnergyBudgets.KooijmanSLAPhotosynthesis,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.photosynthesis",
    "category": "method",
    "text": "photosynthesis(f::KooijmanSLAPhotosynthesis, o1, o2)\n\nReturns carbon assimilated in mols per time.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.uptake_nitrogen-Tuple{DynamicEnergyBudgets.Kooijman_NH4_NO3Assimilation,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.uptake_nitrogen",
    "category": "method",
    "text": "uptake_nitrogen(f::Kooijman_NH4_NO3Assimilation, o1, o2)\n\nReturns total nitrogen, nitrate and ammonia assimilated in mols per time.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.uptake_nitrogen-Tuple{DynamicEnergyBudgets.N_Assimilation,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.uptake_nitrogen",
    "category": "method",
    "text": "uptake_nitrogen(f::N_Assimilation, o1, o2)\n\nReturns nitrogen assimilated in mols per time.\n\n\n\n"
},

{
    "location": "index.html#Assimilation-1",
    "page": "Home",
    "title": "Assimilation",
    "category": "section",
    "text": "Carbon and nitrogen assimilation functions. These are the main points of difference between organs.Modules = [DynamicEnergyBudgets]\nOrder   = [:function]\nPages   = [\"assimilation.jl\"]"
},

{
    "location": "index.html#DynamicEnergyBudgets.catabolic_fluxes-Tuple{Any,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.catabolic_fluxes",
    "category": "method",
    "text": "catabolic_fluxes(ureserve, A_turnover, r)\n\n\nReturns the current catabolic flux at rate r, or the flux as a proportion of u[V], depending on ureserve values.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.half_saturation-Tuple{Any,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.half_saturation",
    "category": "method",
    "text": "half_saturation(max, half, x)\n\n\nHalf satration curve.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.rate_formula-Tuple{Any,Tuple{T,T,T} where T,Tuple{T,T,T} where T,Any,Any,Any,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.rate_formula",
    "category": "method",
    "text": "rate_formula(r, ureserve::NTuple, A_turnover::NTuple, j_E_mai, y_V_E, κsoma)\n\nRate formulas for E, CN or CNE reserves\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.stoich_merge-Tuple{Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.stoich_merge",
    "category": "method",
    "text": "stoich_merge(a, b)\n\n\nMerge two inputs stoichiometrically. The minimum value is limiting, and stochasticity of pairing is simulated so that for any a, b  stoich_merge(a, b) < a, to the limit b → ∞ where stoich_merge(a, b) = a   \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.synthesizing_unit-NTuple{4,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.synthesizing_unit",
    "category": "method",
    "text": "synthesizing_unit(Ja, Jb, ya, yb)\n\n\nMerge fluxes stoichiometrically into general reserve Eab based on yeild  fractions ya and yb. The remainder is returned as unmixed reserves Ea and Eb.\n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.tempcorr-Tuple{Any,Any,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.tempcorr",
    "category": "method",
    "text": "tempcorr(T, T1, A, [L, AL,] [H, AH])\n\nDEB tempcorr function. Uses lower and uppper bounds if they are supplied. Temperatures all in Kelvins.\n\n\n\n"
},

{
    "location": "index.html#Functions-2",
    "page": "Home",
    "title": "Functions",
    "category": "section",
    "text": "Low-level DEB theory functionsModules = [DynamicEnergyBudgets]\nOrder   = [:function]\nPages   = [\"functions.jl\"]"
},

{
    "location": "index.html#DynamicEnergyBudgets.set_cur_records!-Tuple{Any,DynamicEnergyBudgets.Records,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.set_cur_records!",
    "category": "method",
    "text": "update vars and flux to record for time t \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.set_cur_records!-Tuple{Any,Void,Any}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.set_cur_records!",
    "category": "method",
    "text": "Void records so we just stay on the same record \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.setstate!-Tuple{Any,Any,Int64}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.setstate!",
    "category": "method",
    "text": "copy the diffeq state to organs \n\n\n\n"
},

{
    "location": "index.html#DynamicEnergyBudgets.sumflux!-Tuple{Any,Any,Int64}",
    "page": "Home",
    "title": "DynamicEnergyBudgets.sumflux!",
    "category": "method",
    "text": "sum flux matrix \n\n\n\n"
},

{
    "location": "index.html#Setup-and-utilities-1",
    "page": "Home",
    "title": "Setup and utilities",
    "category": "section",
    "text": "Modules = [DynamicEnergyBudgets]\nOrder   = [:function]\nPages   = [\"apply.jl\"]\nPages   = [\"setup.jl\"]"
},

]}
