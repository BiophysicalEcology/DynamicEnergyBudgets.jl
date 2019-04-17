# Traits
struct HasE end
struct HasCN end
struct HasCNE end

has_reserves(o) = 
  if typeof(catabolism_pars(o)) <: AbstractCatabolismCN 
      HasCN()
  elseif typeof(catabolism_pars(o)) <: AbstractCatabolismCNE
      HasCNE()
  elseif typeof(catabolism_pars(o)) <: AbstractCatabolism
      HasE()
  else
      nothing
  end

