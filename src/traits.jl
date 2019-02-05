# Traits
struct HasCN end
struct HasCNE end

has_reserves(o) = 
  if typeof(catabolism_pars(o)) <: CatabolismCN 
      HasCN()
  elseif typeof(catabolism_pars(o)) <: CatabolismCNE
      HasCNE()
  else
      nothing
  end

