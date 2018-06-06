abstract type AbstractState{N,T} <: FieldVector{N,T} end
abstract type AbstractStateE{N,T} <: AbstractState{N,T} end
abstract type AbstractStateCN{N,T} <: AbstractState{N,T} end
abstract type AbstractStateCNE{N,T} <: AbstractState{N,T} end

@traitdef StateHasM{X}
@traitdef StateHasP{X}
SimpleTraits.trait{X1}(::Type{StateHasM{X1}}) = :M in fieldnames(X1) ? StateHasM{X1} : Not{StateHasM{X1}}
SimpleTraits.trait{X1}(::Type{StateHasP{X1}}) = :P in fieldnames(X1) ? StateHasP{X1} : Not{StateHasP{X1}}

@mix @label @with_kw struct P{T} P::T = 0.0u"mol"  | "Production"       end
@mix @label @with_kw struct V{T} V::T = 1e-4u"mol" | "Structure"        end
@mix @label @with_kw struct M{T} M::T = 0.0u"mol"  | "Maturity"         end
@mix @label @with_kw struct C{T} C::T = 1e-4u"mol" | "Carbon Reserve"   end
@mix @label @with_kw struct N{T} N::T = 1e-4u"mol" | "Nitrogen Reserve" end
@mix @label @with_kw struct E{T} E::T = 10.0u"mol" | "General Reserve"  end

@E mutable struct StateE{} <: AbstractStateE{1,T} end
@C @N mutable struct StateCN{} <: AbstractStateCN{2,T} end
@C @N @E mutable struct StateCNE{} <: AbstractStateCNE{3,T} end
@V mutable struct StateV{} <: AbstractStateE{1,T}  end
@V @E mutable struct StateVE{} <: AbstractStateE{2,T} end
@V @C @N mutable struct StateVCN{} <: AbstractStateCN{3,T} end
@V @C @N @E mutable struct StateVCNE{} <: AbstractStateCN{3,T} end
@V @M @E mutable struct StateVME{} <: AbstractStateE{3,T} end
@V @M @C @N mutable struct StateVMCN{} <: AbstractStateE{3,T} end
@V @M @C @N @E mutable struct StateVMCNE{} <: AbstractStateE{3,T} end
@P @V mutable struct StatePV{} <: AbstractState{2,T} end
@P @V @E mutable struct StatePVE{} <: AbstractStateE{3,T} end
@P @V @C @N mutable struct StatePVCN{} <: AbstractStateCN{4,T} end
@P @V @C @N @E mutable struct StatePVCNE{} <: AbstractStateCNE{5,T} end
@P @V @M @E mutable struct StatePVME{} <: AbstractStateE{4,T} end
@P @V @M @C @N mutable struct StatePVMCN{} <: AbstractStateCN{5,T} end
@P @V @M @C @N @E mutable struct StatePVMCNE{} <: AbstractStateCNE{6,T} end

get_state1_names(state::AbstractStateE)   = [:E]
get_state1_names(state::AbstractStateCN)  = [:C, :N, :E]
get_state1_names(state::AbstractStateCNE) = [:EE, :CN, :C, :N, :E]

init_state(field, statetype) = statetype([0.0unit(field) for n in fieldnames(statetype)]...)
