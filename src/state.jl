abstract type AbstractState{N,T} <: FieldVector{N,T} end
abstract type AbstractStateE{N,T} <: AbstractState{N,T} end
abstract type AbstractStateCN{N,T} <: AbstractState{N,T} end
abstract type AbstractStateCNE{N,T} <: AbstractState{N,T} end

@traitdef StateHasM{X}
@traitdef StateHasP{X}
SimpleTraits.trait{X1}(::Type{StateHasM{X1}}) = :M in fieldnames(X1) ? StateHasM{X1} : Not{StateHasM{X1}}
SimpleTraits.trait{X1}(::Type{StateHasP{X1}}) = :P in fieldnames(X1) ? StateHasP{X1} : Not{StateHasP{X1}}

@def P begin P::T = 0.0 end
@def V begin V::T = 1e-4 end
@def M begin M::T = 0.0 end
@def C begin C::T = 1e-4 end
@def N begin N::T = 1e-4 end
@def E begin E::T = 10.0 end

@with_kw mutable struct StateV{T} <: AbstractState{1,T}
    @V
end
@with_kw mutable struct StateE{T} <: AbstractStateE{1,T}
    @E
end
@with_kw mutable struct StatePV{T} <: AbstractState{2,T}
    @P
    @V
end
@with_kw mutable struct StateVE{T} <: AbstractStateE{2,T}
    @V
    @E
end
@with_kw mutable struct StateCN{T} <: AbstractStateCN{2,T}
    @C
    @N
end
@with_kw mutable struct StatePVE{T} <: AbstractStateCN{3,T}
    @P
    @V
    @E
end
@with_kw mutable struct StateVCN{T} <: AbstractStateCN{3,T}
    @V
    @C
    @N
end
@with_kw mutable struct StateVME{T} <: AbstractStateE{3,T}
    @V
    @M
    @E
end
@with_kw mutable struct StateCNE{T} <: AbstractStateCNE{3,T}
    @C
    @N
    @E
end
@with_kw mutable struct StatePVME{T} <: AbstractStateE{4,T}
    @P
    @V
    @M
    @E
end
@with_kw mutable struct StatePVCN{T} <: AbstractStateCN{4,T}
    @P
    @V
    @C
    @N
end
@with_kw mutable struct StatePVMCN{T} <: AbstractStateCN{5,T}
    @P
    @V
    @M
    @C
    @N
end
@with_kw mutable struct StatePVCNE{T} <: AbstractStateCNE{5,T}
    @P
    @V
    @C
    @N
    @E
end
@with_kw mutable struct StatePVMCNE{T} <: AbstractStateCNE{6,T}
    @P
    @V
    @M
    @C
    @N
    @E
end

get_state1_names(state::AbstractStateE)   = [:E]
get_state1_names(state::AbstractStateCN)  = [:C, :N, :E]
get_state1_names(state::AbstractStateCNE) = [:EE, :CN, :C, :N, :E]

init_state(field, statetype) =
    statetype([0.0unit(field) for n in fieldnames(statetype)]...)
