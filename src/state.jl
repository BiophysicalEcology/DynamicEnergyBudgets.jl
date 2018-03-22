abstract type AbstractState{N} <: FieldVector{N,Float64} end
abstract type AbstractStateE{N} <: AbstractState{N} end
abstract type AbstractStateCN{N} <: AbstractState{N} end
abstract type AbstractStateCNE{N} <: AbstractState{N} end

@traitdef StateHasM{X}
@traitdef StateHasP{X}
SimpleTraits.trait{X1}(::Type{StateHasM{X1}}) = :M in fieldnames(X1) ? StateHasM{X1} : Not{StateHasM{X1}}
SimpleTraits.trait{X1}(::Type{StateHasP{X1}}) = :P in fieldnames(X1) ? StateHasP{X1} : Not{StateHasP{X1}}

mutable struct StateV <: AbstractState{2}
    V::Float64
end
mutable struct StatePV <: AbstractState{2}
    P::Float64
    V::Float64
end
mutable struct StateVE <: AbstractStateE{2}
    V::Float64
    E::Float64
end
mutable struct StatePVE <: AbstractStateCN{3}
    P::Float64
    V::Float64
    E::Float64
end
mutable struct StateVCN <: AbstractStateCN{3}
    V::Float64
    C::Float64
    N::Float64
end
mutable struct StateVME <: AbstractStateE{3}
    V::Float64
    M::Float64
    E::Float64
end
mutable struct StatePVME <: AbstractStateE{4}
    P::Float64
    V::Float64
    M::Float64
    E::Float64
end
mutable struct StatePVCN <: AbstractStateCN{4}
    P::Float64
    V::Float64
    C::Float64
    N::Float64
end
mutable struct StatePVMCN <: AbstractStateCN{5}
    P::Float64
    V::Float64
    M::Float64
    C::Float64
    N::Float64
end
mutable struct StatePVCNE <: AbstractStateCNE{5}
    P::Float64
    V::Float64
    C::Float64
    N::Float64
    E::Float64
end
mutable struct StatePVMCNE <: AbstractStateCNE{6}
    P::Float64
    V::Float64
    M::Float64
    C::Float64
    N::Float64
    E::Float64
end

get_state1_names(state::AbstractStateE)   = [:E]
get_state1_names(state::AbstractStateCN)  = [:C, :N, :E]
get_state1_names(state::AbstractStateCNE) = [:EE, :CN, :C, :N, :E]

function init_state(state_type)
    state_type(zeros(length(fieldnames(state_type)))...)
end
