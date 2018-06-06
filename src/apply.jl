export apply, offset_apply

function apply(f::F, a::Tuple{A,Vararg}, args...) where {F,A}
    f(a[1], args...)
    apply(f, Base.tail(a), args...)
end
function apply(f::F, a::Tuple{A,Vararg}, b::Tuple{B,Vararg}, args...) where {F,A,B}
    f(a[1], b[1], args...)
    apply(f, Base.tail(a), Base.tail(b), args...)
end
function apply(f::F, a::Array{T}, b::Tuple{B,Vararg}, args...) where {F,T,B}
    f(a[1], b[1], args...)
    apply(f, a[2:end], Base.tail(b), args...)
end
function apply(f::F, a::Tuple{A,Vararg}, b::Array{T}, args...) where {F,T,A}
    f(a[1], b[1], args...)
    apply(f, Base.tail(a), b[2:end], args...)
end
function apply(f::F, a::Array{T}, b::Array{T}, args...) where {F,T}
    f(a[1], b[1], args...)
    apply(f, a[2:end], b[2:end], args...)
end
apply(f, a::Tuple{}, b::Tuple{}, args...)::Void = nothing
apply(f, a, b::Tuple{}, args...)::Void = nothing
apply(f, a::Tuple{}, b, args...)::Void = nothing
apply(f, a::Tuple{}, args...)::Void = nothing

offset_apply!(f, a::AbstractArray, o::Tuple, args...) = offset_apply!(f, a, 0, o, args...)
offset_apply!(f, a::AbstractArray, offset::Int, o::Tuple{O,Vararg}, args...) where O = begin
    offset = f(a, offset, o[1], args...)
    offset_apply!(f, a, offset::Int, Base.tail(o), args...)
end
offset_apply!(f, a::AbstractArray, offset::Int, o::Tuple{}, args...) = nothing
