export apply

function apply(f::F, a::Tuple{A,Vararg}, args...) where {F,A}
    f(a[1], args...)
    apply(f, Base.tail(a), args...)
end
function apply(f::F, a::Tuple{A,Vararg}, b::Tuple{B,Vararg}, args...) where {F,A,B}
    f(a[1], b[1], args...)
    apply(f, Base.tail(a), Base.tail(b), args...)
end
function apply(f::F, a::AbstractArray{T}, b::Tuple{B,Vararg}, args...) where {F,T,B}
    f(a[1], b[1], args...)
    apply(f, a[2:end], Base.tail(b), args...)
end
function apply(f::F, a::Tuple{A,Vararg}, b::AbstractArray{T}, args...) where {F,T,A}
    f(a[1], b[1], args...)
    apply(f, Base.tail(a), b[2:end], args...)
end
function apply(f::F, a::AbstractArray{T}, b::AbstractArray{T}, args...) where {F,T}
    f(a[1], b[1], args...)
    apply(f, a[2:end], b[2:end], args...)
end
apply(f, a::Tuple{}, b::Tuple{}, args...) = nothing
apply(f, a, b::Tuple{}, args...) = nothing
apply(f, a::Tuple{}, args...) = nothing

offset_apply!(f, o::Tuple{O,Vararg}, a::AbstractArray, offset::Int, args...) where O = begin
    offset = f(o[1], a, offset, args...)
    offset_apply!(f, Base.tail(o), a, offset::Int, args...)
end
offset_apply!(f, a::AbstractArray, o::Tuple{O,Vararg}, offset::Int, args...) where O = begin
    offset = f(a, o[1], offset, args...)
    offset_apply!(f, a, Base.tail(o), offset::Int, args...)
end
offset_apply!(f, a::AbstractArray, o::Tuple{}, offset::Int, args...) = nothing
offset_apply!(f, o::Tuple{}, a::AbstractArray, offset::Int, args...) = nothing
