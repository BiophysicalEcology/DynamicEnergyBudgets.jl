# struct HasE end
# struct HasCN end
# struct HasCNE end

# has_reserves(o) = has_reserves(o.J)
# has_reserves(::AbstractLMatrix{T,A,Syms1,Syms2}) where {T,A,Syms1,Syms2} = 
#     if :E in Syms1
#         if :C in Syms1 && :N in Syms1
#             HasCNE()
#         else
#             HasE()
#         end
#     else
#         HasCN()
#     end
