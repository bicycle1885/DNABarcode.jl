# Abstract Code
# =============

abstract AbstractCode{k} <: AbstractVector{DNAKmer{k}}

function Base.getindex(code::AbstractCode, i::Integer)
    return code.codewords[i]
end

function Base.size(code::AbstractCode)
    return size(code.codewords)
end
