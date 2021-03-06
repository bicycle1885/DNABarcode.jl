# Abstract Code
# =============

# The existence of .codeword::Vector{DNAKmer{k}} is assumed.
abstract AbstractCode{k} <: AbstractVector{DNAKmer{k}}

function Base.getindex(code::AbstractCode, i::Integer)
    return code.codewords[i]
end

function Base.size(code::AbstractCode)
    return size(code.codewords)
end

function Base.filter{k}(f::Function, code::AbstractCode{k})
    return filter!(f, deepcopy(code))
end

function Base.filter!{k}(f::Function, code::AbstractCode{k})
    filter!(f, code.codewords)
    return code
end
