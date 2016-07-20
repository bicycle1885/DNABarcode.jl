# Hamming Code
# ============

type HammingCode{k} <: AbstractCode{k}
    codewords::Vector{DNAKmer{k}}
    mindist::Int
end

function distance{k}(::Type{HammingCode{k}}, x::DNAKmer{k}, y::DNAKmer{k})
    return mismatches(x, y)
end

function genbarcodes{k}(::Type{HammingCode{k}}, mindist::Integer, n::Integer=typemax(Int))
    @assert k < 32
    codewords = DNAKmer{k}[]
    x = typemin(UInt64)
    while length(codewords) < n && x < 4^k
        kmer_x = convert(DNAKmer{k}, x)
        for kmer_y in codewords
            if distance(HammingCode{k}, kmer_x, kmer_y) < mindist
                @goto next
            end
        end
        push!(codewords, kmer_x)
        @label next
        x += 1
    end
    return HammingCode(codewords, mindist)
end

function demultiplex{k}(code::HammingCode{k}, seq::DNASequence)
    barcode = extractkmer(DNAKmer{k}, seq, 1)
    i_min = 0
    dist_min = typemax(Int)
    for (i, codeword) in enumerate(code.codewords)
        dist = distance(HammingCode{k}, codeword, barcode)
        if dist < dist_min
            i_min = i
            dist_min = dist
        end
    end
    return i_min
end

function extractkmer{k}(::Type{DNAKmer{k}}, seq::DNASequence, from::Integer)
    x = UInt64(0)
    for i in from:from+k-1
        x <<= 2
        nt = seq[i]
        if isambiguous(nt)
            x |= UInt64(DNA_A)
        else
            x |= UInt64(nt)
        end
    end
    return DNAKmer{k}(x)
end
