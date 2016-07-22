# Sequence-Levenshtein Code
# =========================
#
# Buschmann, Tilo, and Leonid V. Bystrykh. "Levenshtein error-correcting
# barcodes for multiplexed DNA sequencing." BMC bioinformatics 14.1 (2013): 272.

type SequenceLevenshteinCode{k} <: AbstractCode{k}
    codewords::Vector{DNAKmer{k}}
    mindist::Int
end

function distance{k,l}(::Type{SequenceLevenshteinCode{k}}, x::DNAKmer{k}, y::DNAKmer{l})
    return rundp!(Matrix{Int}(k + 1, l + 1), x, y)
end

function distance{k}(::Type{SequenceLevenshteinCode{k}}, kmer::DNAKmer{k}, seq::DNASequence)
    return rundp!(Matrix{Int}(k + 1, length(seq) + 1), kmer, seq)
end

function rundp!(dist, seq1, seq2)
    k = length(seq1)
    l = length(seq2)
    @assert size(dist) == (k + 1, l + 1)

    for i in 0:k
        dist[i+1,1] = i
    end
    for j in 1:l
        dist[1,j+1] = j
        for i in 1:k
            dist[i+1,j+1] = min(
                dist[i,j+1] + 1,
                dist[i+1,j] + 1,
                dist[i,j]   + ifelse(seq1[i] == seq2[j], 0, 1))
        end
    end

    mindist = dist[end,end]
    for i in 0:k
        mindist = min(mindist, dist[i+1,end])
    end
    for j in 0:l
        mindist = min(mindist, dist[end,j+1])
    end

    return mindist
end

function genbarcodes{k}(::Type{SequenceLevenshteinCode{k}}, mindist::Integer, n::Integer=typemax(Int))
    @assert k < 32
    codewords = DNAKmer{k}[]
    x = typemin(UInt64)
    while length(codewords) < n && x < 4^k
        kmer_x = convert(DNAKmer{k}, x)
        for kmer_y in codewords
            if distance(SequenceLevenshteinCode{k}, kmer_x, kmer_y) < mindist
                @goto next
            end
        end
        push!(codewords, kmer_x)
        @label next
        x += 1
    end
    return SequenceLevenshteinCode(codewords, mindist)
end

function demultiplex{k}(code::SequenceLevenshteinCode{k}, seq::DNASequence)
    i_min = 0
    dist_min = typemax(Int)
    for (i, codeword) in enumerate(code.codewords)
        dist = distance(SequenceLevenshteinCode{k}, codeword, seq)
        if dist < dist_min
            i_min = i
            dist_min = dist
        end
    end
    return i_min
end
