using DNABarcode
using Base.Test
using Bio.Seq

function hamming_distance(x, y)
    @assert length(x) == length(y)
    dist = 0
    for i in 1:endof(x)
        dist += x[i] != y[i]
    end
    return dist
end

function randdna()
    r = rand()
    if r < 0.20
        return DNA_A
    elseif r < 0.40
        return DNA_C
    elseif r < 0.60
        return DNA_G
    elseif r < 0.80
        return DNA_T
    else
        return DNA_N
    end
end

@testset "HammingCode" begin
    for k in 2:6, mindist in 2:k
        barcodes = genbarcodes(HammingCode{k}, mindist)
        @test !isempty(barcodes)
        ok = true
        for i in 1:length(barcodes)
            for j in i+1:length(barcodes)
                if hamming_distance(barcodes[i], barcodes[j]) < mindist
                    ok = false
                end
            end
        end
        @test ok

        for i in rand(1:length(barcodes), 10)
            barcode = barcodes[i]
            seq = DNASequence(barcode) * dna"ACGTAT"
            @test demultiplex(barcodes, seq) == i
            for _ in 1:fld(mindist - 1, 2)
                seq[rand(1:k)] = randdna()
            end
            @test demultiplex(barcodes, seq) == i
        end
    end
end
