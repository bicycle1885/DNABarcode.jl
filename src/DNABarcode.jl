module DNABarcode

export
    HammingCode,
    genbarcodes,
    demultiplex

import Bio.Seq: DNAKmer, DNASequence, DNA_A, mismatches, isambiguous

include("abstractcode.jl")
include("hamming.jl")

end # module
