

"Weighted sampler that gets sampled based on GC-ratio"
function dna_sampler(GC_ratio::Float64)::SamplerWeighted
    SamplerWeighted(
        dna"ACGTN",
        [1-GC_ratio, GC_ratio, GC_ratio, 1-GC_ratio] .* 0.5,
    )
end


function random_dna(len::Integer, sampler::SamplerWeighted)::LongDNA{4}
    LongDNA{4}(rand(sampler, len))
end

function random_dna(len::Integer, GC_ratio::Float64 = 0.5)::LongDNA{4}
    random_dna(len, dna_sampler(GC_ratio))
end

export random_dna