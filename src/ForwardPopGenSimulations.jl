module ForwardPopGenSimulations

export toarray,
       history,
       ca,
       mrca,
       nsegsites,
       distances,
       isidbystate,
       selectparents!,
       selectmutatedsites!

include("chromosomes.jl")
include("sexes.jl")
include("genes.jl")
include("basicdata.jl")
include("populations.jl")
include("postsims.jl")

isidbystate(gdb::GeneDB, id1::Int, id2::Int) = gdb.data[id1].state == gdb.data[id2].state

function hascoalesced(gdb::GeneDB, locus::Int, pops::Tuple{Population, SexType})
    coalesced = true
    ancestor = 0
    for (pop, _) in pops, org = 1:length(pop), chr = 1:ploidy(pop)
        anctmp = gdb[pop[org, locus, chr]].ancestor
        if anctmp != 0
            ancestor = anctmp
            break
        end
    end
    for (pop, _) in pops, org = 1:length(pop), chr = 1:ploidy(pop)
        anctmp = gdb[pop[org, locus, chr]].ancestor
        anctmp == 0 && continue
        if ancestor != anctmp
            coalesced = false
            break
        end
    end
    coalesced
end

function selectparents!(ps, n; replace=false)
    rand!(ps, 1:n)
    lps = length(ps)
    if replace == false && lps > 1
        for i = 2:lps
            while in(ps[i], sub(ps, 1:(i-1)))
                ps[i] = rand(1:n)
            end
        end
    end
    nothing
end

function selectmutatedsites!(mutarray, rates)
    size(mutarray) == size(rates) || error("Dimension mismatch in a mutation array and rate array.")
    for iter = eachindex(mutarray)
        mutarray[iter] = rand() < rates[iter]
    end
    nothing
end

end
