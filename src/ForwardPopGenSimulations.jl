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

isidbystate(gdb::GeneDB, id1::Int, id2::Int) = gdb.data[id1].state == gdb.data[id2].state

function hascoalesced(gdb::GeneDB, gids::AbstractArray)
    coalesced = true
    ancestor = gdb[first(gids)].ancestor
    for gid in gids
        if ancestor != gdb[gids].ancestor
            coalesced = false
            break
        end
    end
    coalesced
end

function history(gdb::GeneDB, idx::Int)
    val = [idx]
    while gdb[idx] != gdb[idx].parent
        idx = gdb[idx].parent.id
        push!(val, idx)
    end
    reverse(val)
end

_getca(gdb::GeneDB, path) = length(path) > 0 ? gdb[maximum(path)] : UndefGene

function ca(gdb::GeneDB, id1::Int, id2::Int)
    path = Set(history(gdb, id1))
    path = intersect(path, Set(history(gdb, id2)))
    _getca(gdb, path)
end

function mrca(gdb::GeneDB, gids::AbstractArray)
    path = Set(history(gdb, first(gids)))
    for gid in gids
        path = intersect(path, Set(history(gdb, gid)))
    end
    _getca(gdb, path)
end

toarray(gdb::GeneDB, gids::AbstractArray, field::Symbol) = [getfield(gdb[gid], field) for gid in gids]

function distances(gdb::GeneDB, gids::AbstractArray)
    ids = toarray(gdb, gids, :id)

    ngenes = length(gids)
    dists = zeros(Int, ngenes, ngenes)

    idx = 1

    tnow = gdb[first(gids)].epoch

    for i = 1:(ngenes - 1), j = (i + 1):ngenes
        dists[i, j] = dists[j, i] = tnow - ca(gdb, ids[i], ids[j]).epoch
    end
    dists
end

function nsegsites(gdb::GeneDB, gids::AbstractArray)
    # The number of segregating sites in the number of mutations in a tree up to the most recent common ancestor.

    states = Set([gdb[h].state for h in history(gdb, first(gids))])
    for gid in gids
        union!(states, Set([gdb[h].state for h in history(gdb, gid)]))
    end

    length(states) - 1
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
