# Functions to analyze population data after simulations finish.
export toarray,
       counts,
       spectra

function toarray(gdb::GeneDB, field::Symbol, pops::Population...)
    initpop = pops[1]
    nl = nloci(initpop)
    pl = ploidy(initpop)
    totalnl = 0
    for pop in pops
        nlpop = nloci(pop)
        nlpop != nl || ploidy(pop) != pl && error("Incompatible populations.")
        totalnl += nlpop
    end

    hcat([getfield(gdb[pop[org, locus, chr]], field)
          for pop in pops,
              org in length(ppo),
              locus = 1:nloci(pop),
              chr = 1:ploidy(pop)])
end

function counts(gdb::GeneDB, pops::Population...)
    alleles = toarray(gdb, :state, pops...)
    # if control gets to this line, it means all population have the same number of loci and
    # ploidy. And, no more tests are needed.
    poptotal = 0
    for pop in pops
        poptotal += length(pop)
    end
    nl = nloci(pops[1][1])
    pl = ploidy(pops[1][1])
    # allele count
    adata = [Dict{Int, Int}() for _ = 1:nl]
    # genotype count
    gdata = [Dict{NTuple{pl, Int}, Int}() for _ = 1:nl]
    # haplotype count
    hdata = Dict{NTuple{nltotal, Int}, Int}()
    for org = 1:poptotal
        for locus = 1:nl
            g = tuple(sort(vec(alleles[org, locus, :]))...)
            gdata[locus][g] = get(gdata[locus], g, 0) + 1
            for chr = 1:2
                a = alleles[org, locus, chr]
                adata[locus][a] = get(adata[locus], a, 0) + 1
            end
        end
        for chr = 1:pl
            h = tuple(vec(alleles[org, :, chr])...)
            hdata[h] = get(hdata, h, 0) + 1
        end
    end
    adata, gdata, hdata
end

function spectra(gdb::GeneDB, pops::Population...)
    adata, gdata, hdata = counts(gdb, pops)

    nl = nloci(pops[1][1])
    afs = [Dict{Int, Int}() for _ = 1:nl]
    gfs = [Dict{Int, Int}() for _ = 1:nl]
    hfs = Dict{Int, Int}()
    for locus = 1:nl
        for v in values(adata[locus])
            afs[locus][v] = get(afs[locus], v, 0) + 1
        end
    end
    for locus = 1:nl
        for v in values(gdata[locus])
            gfs[locus][v] = get(gfs[locus], v, 0) + 1
        end
    end
    for v in values(hdata)
        hfs[v] = get(hfs, v, 0) + 1
    end
    afs, gfs, hfs
end


function distances(gdb::GeneDB, locus, pops::Population...)
    ids = toarray(gdb, :id, pops)

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

function history(gdb::GeneDB, idx::Int)
    val = [idx]
    while gdb[idx].parent !== UndefGene
        idx = gdb[idx].parent.id
        push!(val, idx)
    end
    reverse(val)
end

