module SymmetricDominance

export ModelParameters, simulate, tmrca, nsegsites, distances, spectra

typealias WithoutNewAllele Val{0}
typealias WithNewAllele Val{1}

type GeneDB
    currentid::Int
    currentstate::Int
    id::Vector{Int}
    epoch::Vector{Int}
    parent::Vector{Int}
    state::Vector{Int}

    GeneDB() = new(0, 0, Int[], Int[], Int[], Int[])
end

function insert!(::Type{WithoutNewAllele}, gdb::GeneDB, epoch::Int, parentid::Int)
    gdb.currentid += 1
    push!(gdb.id, gdb.currentid)
    push!(gdb.epoch, epoch)
    push!(gdb.parent, parentid)
    push!(gdb.state, gdb.state[parentid])
    gdb.currentid
end

function insert!(::Type{WithoutNewAllele}, gdb::GeneDB, epoch::Int, parentid::Int, state::Int)
    gdb.currentid += 1
    push!(gdb.id, gdb.currentid)
    push!(gdb.epoch, epoch)
    push!(gdb.parent, parentid)
    push!(gdb.state, state)
    state > gdb.currentstate && (gdb.currentstate = state)
    gdb.currentid
end

function insert!(::Type{WithNewAllele}, gdb::GeneDB, epoch::Int, parentid::Int)
    gdb.currentid += 1
    push!(gdb.id, gdb.currentid)
    push!(gdb.epoch, epoch)
    push!(gdb.parent, parentid)
    gdb.currentstate += 1
    push!(gdb.state, gdb.currentstate)
    gdb.currentid
end

zerof(x) = 0
function select(gdb::GeneDB, geneid::Int, fields::Symbol...)
    if geneid == 0
        ntuple(length(fields), zerof)
    else
        tuple([getfield(gdb, field)[geneid] for field in fields]...)
    end
end

immutable ModelParameters
    popsize::Int
    numberofloci::Int
    heterozygousfitness::Float64
    homozygousfitness::Float64
    recombinationrates::Vector{Float64}
    mutationrates::Vector{Float64}
end

immutable Organism
    geneids::Matrix{Int}
    Organism(nloci::Int) = new(Array(Int, nloci, 2))
end

immutable Population
    data::Vector{Organism}
    Population(n::Int, nloci::Int) = new([Organism(nloci) for _ = 1:n])
end
Base.start(pop::Population) = start(pop.data)
Base.next(pop::Population, iter) = next(pop.data, iter)
Base.done(pop::Population, iter) = done(pop.data, iter)
Base.length(pop::Population) = length(pop.data)

nloci(pop::Population) = size(pop[1].geneids, 1)
nloci(org::Organism) = size(org.geneids, 1)

Base.getindex(org::Organism, locus::Int, chr::Int) = org.geneids[locus, chr]
Base.setindex!(org::Organism, val::Int, locus::Int, chr::Int) = org.geneids[locus, chr] = val
Base.getindex(pop::Population, orgid::Int) = pop.data[orgid]
Base.getindex(pop::Population, orgid::Int, locus::Int, chr::Int) = pop.data[orgid].geneids[locus, chr]
Base.setindex!(pop::Population, val::Int, orgid::Int, locus::Int, chr::Int) = pop[orgid][locus, chr] = val

isidbystate(gdb::GeneDB, g1::Int, g2::Int) = select(gdb, g1, :state) == select(gdb, g2, :state)

mutate!(gdb::GeneDB, epoch::Int, geneid::Int) = insert!(WithNewAllele, gdb, epoch, geneid)

function getparentids!(ps, n)
    rand!(ps, 1:n)
    while ps[1] == ps[2]
        ps[2] = rand(1:n)
    end
    nothing
end

function hascoalesced(gdb::GeneDB, pop::Population, locus::Int, lastcoal::Int)
    coalesced = true
    ancestor = getancestor(gdb, pop[1, locus, 1], lastcoal)
    for org in pop , chr = 1:2
        if ancestor != getancestor(gdb, org[locus, chr], lastcoal)
            coalesced = false
            break
        end
    end
    coalesced
end

function getancestor(gdb::GeneDB, geneid::Int, lastcoal::Int)
    epoch = select(gdb, geneid, :epoch)[1]
    while epoch > lastcoal
        geneid = select(gdb, geneid, :parent)[1]
        epoch = select(gdb, geneid, :epoch)[1]
    end
    geneid
end

function evolve!(gdb::GeneDB, parpop::Population, params::ModelParameters, t::Int, termon::Int)

    # unpacking parameters
    n = params.popsize
    heterofit = params.heterozygousfitness
    homofit = params.homozygousfitness
    # The last element is just there as a placeholder, and its value does not affect runs.
    recombs = [params.recombinationrates; 0.0]
    muts = params.mutationrates
    nloci = params.numberofloci

    # normalize mutation rates
    maxfit = max(heterofit, homofit)
    heterofit /= maxfit
    homofit /= maxfit

    mutarray = Array(Bool, nloci, 2) # boolean value for each gene if it will be mutated.
    ps = Array(Int, 2) # indices of parents of an offspring.
    parchrs = Array(Int, 2) # a gene from which chromosome is passed on to offspring.

    ncoals = [0 for _ = 1:nloci]
    lastcoals = [1 for _ = 1:nloci]

    chpop = Population(n, nloci)

    gen = 1 # current generation
    for gen = 1:t
        for i = 1:n # iterate over offspring
            while true
                getparentids!(ps, n)
                # determine if mutations occur.
                for locus = 1:nloci, chr = 1:2
                    mutarray[locus, chr] = rand() < muts[locus] ? true : false
                end
                # process the first locus, which is under selection. A offspring is homozygous only when it
                # inherits identical-by-state genes from both parents without mutation. Otherwise, the offspring
                # is heterozygous.
                rand!(parchrs, 1:2)
                if isidbystate(gdb, parpop[ps[1], 1, parchrs[1]], parpop[ps[2], 1, parchrs[2]]) &&
                    mutarray[1,1] == mutarray[1,2] == false

                    if rand() > homofit
                        continue
                    end

                else
                    if rand() > heterofit
                        continue
                    end
                end

                for par = 1:2,  locus = 1:nloci
                    if mutarray[locus, par]
                        chpop[i, locus, par] = mutate!(gdb, gen, parpop[ps[par], locus, parchrs[par]])
                    else
                        chpop[i, locus, par] = parpop[ps[par], locus, parchrs[par]]
                    end
                    parchrs[par] = rand() < recombs[locus] ? 3 - parchrs[par] : parchrs[par]
                end
                break
            end
        end
        parpop, chpop = chpop, parpop
        for locus in 1:nloci
            if hascoalesced(gdb, parpop, locus, lastcoals[locus])
                lastcoals[locus] = gen
                ncoals[locus] += 1
                reset!(gdb, parpop, gen, locus)
                write(STDERR, string("INFO: turn over ", ncoals[locus], " of locus ", locus, " on generation ", gen, "\n"))
            end
        end
        if termon == minimum(ncoals)
            break
        end
    end
    parpop, gen
end

function initialize!(pop::Population)
    # Initialize a parental population. Genes are distinct.
    gdb = GeneDB()
    for org in pop, locus = 1:nloci(pop), chr = 1:2
        org[locus, chr] = insert!(WithNewAllele, gdb, 1, 0)
    end
    gdb
end

function reset!(gdb::GeneDB, pop::Population, gen::Int, locus::Int)
    # Reset lineage info of genes so that all genes can be distinguished even for identical allelic state.
    for org in pop, chr = 1:2
        org[locus, chr] = insert!(WithoutNewAllele, gdb, gen, org[locus, chr])
    end
    nothing
end

function reinitialize!(pop::Population, oldgdb::GeneDB)
    gdb = GeneDB()
    for org in pop, locus = 1:nloci(pop), chr = 1:2
        org[locus, chr] = insert!(WithoutNewAllele, gdb, 1, 0, select(oldgdb, org[locus, chr], :state)[1])
    end
    gdb
end

function simulate(params::ModelParameters, burnin::Int, t::Int, turmon::Int)

    # This is a parental population, a population of offspring is created within evolve! function.
    pop = Population(params.popsize, params.numberofloci)

    # Burnin
    # Execute the exact-same sequence as main-loop of evolution and throws out lineage information afterwords.
    # This loop runs exacctly "burnin" generations regardless of the presence of coalescence.
    gdb = initialize!(pop) # all genes are distinct
    pop, t = evolve!(gdb, pop, params, burnin, -1)

    # Main loop of evolution
    # This loop terminates upon the first coalescence or after "t" generations.
    gdb = reinitialize!(pop, gdb)
    pop, t = evolve!(gdb, pop, params, t, turmon)
    pop, gdb, t
end

function toarray(gdb::GeneDB, pop::Population, field::Symbol)
    Int[select(gdb, org[locus, chr], field)[1] for org in pop, locus = 1:nloci(pop), chr = 1:2]
end

function counts(gdb::GeneDB, pop::Population)
    alleles = toarray(gdb, pop, :state)
    nl = nloci(pop)
    # allele count
    adata = [Dict{Int, Int}() for _ = 1:nl]
    # genotype count
    gdata = [Dict{NTuple{2, Int}, Int}() for _ = 1:nl]
    # haplotype count
    hdata = Dict{NTuple{nl, Int}, Int}()
    for org in 1:length(pop)
        for locus = 1:nl
            g = tuple(sort(vec(alleles[org, locus, :]))...)
            gdata[locus][g] = get(gdata[locus], g, 0) + 1
            for chr = 1:2
                a = alleles[org, locus, chr]
                adata[locus][a] = get(adata[locus], a, 0) + 1
            end
        end
        for chr = 1:2
            h = tuple(vec(alleles[org, :, chr])...)
            hdata[h] = get(hdata, h, 0) + 1
        end
    end
    adata, gdata, hdata
end

function spectra(gdb::GeneDB, pop::Population)
    adata, gdata, hdata = counts(gdb, pop)

    nl = nloci(pop)
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

function history(gdb::GeneDB, idx::Int, lastcoal::Int)
    val = [idx]
    epoch = select(gdb, idx, :epoch)[1]
    while epoch > lastcoal
        idx = select(gdb, idx, :parent)[1]
        epoch = select(gdb, idx, :epoch)[1]
        push!(val, idx)
    end
    reverse(val)
end

function distances(gdb::GeneDB, pop::Population, lastcoal::Int)
    alleles = toarray(gdb, pop, :id)
    states = Dict{Int, Int}()
    for org in pop, locus = 1:nloci(pop), chr = 1:2
        oid = org[locus, chr]
        states[oid] = select(gdb, oid, :state)[1]
    end

    dists = Array(Int, 4, 0)

    idx = 1

    for locus = 1:nloci(pop)
        lineages = sort(unique(alleles[:,locus,:]))
        dists = hcat(dists, Array(Int, 4, binomial(length(lineages), 2)))
        nl = length(lineages)
        h = Dict{Int, Vector{Int}}()
        for i = 1:nl
            l = lineages[i]
            h[l] = history(gdb, l, lastcoal)
        end
        for i = 1:(nl-1), j = (i+1):nl
            h1, h2 = h[lineages[i]], h[lineages[j]]
            ca = intersect(h1, h2)
            dists[1, idx] = locus
            dists[2, idx] = states[lineages[i]]
            dists[3, idx] = states[lineages[j]]
            if length(ca) > 0
                locca1 = findfirst(h1, ca[end])
                locca2 = findfirst(h2, ca[end])
                dists[4, idx] = length(h1) - locca1 + length(h2) - locca2
            else
                dists[4, idx] = -1
            end
            idx += 1
        end
    end
    dists'
end

function nsegsites(gdb::GeneDB, pop::Population, lastcoal::Int)
    # The number of segregating sites in the number of mutations in a tree up to the most recent common ancestor.
    lineages = toarray(gdb, pop, :id)

    nl = nloci(pop)
    nms = Array(Int, nl)
    for locus = 1:nl
        # First, check if every loci have MRCA. If not, number of mutations is not well-defined. In this case, simply return -1.
        if !hascoalesced(gdb, pop, locus, lastcoal)
            nms[locus] = -1
            continue
        end

        # Then, construct a set of mutations along a tree up until MRCA. This will include
        # MRCA. However, this should not be in the final number, as all samples have the same allele.
        ids = unique(lineages[:, locus, :])
        hists = Vector{Int}[history(gdb, id, lastcoal) for id in ids]
        commons = copy(hists[1])
        for h in hists
            commons = intersect(commons, h)
        end
        ca = commons[end]
        muts = Int[]
        for h in hists
            muts = union(muts, h[findfirst(h, ca):end])
        end
        nms[locus] = length(muts) - 1
    end
    nms
end

function tmrca(gdb::GeneDB, pop::Population)
    # Computes time until the most recent common ancestor.
    nl= nloci(pop)
    ts = Array(Int, nl)

    for locus = 1:nl
        path = IntSet(history(gdb, pop[1, locus, 1], 0))
        for org in pop, chr = 1:2
            intersect!(path, IntSet(history(gdb, org[locus, chr], 0)))
        end
        ts[locus] = select(gdb, maximum(path), :epoch)[1]
    end
    ts
end

end
