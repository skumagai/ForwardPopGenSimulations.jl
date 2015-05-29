module SymmetricDominance

export Population, ModelParameters, evolve!, initialize!, reinitialize!, tmrca, nsegsites, distances, spectra

include("genes.jl")

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

function getgenes!(gids::AbstractArray, pop::Population, loci::Int)
    length(gids) == length(pop) * 2 || error("length(gene array) != 2 * population size")
    i = 1
    for org in pop, chr = 1:2
        gids[i] = org[locus, chr]
    end
    nothing
end

function cleandb!(gdb::GeneDB, pop::Population, gids::AbstractArray)
    n = length(pop)
    nl = nloci(pop)
    blocksize = 2 * n
    length(gids) == 2 * n * nl || error("array(# genes) != 2 * population size * number of loci")
    for locus = 1:nloci
        fst, lst = (1 + (locus - 1) * blocksize):(locus * blocksize)
        getgenes!(sub(gids, fst:last), pop, locus)
    end
    clean!(gdb, gids)
    nothing
end

function getparentids!(ps, n)
    rand!(ps, 1:n)
    while ps[1] == ps[2]
        ps[2] = rand(1:n)
    end
    nothing
end

function evolve!(gdb::GeneDB, parpop::Population, params::ModelParameters, state::Int, t::Int, termon::Int, tclean::Int)
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
    lastcoals = [0 for _ = 1:nloci]

    chpop = Population(n, nloci)

    gids = Array(Int, 2 * n)
    allgids = Array(Int, 2 * n * nloci)

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

                    rand() > homofit && continue
                else
                    rand() > heterofit && continue
                end

                for par = 1:2,  locus = 1:nloci
                    if mutarray[locus, par]
                        chpop[i, locus, par] = transmit!(gdb, gen, state, parpop[ps[par], locus, parchrs[par]])
                        state += 1
                    else
                        parent = parpop[ps[par], locus, parchrs[par]]
                        chpop[i, locus, par] = transmit!(gdb, gen, parent, gdb[parent].state)
                    end
                    parchrs[par] = rand() < recombs[locus] ? 3 - parchrs[par] : parchrs[par]
                end
                break
            end
        end
        parpop, chpop = chpop, parpop
        for locus in 1:nloci
            getgenes!(gids, parpop, locus)
            anc = mrca(gdb, gids)
            if anc.epoch > lastcoals[locus]
                lastcoals[locus] = anc.epoch
                ncoals[locus] += 1
                info("turn over ", ncoals[locus], " of locus ", locus, " on generation ", gen)
            end
        end
        if termon == minimum(ncoals)
            break
        end
        gen % tclean == 0 && cleandb!(gdb, parpop, allgids)
    end
    cleandb!(gdb, parpop, allgids)
    parpop, state, gen
end

function initialize!(pop::Population)
    # Initialize a parental population. Genes are distinct.
    gdb = GeneDB()
    state = 1
    for org in pop, locus = 1:nloci(pop), chr = 1:2
        org[locus, chr] = insert!(gdb, GeneRecord(0, state))
        state += 1
    end
    gdb, state
end

function reinitialize!(oldgdb::GeneDB, pop::Population)
    gdb = GeneDB()
    smap = Dict{Int, Int}()
    smax = 1
    for org in pop, locus = 1:nloci(pop), chr = 1:2
        state = oldgdb[org[locus, chr]].state
        if !haskey(smap, state)
            smap[state] = smax
            smax += 1
        end
        org[locus, chr] = insert!(gdb, GeneRecord(0, smap[state]))
    end
    gdb
end

function simulate(params::ModelParameters, burnin::Int, t::Int, termon::Int)

    # This is a parental population, a population of offspring is created within evolve! function.
    pop = Population(params.popsize, params.numberofloci)

    # Burnin
    # Execute the exact-same sequence as main-loop of evolution and throws out lineage information afterwords.
    # This loop runs exacctly "burnin" generations regardless of the presence of coalescence.
    gdb, state = initialize!(pop) # all genes are distinct
    pop, state, t = evolve!(gdb, pop, params, state, burnin, -1)

    # Main loop of evolution
    # This loop terminates upon the first coalescence or after "t" generations.
    gdb = reinitialize!(pop, gdb)
    pop, state, t = evolve!(gdb, pop, params, state, t, termon)
    pop, gdb, t
end

function toarray(gdb::GeneDB, pop::Population, field::Symbol)
    [getfield(gdb[org[locus, chr]], field) for org in pop, locus = 1:nloci(pop), chr = 1:2]
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

end
