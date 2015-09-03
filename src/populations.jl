export Population,
       initialize!,
       reinitialize!,
       nloci,
       ploidy,
       offset,
       listgenes!


immutable Population
    data::Vector{Int}
    nloci::Int
    ploidy::Int
    popsize::Vector{Int}
    offset::Vector{Int}
    Population(n::Int, nloci::Int, ploidy::Int) =
        new(Vector{Int}(n * nloci * ploidy), nloci, ploidy, Int[n], Int[0])
    function Population(n::Vector{Int}, nloci::Int, ploidy::Int)
        offset = cumsum([0; n])[1:end-1]
        factor = nloci * ploidy
        for i in eachindex(offset)
            offset[i] *= factor
        end
        new(Vector{Int}(sum(n) * nloci * ploidy), nloci, ploidy, n, offset)
    end
end

Base.start(pop::Population) = start(pop.data)
Base.next(pop::Population, iter) = next(pop.data, iter)
Base.done(pop::Population, iter) = done(pop.data, iter)
Base.length(pop::Population) = sum(pop.popsize)
Base.size(pop::Population) = pop.popsize
Base.eachindex(pop::Population) = 1:length(pop.data)

@inline nloci(pop::Population) = pop.nloci
@inline ploidy(pop::Population) = pop.ploidy
@inline offset(pop::Population) = pop.offset
@inline offset(pop::Population, deme) = pop.offset[deme]

function toindex(pop::Population, orgid, locus, chr, deme)
    nl = nloci(pop)
    pl = ploidy(pop)
    # First, get an index local to a deme.
    pos = nl * pl * (orgid - 1) + pl * (locus - 1) + chr
    # Now converting to global index.
    pos + offset(pop, deme)
end

Base.getindex(pop::Population, i::Int) = pop.data[i]
function Base.getindex(pop::Population, orgid::Int, locus::Int, chr::Int, deme::Int=1)
    pop.data[toindex(pop, orgid, locus, chr, deme)]
end

Base.setindex!(pop::Population, val::Int, i::Int) = pop.data[i] = val
function Base.setindex!(pop::Population, val::Int, orgid::Int, locus::Int, chr::Int, deme::Int=1)
    pop.data[toindex(pop, orgid, locus, chr, deme)] = val
end

# default version to initialize populations.
function initialize!(core, chrtypes::Vector{ChromosomeType}, pops::Tuple{Population, SexType}...)
    for (pop, sex) in pops
        pp = length(pop)
        nl = nloci(pop)
        pl = ploidy(pop)
        for org = 1:pp, loc = 1:nl, chr = 1:pl
            pop[org, loc, chr] = initgenotype!(core, Val{sex}, Val{chr}, Val{chrtypes[loc]})
        end
    end
    nothing
end

# When simulating a population of asexual diploids, all organisms are females. Moreover,
# all simulated loci are autosomal.
initialize!(core, pop::Population) = initialize!(core, fill(Autosome, nloci(pop)), (pop, Female))

# default case
initgenotype!(core, locus, chr, ctype) = 0
# specific case; all cases perform identical.s
initgenotype!(core, sex, chr, ::Type{Val{Autosome}}) =
    insert!(db(core), GeneRecord(0, nextstate!(core)))
initgenotype!(core, ::Type{Val{Female}}, chr, ::Type{Val{XChromosome}}) =
    insert!(db(core), GeneRecord(0, nextstate!(core)))
initgenotype!(core, ::Type{Val{Male}}, ::Type{Val{1}}, ::Type{Val{XChromosome}}) =
    insert!(db(core), GeneRecord(0, nextstate!(core)))
initgenotype!(core, ::Type{Val{Male}}, ::Type{Val{2}}, ::Type{Val{YChromosome}}) =
    insert!(db(core), GeneRecord(0, nextstate!(core)))
initgenotype!(core, ::Type{Val{Female}}, ::Type{Val{1}}, ::Type{Val{Mitochondrion}}) =
    insert!(db(core), GeneRecord(0, nextstate!(core)))

function reinitialize!(oldcore, pops::Tuple{Population, SexType}...)
    oldgdb = db(oldcore)
    core = BasicData()
    gdb = db(core)
    # Preemptively map 0 to 0 as this value has a special meaning.
    smap = Dict{Int, Int}(0 => 0)
    smax = 1
    for (pop, _) in pops
        for i in eachindex(pop)
            state = oldgdb[pop[i]].state
            if !haskey(smap, state)
                smap[state] = smax
                smax += 1
            end
            pop[i] = insert!(gdb, GeneRecord(0, smap[state]))
        end
    end
    core.state = maximum(keys(db(core)))
    core
end

reinitialize!(oldcore, pop::Population) = reinitialize!(oldcore, (pop, Female))

function listgenes!(gids::AbstractArray, locus, pops::Tuple{Population, SexType}...)
    gi = 1
    for (pop, _) in pops
        popsize = sum(size(pop))
        nchr = ploidy(pop)
        for org in 1:popsize, chr = 1:nchr
            geneid = pop[org, locus, chr]
            if geneid != 0
                gids[gi] = geneid
                gi += 1
            end
        end
    end
end

listgenes!(gids::AbstractArray, locus, pop::Population) = listgenes!(gids, locus, (pop, Female))
