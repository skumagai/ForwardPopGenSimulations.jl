export Population

type Population{P<:PopulationParameterSet}
    data::Array{Int, 2}
    popsize::Int
    nloci::Int
    params::P
end

get(p::Population, s::Symbol) = getfield(p.params, s)
get(p::Population, s::Symbol, sex::SexType) = getfield(p.params, s)[value(sex)]
get(p::Population, s::Symbol, i::Integer) = getfield(p.params, s)[i]

function Population{P<:PopulationParameterSet}(popsize::Int, nloci::Int, params::P)
    d = zeros(Int, popsize * nloci, 2)
    Population(d, popsize, nloci, params)
end

@inline populationsize(p::Population) = p.popsize
@inline numberofloci(p::Population) = p.nloci
@inline toindex(p::Population, ind, locus) = numberofloci(p) * (ind - 1) + locus
