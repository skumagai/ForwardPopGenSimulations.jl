export Migration, GeneDB, GeneRecord, UndefGene
export transmit!, toarray, history, ca, mrca, nsegsites, distances, nextid!

immutable Transmission end
immutable Mutation end
immutable Migration
    from::Int32
    to::Int32
end

type GeneRecord
    # Primary information of lineages: Immutable throughout a simulation.
    epoch::Int
    state::Int
    event::Union(Transmission, Mutation, Migration)
    # id = 0 is a special type that the gene does not appear explicitly in a simulation.
    # In other words, id = 0 indicates the identity of a gene is unknown, so I can't say
    # gene(id = 0) == gene(id = 0).
    id::Int

    # Secondary information of lineages: This will be updated during a simulation.
    parent::GeneRecord

    # id is intentionally left unspecified. That field is specified upon insertion in to GeneDB.
    function GeneRecord(epoch::Int, state::Int)
        self = new(epoch, state, Transmission())
        self.parent = self
        self
    end
    function GeneRecord(epoch::Int, parent::GeneRecord)
        self = new(epoch, parent.state, Transmission())
        self.parent = parent
        self
    end
    function GeneRecord(epoch::Int, state::Int, parent::GeneRecord)
        self = new(epoch, state, Mutation())
        self.parent = parent
        self
    end
    function GeneRecord(epoch::Int, state::Int, m::Migration, parent::GeneRecord)
        self = new(epoch, state, m)
        self.parent = parent
        self
    end
end

const UndefGene = (ug = GeneRecord(0, 0); ug.id = 0; ug)

type GeneDB
    currentid::Int
    data::Dict{Int, GeneRecord}

    GeneDB() = new(0, Dict{Int, GeneRecord}())
end

Base.getindex(gdb::GeneDB, id::Int) = gdb.data[id]
Base.setindex!(gdb::GeneDB, record::GeneRecord, id::Int) = gdb.data[id] = record
Base.haskey(gdb::GeneDB, id::Int) = haskey(gdb.data, id)
Base.delete!(gdb::GeneDB, id::Int) = delete!(gdb.data, id)
Base.keys(gdb::GeneDB) = keys(gdb.data)

update!(gdb::GeneDB, id::Int, field::Symbol, value) = setfield!(gdb[id], field, value)
update!(gdb::GeneDB, record::GeneRecord, field::Symbol, value) = setfield!(record, field, value)

drop!(gdb::GeneDB, id::Int) = delete!(gdb.data, id)
drop!(gdb::GeneDB, record::GeneRecord) = delete!(gdb.data, record.id)

nextid!(gdb::GeneDB) = (gdb.currentid += 1; gdb.currentid)

function insert!(gdb::GeneDB, record::GeneRecord)
    id = nextid!(gdb)
    record.id = id
    haskey(gdb, id) && error("Attempted to insert an existing record.")
    gdb[id] = record
    id
end

isidbystate(gdb::GeneDB, id1::Int, id2::Int) = gdb.data[id1].state == gdb.data[id2].state

# Without mutation
transmit!(gdb::GeneDB, epoch::Int, pid::Int) = insert!(gdb, GeneRecord(epoch, gdb[pid]))
# With mutation
transmit!(gdb::GeneDB, epoch::Int, state::Int, pid::Int) = insert!(gdb, GeneRecord(epoch, state, gdb[pid]))

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

function _clean!(gdb::GeneDB, gids::AbstractArray)
    cs = Dict{Int, Array{Int, 1}}()
    # Organisms in yonger generations are guaranteed to have higher ID.
    # Reverse iteration visits offspring before parents. This makes it possible to construct a familial tree
    # where parents know offspring.
    ids = collect(keys(gdb))
    sort!(ids, rev=true)
    for id in ids
        pid = gdb[id].parent.id
        id == pid && continue
        if haskey(cs, pid)
            push!(cs[pid], id)
        else
            cs[pid] = [id]
        end
    end

    # After the familial tree is built, identify organisms with zoero or one descendants. In case of
    # no descendant, an organism is simply removed from the database regardless of its event type (
    # Transmission, Mutation, Migration). On the other hand when there is only one offspring, a organism
    # gets removed only if its event type is Transmission.
    for id in ids
        p = gdb[id].parent
        pid = p.id
        if !haskey(cs, id) && !in(id, gids)
            # This organism neither has offspring nor extant.
            id != pid && deleteat!(cs[pid], findfirst(id))
            drop!(gdb, id)
        elseif haskey(cs, id) && length(cs[id]) == 1 && isa(gdb[id].event, Transmission)
            c = gdb[cs[id][1]]
            if id == pid
                c.parent = c
            else
                c.parent = p
                cs[pid][findfirst(cs[pid], id)] = cs[id][1]
            end
            delete!(cs, id)
            drop!(gdb, id)
        end
    end
end

function _trimancestors!{T}(gdb::GeneDB, gids::AbstractArray{T, 1})
    # Up to this point, Mutation and Migration nodes ancestoral to MRCA remain in the database.
    # As the last step, remove those nodes.
    ca = mrca(gdb, gids)
    if ca != ca.parent
        anc = ca.parent
        ca.parent = ca
        while true
            old = anc
            anc = anc.parent
            drop!(gdb, old)
            old == anc && break
        end
    end
    nothing
end

function clean!{T}(gdb::GeneDB, gids::AbstractArray{T, 1})
    _clean!(gdb, gids)
    _trimancestors!(gdb, gids)
    nothing
end

function clean!{T}(gdb::GeneDB, gids::AbstractArray{T, 2})
    _clean!(gdb, gids)
    for locus = 1:size(gids, 2)
        _trimancestors(gdb, sub(gids, :, locus))
    end
    nothing
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
