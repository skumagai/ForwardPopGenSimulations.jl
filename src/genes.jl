type GeneRecord
    # Primary information of lineages: Immutable throughout a simulation.
    epoch::Int
    state::Int
    locus::Int
    # id = 0 is a special type that the gene does not appear explicitly in a simulation.
    # In other words, id = 0 indicates the identity of a gene is unknown, so I can't say
    # gene(id = 0) == gene(id = 0).
    id::Int

    # Secondary information of lineages: This will be updated during a simulation.
    parent::GeneRecord

    # id is intentionally left unspecified. That field is specified upon insertion in to GeneDB.
    function GeneRecord(epoch::Int, state::Int, locus::Int)
        self = new(epoch, state, locus)
        self.parent = self
        self
    end
    function GeneRecord(epoch::Int, state::Int, parent::GeneRecord)
        self = new(epoch, state, parent.locus)
        self.parent = parent
        self
    end
end

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

function clean!(gdb::GeneDB, gids::AbstractArray)
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

    # After the familial tree is built, identify organisms with one (or less) offspring. Those organisms are
    # removed from the tree, and their offspring gets assigned to a parent of that offspring.
    for id in ids
        p = gdb[id].parent
        pid = p.id
        if !haskey(cs, id) && !in(id, gids)
            # This organism neither has offspring nor extant.
            id != pid && deleteat!(cs[pid], findfirst(id))
            drop!(gdb, id)
        elseif haskey(cs, id) && length(cs[id]) == 1
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

_getca(gdb::GeneDB, path) = length(path) > 0 ? gdb[maximum(path)] : (g = GeneRecord(0, 0, 0); g.id = 0; g)

function ca(gdb::GeneDB, id1::Int, id2::Int)
    path = IntSet(history(gdb, id1))
    intersect!(path, IntSet(history(gdb, id2)))
    _getca(gdb, path)
end

function mrca(gdb::GeneDB, gids::AbstractArray)
    path = IntSet(history(gdb, first(gids)))
    for gid in gids
        intersect!(path, IntSet(history(gdb, gid)))
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

    states = IntSet([gdb[h].state for h in history(gdb, first(gids))])
    for gid in gids
        union!(states, IntSet([gdb[h].state for h in history(gdb, gid)]))
    end

    length(states) - 1
end
