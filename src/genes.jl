export GeneDB,
       GeneRecord,
       UndefGene

abstract AbstractEvent
immutable Transmission <: AbstractEvent end
immutable Mutation <: AbstractEvent end
immutable Migration <: AbstractEvent
    from::Int32
    to::Int32
end
immutable MigrationAndMutation <: AbstractEvent
    from::Int32
    to::Int32
end

type GeneRecord
    # Primary information of lineages: Immutable throughout a simulation.
    epoch::Int
    state::Int
    event::AbstractEvent
    # id = 0 is a special type that the gene does not appear explicitly in a simulation.
    # In other words, id = 0 indicates the identity of a gene is unknown, so I can't say
    # gene(id = 0) == gene(id = 0).
    id::Int

    # Secondary information of lineages: This will be updated during a simulation.
    parent::GeneRecord
    ndescs::Int

    # id is intentionally left unspecified. That field is specified upon insertion in to GeneDB.
    GeneRecord(
        e::Int,
        s::Int,
        ev::AbstractEvent) = new(e, s, ev)
end

# a gene of unknown identity. This is also used to represent a parent of roots.
const UndefGene = GeneRecord(0, 0, Transmission())

## constructors
# initial genes
function GeneRecord(epoch::Int, state::Int)
    self = GeneRecord(epoch, state, Transmission())
    self.parent = UndefGene
    self.ndescs = 0
    self
end

# anything else
function GeneRecord(epoch::Int, parent::GeneRecord; state=-1, src=-1, dest=-1)
    if state <= 0 || state == parent.state
        if src <= 0 || dest <= 0
            self = GeneRecord(epoch, parent.state, Transmission())
        else
            self = GeneRecord(epoch, parent.state, Migration(src, dest))
        end
    else # state is different from its parent's.
        if src <= 0 || dest <= 0
            self = GeneRecord(epoch, state, Mutation())
        else
            self = GeneRecord(epoch, state, MigrationAndMutation(src, dest))
        end
    end
    self.parent = parent
    self.parent.ndescs += 1
    self.ndescs = 0
    self
end

type GeneDB
    currentid::Int
    data::Dict{Int, GeneRecord}

    GeneDB() = new(0, Dict{Int, GeneRecord}(0 => UndefGene))
end

Base.getindex(gdb::GeneDB, id::Int) = gdb.data[id]
Base.setindex!(gdb::GeneDB, record::GeneRecord, id::Int) = gdb.data[id] = record
Base.haskey(gdb::GeneDB, id::Int) = haskey(gdb.data, id)
function Base.delete!(gdb::GeneDB, id::Int)
    # delete GeneRecord from GeneDB only if it has no descendent.
    parent = gdb[id].parent
    parent.ndescs -= 1
    delete!(gdb.data, id)
end

Base.keys(gdb::GeneDB) = keys(gdb.data)

update!(gdb::GeneDB, id::Int, field::Symbol, value) = setfield!(gdb[id], field, value)
update!(gdb::GeneDB, record::GeneRecord, field::Symbol, value) = setfield!(record, field, value)

drop!(gdb::GeneDB, id::Int) = delete!(gdb.data, id)
drop!(gdb::GeneDB, record::GeneRecord) = delete!(gdb.data, record.id)

nextid!(gdb::GeneDB) = (gdb.currentid += 1; gdb.currentid)

function Base.insert!(gdb::GeneDB, record::GeneRecord)
    id = nextid!(gdb)
    record.id = id
    haskey(gdb, id) && error("Attempted to insert an existing record.")
    id < record.parent.id && error("Number of record exceeds supported maximum.")
    gdb[id] = record
    id
end

function transmit!(gdb::GeneDB, epoch::Int, pid::Int; state=-1, src=-1, dest=-1)
    pid = insert!(gdb, GeneRecord(epoch, gdb[pid], state=state, src=src, dest=dest))
end

function registerdescendant!(db, pid, did)
    if haskey(db, pid)
        push!(db[pid], did)
    else
        db[pid] = [did]
    end
end

function clean!{T}(gdb::GeneDB, cmin::T, cmax::T)
    ids = filter(x-> x < cmin || x > cmax, collect(keys(gdb)))
    sort!(ids, rev=true)
    # remove an element for UndefGene.
    pop!(ids)
    descs = Dict{T,Array{T,1}}()

    for id in cmin:cmax
        registerdescendant!(descs, gdb[id].parent.id, id)
    end

    for id in ids
        org = gdb[id]
        parent = org.parent
        ndescs = org.ndescs
        if ndescs == 0
            parent === UndefGene || (parent.ndescs -= 1)
            update!(gdb, org, :parent, org)
            drop!(gdb, id)
        elseif ndescs == 1 && isa(org.event, Transmission)
            length(descs[id]) == 1 || error()
            desc = gdb[descs[id][1]]
            update!(gdb, desc, :parent, parent)
            parent === UndefGene || registerdescendant!(descs, parent.id, desc.id)

            delete!(descs, id)
            drop!(gdb, id)
        else
            registerdescendant!(descs, parent.id, id)
        end
    end

    reverse!(ids)
    for id in ids
        haskey(gdb, id) || continue
        org = gdb[id]
        parent = org.parent
        while parent === UndefGene && org.ndescs == 1
            length(descs[org.id]) == 1 || error()
            desc = gdb[descs[org.id][1]]
            desc.parent = parent
            drop!(gdb, org.id)
            parent = org
            org = desc
        end
    end
end
