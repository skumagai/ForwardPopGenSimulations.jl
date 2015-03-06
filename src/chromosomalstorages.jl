import Base: getindex, setindex!, enumerate, start, next, done

type ChromosomalStorage{T<:Union(AbstractVector, Associative)}
    data::T
    length::Int
end

ChromosomalStorage{T<:DenseVector}(::Type{T}, len::Integer) = ChromosomalStorage{T}(T(len), len)
ChromosomalStorage{T<:AbstractSparseVector}(::Type{T}, len::Integer) = ChromosomalStorage{T}(

getindex(chr::ChromosomalStorage, i::Integer) = data[i]
setindex!{Gene}(chr::ChromosomalStorage, g::Gene, i::Integer) = data[i] = g

enumerate(chr::ChromosomalStorage) = enumerate(chr.data)

enumeratenz(chr::ChromosomalStorage{DenseVector}) = enumerate(chr)
enumeratenz(chr::ChromosomalStorage{AbstractSparseVector}) = zip(chr.data.rowval, chr.data.nzval)
enumeratenz(chr::ChromosomalStorage{Associative}) = chr

length(chr::ChromosomalStorage) = chr.length

function deletedata!(c::ChromosomalStorage{AbstractSparseArray})
    # delete old data
    for i in findnz(d.data)[1]
        d.data[i] = 0
    end
end

function deletedata!(c::ChromosomalStorage{Associative})
    # delete old data
    for i in keys(d.data)
        delete!(d.data, i)
    end
end

function writedata!{C<:ChromosomalStorage}(d::C, p::C)
    # write new data
    for (i, j) in enumerate(p)
        d.data[i] = j
    end
    d.length = p.length
end

daughter!{C<:ChromosomalStorage{DenseArray}}(d::C, p::C) = writedata!(d, p)

function daughter!{C<:ChromosomalStorage}(d:C, p::C)
    deletedata!(d)
    writedata!(d, p)
end

function daughter{C<:ChromosomalStorage, I<:Integer}(p1::C, p2::C, Vector{I})
    ps = (p1, p2)
    c = 1

end

