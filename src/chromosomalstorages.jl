import Base: getindex, setindex!, enumerate, start, next, done, call
import DataStructures: SortedDict

export Chromosome, DenseChromosome, IntervalChromosome

abstract Chromosome{Gene}

type DenseChromosome{Gene} <: Chromosome{Gene}
    data::Vector{Gene}
end
call{T}(::Type{DenseChromosome{T}}, len::Integer) = DenseChromosome{T}(Vector{T}(len))

type IntervalChromosome{Gene, Tk<:Real} <: Chromosome{Gene}
    data::SortedDict{Tk, Gene}
    n::Int
end
call{T, Tk<:Real}(::Type{IntervalChromosome{T, TK}}, len::Integer) =
    IntervalChromosome{T,Tk}(SortedDict{Tk,T}())

getindex(chr::Chromosome, i::Integer) = chr.data[i]
setindex!{Gene}(chr::Chromosome{Gene}, g::Gene, i::Integer) = chr.data[i] = g

enumerate(chr::Chromosome) = enumerate(chr.data)

length(chr::Chromosome) = length(chr)
length(chr::IntervalChromosome) = chr.n

function _writedata!{C<:Chromosome}(d::C, p::C)
    # write new data
    for (i, j) in enumerate(p)
        d.data[i] = j
    end
end

writedata!{C<:Chromosome}(d::C, p::C) = _writedata(d, p)

function writedata!{C<:IntervalChromosome}(d::C, p::C)
    _writedata!(d, p)
    d.n = p.n
end

function daughter!{C<:Chromosome}(d:C, p::C)
    empty!(d.data)
    writedata!(d, p)
end

function daughter!{C<:Chromosome, I<:Integer}(d::C, ps::(C, C), recs::Vector{I})
    empty!(d.data)
end

