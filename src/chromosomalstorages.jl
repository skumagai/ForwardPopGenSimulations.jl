import Base: getindex, setindex!, enumerate, start, next, done, call, length
import DataStructures: SortedDict

export Chromosome, DenseChromosome, IntervalChromosome, setoffspring!

abstract Chromosome{Gene}

type DenseChromosome{Gene} <: Chromosome{Gene}
    data::Vector{Gene}
end
call{T}(::Type{DenseChromosome{T}}, len::Integer) = DenseChromosome{T}(Vector{T}(len))

type IntervalChromosome{Gene, Tk<:Real} <: Chromosome{Gene}
    data::SortedDict{Tk, Gene}
    n::Int
end
call{T, Tk<:Real}(::Type{IntervalChromosome{T, Tk}}, len::Integer) =
    IntervalChromosome{T,Tk}(SortedDict{Tk,T}())

getindex(chr::Chromosome, i::Integer) = chr.data[i]
setindex!{Gene}(chr::Chromosome{Gene}, g::Gene, i::Integer) = chr.data[i] = g

enumerate(chr::Chromosome) = enumerate(chr.data)

length(chr::DenseChromosome) = length(chr.data)
length(chr::IntervalChromosome) = chr.n

function _writedata!{C<:Chromosome}(d::C, p::C)
    # write new data
    for (i, j) in enumerate(p)
        d.data[i] = j
    end
end

writedata!{C<:DenseChromosome}(d::C, p::C) = _writedata!(d, p)

function writedata!{C<:IntervalChromosome}(d::C, p::C)
    empty!(d.data)
    _writedata!(d, p)
    d.n = p.n
end

setoffspring!{C<:Chromosome}(d::C, p::C) = writedata!(d, p)

function setoffspring!{C<:Chromosome, I<:Integer}(d::C, ps::(C, C), recs::Vector{I})
    i0 = 1
    selected = 1
    for recblock = 1:length(recs)
        for i = i0:recs[recblock]
            d.data[i] = ps[selected].data[i]
        end
        selected = 3 - selected
        i0 = recs[recblock] + 1
    end
    for i = i0:length(d)
        d.data[i] = ps[selected].data[i]
    end
end

