import Base: getindex, setindex!, enumerate, call, length
import DataStructures: SortedDict, pastendtoken, deref, deref_value, startof, advance

export Chromosome, DenseChromosome, IntervalChromosome, setoffspring!, inherit

abstract Chromosome{Gene}

type DenseChromosome{Gene} <: Chromosome{Gene}
    data::Vector{Gene}
end
call{T}(::Type{DenseChromosome{T}}, len::Integer) = DenseChromosome{T}(Vector{T}(len))

type IntervalChromosome{Gene, Tk<:Real} <: Chromosome{Gene}
    data::SortedDict{Tk, Gene}
    n::Tk
end
call{T, Tk<:Real}(::Type{IntervalChromosome{T, Tk}}, len::Tk) =
    IntervalChromosome{T,Tk}(SortedDict(Dict{Tk,T}()), len)

getindex(chr::DenseChromosome, i::Integer) = chr.data[i]
function getindex(chr::IntervalChromosome, i::Real)
    val = find(chr.data, i)
    if val == pastendtoken(chr.data)
        Nullable{typeof(i)}()
    else
        Nullable(deref_value(val))
    end
end
setindex!{Gene}(chr::DenseChromosome{Gene}, g::Gene, i::Integer) = chr.data[i] = g
setindex!{Gene}(chr::IntervalChromosome{Gene}, g::Gene, i::Real) = chr.data[i] = g

enumerate(chr::DenseChromosome) = enumerate(chr.data)

enumerate(chr::IntervalChromosome) = chr.data

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

function setoffspring!{C<:DenseChromosome, I<:Integer}(d::C, ps::(C, C), selected::Int, recs::Vector{I})
    i0 = 1
    for recblock = 1:length(recs)
        for i = i0:recs[recblock]
            d.data[i] = inherit(ps[selected].data[i])
        end
        selected = 3 - selected
        i0 = recs[recblock] + 1
    end
    for i = i0:length(d)
        d.data[i] = inherit(ps[selected].data[i])
    end
end

function setoffspring!{C<:IntervalChromosome, I<:Real}(d::C, ps::(C, C), selected::Int, recs::Vector{I})
    empty!(d.data)
    # Assuming position is always positive.
    pos0 = -1.0
    iends = (pastendtoken(ps[1].data), pastendtoken(ps[2].data))
    is = [startof(ps[1].data), startof(ps[2].data)]
    for recsite in [recs; length(ps[1])]
        i = is[selected]
        while i != iends[selected]
            pos, val = deref(i)
            if pos > recsite
                # If the current position is past the next recombination site, switch to another parantal chromosome.

                # Cache the current state
                is[selected] = i
                selected = 3 - selected
                i = is[selected]
                pos0 = recsite
                break
            elseif pos <= pos0
                i = advance(i)
            else
                d.data[pos] = inherit(val)
                i = advance(i)
            end
        end
    end
end

inherit(i::Number) = i
