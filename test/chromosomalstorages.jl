"""
Testing storage of genetic materials.

ChromosomalStorage provides linearly arranged (vector) storage for genetic information.
Such information can be either dense or sparse, and the type of the information can be
arbitrary.

In terms of design, ChromosomalStorage is a type taking two parametr. The first type argument
specifies actual storage mode (e.g., Vector), and all backends are 1-dimensional. The second type
is the actual type of genetic information. There is no restriction as to what type can be used
here.

The three linear storage backends are provided for ChromosomalStorage. The first and second
storages are, respectively, backed by a dense and sparse vector. These backends use 1-based
integer indexing. The last backends is sparse and allows floating point positions. This is
supported by Dict{R<:Real,T}.
"""

import ForwardPopGenSimulations: Chromosome, DenseChromosome, IntervalChromosome

facts("Chromosomal has two (parameterized) subtypes.") do
    @fact DenseChromosome <: Chromosome => true
    @fact IntervalChromosome <: Chromosome => true
    @fact DenseChromosome{Int} <: Chromosome => true
    @fact IntervalChromosome{Int} <: Chromosome{Int} => true
    @fact IntervalChromosome{Int, Float64} <: Chromosome{Int} => true
    @fact DenseChromosome{Int} <: Chromosome{Int} => true
    @fact IntervalChromosome{Int} <: Chromosome{Int} => true
    @fact IntervalChromosome{Int, Float64} <: Chromosome{Int} => true
end

facts("DenseChromosome implements these operations.") do
    nelems = 10
    chr = DenseChromosome{Int}(nelems)
    for i = 1:nelems
        chr[i] = i
    end

    context("It can be indexed into.") do
        for i = 1:nelems
            @fact chr[i] => i
        end
    end

    context("It can be iterated over all genetic information.") do
        for (i, site) in enumerate(chr)
            @fact site => i
        end
    end

    context("It can be queried about the length.") do
        for i = 2:2:10
            c = DenseChromosome{Int}(i)
            @fact length(c) => i
        end
    end

    context("It can be created from a parental chromosome.") do
        o = DenseChromosome{Int}(nelems)

        setoffspring!(o, chr)
        for i = 1:nelems
            @fact o[i] => i
        end
        o[1] = 11
        @fact o[1] => not(chr[1])
    end

    context("""
    It can be created from a pair of chromosomes and a list of recombination sites. By
    convension, a offspring's chromosome begins with the first parental chromosome.
    """) do
        recsites = [3, 6, 9]
        par1 = DenseChromosome{Int}(nelems)
        par2 = DenseChromosome{Int}(nelems)
        for i = 1:nelems
            par1[i] = i
            par2[i] = i + nelems
        end
        ex1 = [1, 2, 3, 14, 15, 16, 7, 8, 9, 20]
        ex2 = [11, 12, 13, 4, 5, 6, 17, 18, 19, 10]

        for (pars, ex) in zip(((par1, par2), (par2, par1)), (ex1, ex2))
            o = DenseChromosome{Int}(nelems)
            setoffspring!(o, pars, recsites)
            vals = [o[i] for i = 1:nelems]
            @fact vals => ex
        end
    end

end

facts("IntervalChromosome backends implements these operations.") do
    len = 1.5
    sites = [1/3, 2/3, 6/5, 7/5]

    chr = IntervalChromosome{Float64, Float64}(len)
    for i in sites
        chr[i] = i
    end

    context("It can be indexed into. Return values are wrapped in Nullable.") do
        for i in sites
            @fact chr[i] => not(isnull)
            @fact get(chr[i]) => i
        end
    end

    context("It can be iterated over all explicitly stored genetic information.") do
        for (i, site) in enumerate(chr)
            @fact site => i
        end
    end

    context("It can be queried about the length.") do
        for i = 2:2:10
            c = IntervalChromosome{Float64, Float64}(i * len)
            @fact length(c) => i * len
        end
    end

    context("It can be created from a parental chromosome.") do
        o = IntervalChromosome{Float64, Float64}(len)
        setoffspring!(o, chr)
        for i in sites
            @fact o[i] => not(isnull)
            @fact get(o[i]) => i
        end
        o[1.0] = 11.0
        @fact o[1.0] => not(isnull)
        @fact get(o[1.0]) => 11.0
        @fact chr[1.0] => isnull
    end

    context("""
    It can be created from a pair of chromosomes and a list of recombination sites. By
    convension, offspring's chromosome begins with the first parental chromosome.
    """) do
        recsites = [0.49, 1.01]
        # sites = [1/3, 2/3, 6/5, 7/5]
        sites2 = [1/4, 2/4, 3/4, 5/4]
        par1 = IntervalChromosome{Float64, Float64}(len)
        par2 = IntervalChromosome{Float64, Float64}(len)
        for i in sites
            par1[i] = i
        end
        for i in sites2
            par2[i] = i
        end
        ex1 = [1/3, 2/4, 3/4, 6/5, 7/5]
        ex2 = [1/4, 2/3, 5/4]

        for (pars, ex) in zip(((par1, par2), (par2, par1)), (ex1, ex2))
            o = IntervalChromosome{Float64, Float64}(len)
            setoffspring!(o, pars, recsites)
            vals = Vector{Float64}(0)
            for (pos, val) in enumerate(o)
                @fact pos => val
                push!(vals, pos)
            end
            @fact vals => ex
        end
    end
end
