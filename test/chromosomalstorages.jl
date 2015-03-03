"""
Test for actual storage of genetic materials.

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

facts("ChromosomalStorage with dense backends implements these operations.") do
    nelems = 10

    context("It can be indexed into.") do
        chr = ChromosomalStorage{Vector, Int}(nelems)
        for i = 1:nelems
            chr[i] = i
            @fact chr[i] => i
        end
    end

    context("It can be iterated over all sites.") do
        chr = ChromosomalStorage{Vector, Int}(nelems)
        for i = 1:nelems
            chr[i] = i
        end
        for (i, site) in enumerate(chr)
            @fact site => i
        end
    end

    context("It can be queried about the length.") do
        for i = 2:2:10
            chr = ChromosomalStorage{Vector, Int}(i)
            @fact length(chr) => i
        end
    end

    context("It can be created from a parental chromosome.") do
        chr = ChromosomalStorage{Vector, Int}(nelems)
        for i = 1:nelems
            chr[i] = i
        end
        o = offspring(chr)
        for i = 1:nelems
            @fact o[i] => i
        end
        o[1] = 11
        @fact o[1] => not(chr[1])
    end

    context("It can be created from a pair of chromosomes and a list of recombination sites.") do
        recsites = [3, 6, 9]
        par1 = ChromosomalStorage{Vector, Int}(nelems)
        par2 = ChromosomalStorage{Vector, Int}(nelems)
        for i = 1:nelems
            par1[i] = i
            par2[i] = i + nelems
        end
        o = offspring(par1, par2, recsites)

        e1 = [1, 2, 3, 14, 15, 16, 7, 8, 9, 20]
        e2 = [11, 12, 13, 4, 5, 6, 17, 18, 19, 10]
        results = [o[i] for i = 1:nelems]
        @fact results => anyof(e1, e2)
    end

end

facts("""
ChromosomalStorage with sparse backends implements these operations. Because, some indices don't
have data associated with it, the type of returned values is nullable.
""") do
    nelems = 10
    sites = [2:2:10]

    context("It can be indexed into.") do
        chr = ChromosomalStorage{SparseMatrixCSC, Int}(nelems)
        for i = 1:sites
            if i in sites
                chr[i] = i
                @fact chr[i] => not(isnull)
                @fact get(chr[i]) => i
            else
                chr[i] = i
                @fact chr[i] => isnull
            end
        end
    end

    context("It can be iterated over all sites.") do
        chr = ChromosomalStorage{SparseMatrixCSC, Int}(nelems)
        for i =
            chr[i] = i
        end
        for (i, site) in enumerate(chr)
            @fact site => not(isnull)
            @fact get(site) => i
        end
    end

    context("It can be queried about the length.") do
        for i = 2:2:10
            chr = ChromosomalStorage{SparseMatrixCSC, Int}(i)
            @fact length(chr) => i
        end
    end

    context("It can be created from a parental chromosome.") do
        chr = ChromosomalStorage{SparseMatrixCSC, Int}(nelems)
        for i = 1:sites
            chr[i] = i
        end
        o = offspring(chr)
        for i = 1:nelems
            if i in sites
                @fact o[i] => not(isnull)
                @fact get(o[i]) => i
            else
                @fact o[i] => isnull
            end
        end
        o[1] = 11
        @fact o[1] => not(isnull)
        @fact chr[1] => isnull
        o[2] = 12
        @fact o[2] => not(isnull)
        @fact chr[2] => not(isnull)
        @fact get(chr[2]) => 2
        @fact get(o[2]) => get(chr[2])
    end

    context("It can be created from a pair of chromosomes and a list of recombination sites.") do
        recsites = [3, 6, 9]
        par1 = ChromosomalStorage{Vector, Int}(nelems)
        par2 = ChromosomalStorage{Vector, Int}(nelems)
        for i = 1:nelems
            if i in sites
                par1[i] = i
            else
                par2[i] = i + nelems
            end
        end
        o = offspring(par1, par2, recsites)

        e1 = [2, 15, 8, 9]
        e2 = [11, 13, 4, 6, 17, 19, 10]
        val = Vector{Int}(0)
        for 1 = 1:nelems
            isnull(o[i]) || push!(val, get(o[i]))
        end
        @fact vals => anyof(e1, e2)
    end

end

facts("ChromosomalStorage with dictionary backends implements these operations.") do
    len = 1.5
    sites = [1/3, 2/3, 6/5, 7/5]

    context("It can be indexed into.") do
        chr = ChromosomalStorage{SortedDict{Float64}, Float64}(len)
        for i = 1:sites
            chr[i] = i
            @fact chr[i] => not(isnull)
            @fact get(chr[i]) => i
        end
    end

    context("It can be iterated over all sites.") do
        chr = ChromosomalStorage{SortedDict{Float64}, Float64}(len)
        for i = 1:sites
            chr[i] = i
        end
        for (i, site) in enumerate(chr)
            @fact site => not(isnull)
            @fact get(site) => sites[i]
        end
    end

    context("It can be queried about the length.") do
        for i = 2:2:10
            chr = ChromosomalStorage{SortedDict{Float64}, Float64}(i * len)
            @fact length(chr) => i * len
        end
    end

    context("It can be created from a parental chromosome.") do
        chr = ChromosomalStorage{SortedDict{Float64}, Float64}(len)
        for i = 1:sites
            chr[i] = i
        end
        o = offspring(chr)
        for i = 1:sites
            @fact o[i] => not(isnull)
            @fact get(o[i]) => i
        end
        o[1.0] = 11.0
        @fact o[1.0] => not(isnull)
        @fact get(o[1.0]) => 11.0
        @fact chr(1.0) => isnull
    end

    context("It can be created from a pair of chromosomes and a list of recombination sites.") do
        recsites = [0.49, 1.01]
        sites = [1/3, 2/3, 6/5, 7/5]
        sites2 = [1/4, 2/4, 3/4, 5/4]
        par1 = ChromosomalStorage{SortedDict{Float64}, Float64}(len)
        par2 = ChromosomalStorage{SortedDict{Float64}, Float64}(len)
        for i = 1:sites
            par1[i] = i
        end
        for i = 1:sites2
            par2[i] = i
        end
        o = offspring(par1, par2, recsites)

        e1 = [1/3, 2/4, 3/4, 6/5, 7/5]
        e2 = [1/4, 2/3, 5/4]
        val = Vector{Float64}(0)
        for i in o
            @fact i => not(isnull)
            push!(val, get(i))
        end
        @fact val => anyof(e1, e2)
    end

end
