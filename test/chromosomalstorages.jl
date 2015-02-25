"""
Test for actual storage of genetic materials.

Subtypes of abstract type ChromosomalStorage provides storage of genetic information. These
subtypes are parameterized by the type of genetic information, and genetic information can
be of any type without restriction (e.g., a linked list to track history along a lineage, or
simple integer to denote allele under the infinite-alleles model.

The first subtype of ChromosomalStorage is DiscreteStorage{T}. This type holds a fixed number
of genetic information stored in a vector. This type is dense in that all organisms store
information at the same sets of genetic sites.

The second subtype is ContinousStorage{T}. This type holds variable number of information.
Unlike the first subtype, different organisms can have information at different sets of sites.
This is implemented as SortedDict with keys representing positions.
"""

facts("Test DiscreteStorage{T}") do
    context("can be initialized with an initializer function.") do
    end

    context("can be created from a parental chromosome.") do
    end

    context("can be created from a pair of chromosomes and a function dealing recombination.") do
    end

    context("can be iterated.") do
    end

    context("can be indexed into.") do
    end

    context("can be queried about the length.") do
    end
end

facts("ContinuousStorage{T} implements ") do
    @pending :all_test => :not_implemented
end
