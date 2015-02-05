"""
Tests for chromosomes and their associated methods.
"""

facts("""
    There are three types of chromosomes.  All of them are subtypes of an
    abstract type called Chromosome. The three subtypes are SegmentedChromosome,
    ContinuousChromosome, and SegmentedContinousChromosome.

    In addition to constructors,
    """) do

    @facts issubtype(SegmentedChromosome, Chromosome)
    @facts issubtype(ContinousChromosome, Chromosome)
    @facts issubtype(SegmentedContinousChromosome, Chromosome)
end

facts("""
    The first subtype of Chromosome is SegmentedChromosome. A object of this
    type consists of a limited number of loci. A locus then stores a subtype
    of Gene. A SegmentedChromosome permits only inter-locus recombination,
    so it is suitable for simulating relatively short segment of chromosome
    intersparsed by much longer intervals.
    """) do

    contexts("""
        At the beginning of a simulation, chromosomes are initialized with
        desired number of loci.
        """) do

        chr = SegmentedChromosome(BaseGene{Int}, 2)
        @fact typeof(chr) => SegmentedChromosome{BaseGene}
        @fact length(chr) => 2
    end

    context("""
        Afterward,

        """) do
    end
end

facts("""
    The second subtype of Chromosome is ContinousChromosome. A object of this
    type can have unlimited Genes. In fact, a gene is misnomer here. What genes
    represent is sites of mutation. Recombinations can occur between between
    any adjacent pairs of sites. Therefore, this type is suitable for the infinite-
    sites model. Internally, as the name of type suggests, sites are points on
    a real interval. Because of this design, it's impossible to preallocate fixed
    number of Genes at the beginning simulation. Likewise, it does not make sense
    to talk about the length of a chromosome.
    """) do

