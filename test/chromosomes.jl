"""
Tests for chromosomes and their associated methods.
"""

facts("""
    There are three types of chromosomes.  All of them are subtypes of an
    abstract type called Chromosome. The three subtypes are SegmentedChromosome,
    ContinuousChromosome, and SegmentedContinousChromosome.
    """) do

    @facts issubtype(SegmentedChromosome, Chromosome)
    @facts issubtype(ContinousChromosome, Chromosome)
    @facts issubtype(SegmentedContinousChromosome, Chromosome)
end

facts("""
    The subtype of Chromosome is SegmentedChromosome. A object of this
    type consists of a limited number of loci. A locus then stores a
    subtype of Gene. A SegmentedChromosome permits only inter-locus
    recombination, so it is suitable for simulating relatively short
    segment of chromosome intersparsed by much longer intervals.
    """) do

    contexts("""


    """) do

        chr = SegmentedChromosome(BaseGene{Int}, 2)
        @fact typeof(chr) => SegmentedChromosome{BaseGene}
        @fact length(chr) => 2
    end
end
