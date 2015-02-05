"""
Tests for genes and its associated methods.

At minimum, a gene has state, of an arbitrary type. Unlike fwdpp but like
simuPOP, genes do not hold fitness.  A main reason for this design is that
fitness is a property of an individual in the context of specific environment.
For example, an allele may be beneficial for an individual in a subpopulation.
But the same allele may be detrimental to the same individual in a different
subpopulation.  If fitness is encoded within a gene/allele, a value of fitness
may be needed to update upon migration. This introduces unnecessary complexity
between demographic event (migration) and a fundamental building block of
an organism.

Furthermore, fitness of individuals may not be cleanly decomposed to loci.  If
a part of fitness is attributed to interactions among loci, how can a gene store
the locus-specific fitness. There is no clear answer for this.

Instead of storing fitness in gene, users of this package are required to provide
mapping between states of gene to fitness externally.
"""

facts("""
    The most basic subtype of an abstract type Gene is BaseGene{T}, and its
    type parameter determines type of "state" field. The only other field
    this type has is "site" field.
    """) do

    @fact issubtype(BaseGene, Gene) => true
    gene = BaseGene(0, 0)
    @fact state(gene) => 0
    @fact typeof(state(gene)) => Int
    @fact site(gene) => 0
    @fact typeof(site(gene)) => Int
    gene = BaseGene(1, 2)
    @fact state(gene) => 1
    @fact site(gene) => 2
    @fact_throws epoch(gene)
    @fact_throws ancestor(gene)
    @fact_throws source(gene)
    @fact_throws destination(gene)
end

facts("""
    Another subtype of Gene can record when the gene originated. In addition to
    "state" and "site" fields, TimedGene{T} also has "epoch" field of type
    Int.
    """) do
    @fact issubtype(TimedGene, Gene) => true
    @fact issubtype(TimedGene, Gene) => true
    @fact issubtype(TimedGene, TimedGene) => true
    gene = TimedGene(2, 0, 3)
    @fact typeof(gene) => TimedGene{Int}
    @fact state(gene) => 2
    @fact typeof(state(gene)) => Int
    @fact site(gene) => 0
    @fact typeof(site(gene)) => Int
    @fact epoch(gene) => 3
    @fact_throws ancestor(gene)
    @fact_throws source(gene)
    @fact_throws destination(gene)
end

facts("""
    TrackedGene can record ancestral gene in addition to "state", "site",
    and "epoch" fields. TrackedGene itself is an abstract type currently three
    subtypes.
    """) do
    @fact issubtype(TrackedGene, Gene)
    context("""The first gene does not ancestral gene.
               This is represented as a special type OriginalGene.""") do
        @fact issubtype(OriginalGene, TrackedGene) => true

        gene = OriginalGene(2, 0, 3)
        @fact state(gene) => 2
        @fact typeof(state(gene)) => Int
        @fact site(gene) => 0
        @fact typeof(site(gene)) => Int
        @fact epoch(gene) => 3
        @fact ancestor(gene) => gene
        @fact typeof(ancestor(gene)) => OriginalGene{Int, Int}
        @fact_throws source(gene)
        @fact_throws destination(gene)
    end
    context("""
        If the last event is mutation, a gene is represented as MutatedGene.
        The new and ancestral genes must share the site. In practice,
        this type is useless for simulations using the infintie-sites model
        of mutation as a new mutation always arises in a new site.
        """) do
        @fact issubtype(MutatedGene, TrackedGene)
        orig = OriginalGene(2, 0, 3)
        gene = MutatedGene(3, 4, orig)
        @fact state(gene) => 3
        @fact typeof(state(gene)) => Int
        @fact site(gene) => 0
        @fact site(gene) => site(orig)
        @fact typeof(site(gene)) => Int
        @fact epoch(gene) => 4
        @fact ancestor(gene) => orig
        @fact issubtype(typeof(ancestor(gene)), TrackedGene)
        @fact_throws source(gene)
        @fact_throws destination(gene)
    end
    context("""
        If the last event is migration, a gene is represented as MigratedGene.
        This gene must have the same site and state as its direct ancestor.
        """) do

        @fact issubtype(MigratedGene, TrackedGene)
        orig = OriginalGene(2, 1, 3)
        gene = MigratedGene(4, orig, 1, 2)
        @fact state(gene) => state(orig)
        @fact typeof(state(gene)) => Int
        @fact site(gene) => site(orig)
        @fact typeof(site(gene)) => Int
        @fact epoch(gene) => 4
        @fact ancestor(gene) => orig
        @fact issubtype(typeof(ancestor(gene)), TrackedGene)
        @fact source(gene) => 1
        @fact destination(gene) => 2
    end
end

facts("""
    Generic function `mutate` returns a new mutated gene if applicable. This
    function should not affect site..
    """) do
    contexts("For BaseGene, it simply returns a new gene.") do
        gene = BaseGene(0, 0)
        @fact state(gene) => 0
        newgene = mutate(gene, 1)
        @fact typeof(newgene) => BaseGene
        @fact state(newgene) => 1
        @fact site(newgene) => site(gene)
    end
    context("For TimedGene, it returns a new gene with new timing.") do
        gene = TimedGene(0, 2, 1)
        @fact state(gene) => 0
        @fact site(gene) => 2
        @fact epoch(gene) => 1
        newgene = mutate(gene, 2, 3)
        @fact issubtype(typeof(newgene), TimedGene) => true
        @fact state(newgene) => 2
        @fact site(newgene) => site(gene)
        @fact epoch(newgene) => 3
    end
    context("For OriginalGene, mutate must return MutatedGene.") do
        gene = OriginalGene(0, 0, 1)
        @fact state(gene) => 0
        @fact site(gene) => 0
        @fact epoch(gene) => 1
        newgene = mutate(gene, 2, 3)
        @fact typeof(newgene) => MutatedGene
        @fact state(newgene) => 2
        @fact site(newgene) => site(gene)
        @fact epoch(newgene) => 3
        @fact ancestor(newgene) => gene
    end
    context("For MutatedGene, mutate must return MutatedGene.") do
        orig = OriginalGene(0, 3, 1)
        @fact state(orig) => 0
        @fact site(orig) => 3
        @fact epoch(orig) => 1
        gene = MutatedGene(2, 3, orig)
        @fact state(gene) => 2
        @fact site(gene) => site(orig)
        @fact epoch(gene) => 3
        @fact ancestor(gene) => orig
        newgene = mutate(gene, 4, 5)
        @fact issubtype(typeof(newgene), MutatedGene) => true
        @fact state(genenew) => 4
        @fact site(newgene) => site(gene)
        @fact epoch(newgene) => 5
        @fact ancestor(newgene) => gene
    end
    context("For MigratedGene, mutate must return MutatedGene.") do
        orig = OriginalGene(0, 3, 1)
        @fact state(orig) => 0
        @fact site(orig) => 3
        @fact epoch(orig) => 1
        gene = MigratedGene(2, orig, 1, 2)
        @fact state(gene) => 0
        @fact site(gene) => site(orig)
        @fact epoch(gene) => 2
        @fact source(gene) => 1
        @fact destination(gene) => 2
        newgene = mutate(gene, 3, 4)
        @fact issubtype(typeof(newgene), MutatedGene) => true
        @fact state(newgene) => 3
        @fact site(newgene) => site(gene)
        @fact epoch(newgene) => 4
    end
end

facts("Generic function `migrate` returns a new mutated gene if applicable.") do
    contexts("For BaseGene, this function returns the gene unmodified.") do
        gene = BaseGene(0)
        newgene = migrate(gene, 4, 1, 2)
        @fact typeof(newgene) => BaseGene
        @fact newgene = exactly(gene)
    end
    context("For TimedGene, this function returns the gene unmodified.") do
        gene = TimedGene(0, 1)
        newgene = migrate(gene, 2, 1, 2)
        @fact typeof(newgene) => TimedGene
        @fact newgene => exactly(gene)
    end
    context("For OriginalGene, this function returns MigratedGene.") do
        gene = OriginalGene(0, 1, 2)
        newgene = migrate(gene, 4, 1, 2)
        @fact issubtype(typeof(newgene), MigratedGene)
        @fact state(newgene) => state(gene)
        @fact site(newgene) => site(gene)
        @fact epoch(newgene) => 4
        @fact sourece(newgene) => 1
        @fact destination(newgene) => 2
    end
    context("For MutatedGene, this function returns MigratedGene.") do
        orig = OriginalGene(0, 1, 2)
        gene = MutatedGene(3, 4, orig)
        @fact state(gene) => 3
        @fact site(gene) => site(orig)
        @fact epoch(gene) => 4
        @fact ancestor(gene) => orig
        newgene = migrate(gene, 5, 6, 7)
        @fact state(newgene) => state(gene)
        @fact site(newgene) => site(orig)
        @fact epoch(newgene) => 5
        @fact ancestor(newgene) => gene
        @fact source(newgene) => 6
        @fact destination(newgene) => 7
    end
end

