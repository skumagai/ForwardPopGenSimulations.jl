"""
Test for various types of chromosomes.

These types handles sex specific transmission of genetic materials.  The currently
implemented are Autosome, XChromosome, YChromosome, and Mitochondrion. Polyploidy
with n > 2 is not supported yet.

Actual storage of genertic data is delegated to ChromosomalStorage type.

In order to make holding organisms in a population simpler, chromosomes are not
parametrized by sex. This incurs unused ChromosomalChromosomes in sex chromosomes and
chromosomes.
"""

facts("""
    All organisms have two autosomes, and the pattern of transmissioin does not depend on
    sex.
    """) do
    context("A new chromosome can be initialized without parental chromosomes.") do

    end
    context("A new chromosome can be created with one maternal and one paternal chromosomes.") do

    end
end

facts("""
    A female contains two X-chromosome. A female offspring inherits each of X-chromosomes from
    its mother and father.
    """) do
    context("A new chromosome can be initialized without parental chromosomes.") do
    end

    context("A new chromosome can be created from one maternal and one paternal chromosomes.") do
    end
end

facts("A male only contains maternally inherited X-chromosome.") do
    context("A new chromosome can be initialized without parental chromosomes.") do
    end

    context("""
        A new chromosome can be created from one maternal and one paternal chromosomes.
        However, the paternal chromosome will be ignored.
        """) do
    end
end


facts("A female does not contain Y-chromosome.") do
    context("""
        A new chromosome can be initialized without parental chromosomes. However,
        the chromosome is intentionally left uninitialized.
    """) do
    end

    context("""
        A new chromosomes can be created from one maternal and one paternal chromosomes.
        However, both chromosomes will be ignored.
    """) do
    end
end

facts("A male only contains paternally inheritec Y-chromosome.") do
    context("A new chromosome can be initialized without parental chromosomes.") do
    end

    context("""
        A new chromosome can be created from one maternal and one paternal chromosomes.
        However, the maternal chromosome will be ignored.
    """) do
    end
end

facts("A female contains maternally inherited micohondrion.") do
    context("A new chromosome can be initialized without parental chromosomes.") do
    end

    context("""
        A new chromosome can be created from one maternal and one paternal chromosomes.
        However, the paternal chromosome will be ignored.
    """) do
    end
end

facts("""
    Although male contains maternally inherited mitochondria, a male parent never passes it
    to any offspring.
    """) do
    context("""
        A new chromosome can be initialized without parental chromosomes. However, the
        chromosome is intentionally left uninitialized.
    """) do
    end

    context("""
        A new chromosomes can be created from one maternal and one paternal chromosomes.
        However, both chromosomes will be ignored.
    """) do
    end
end
