export ChromosomeType,

       Autosome,
       XChromosome,
       YChromosome,
       Mitochondrion,

       value

immutable ChromosomeType
    data::Int
end

@inline value(c::ChromosomeType) = c.data

const Autosome = ChromosomeType(1)
const XChromosome = ChromosomeType(2)
const YChromosome = ChromosomeType(3)
const Mitochondrion = ChromosomeType(4)
