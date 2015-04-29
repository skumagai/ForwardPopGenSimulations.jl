# SymmetricDominance

SymmetricDominance performs forward-in-time population genetics simulations under
symmetric overdominance and underdominance.
Notable characteristics of this simulator include:
- simulations with a selected locus under symmetric overdominance and underdominance,
- arbitrary number of neutral loci,
- recording and constructing a tree of mutations,
- terminating simulations after k-th coalescence of the entire population (turnover), and
- constant population size.

## Symmetric overdominance and underdominance

Simulations are performed with a single selected locus under either symmetric overdominance or underdominance.
This selected locus is located at one end of an autosomal chromosome.
Overdominance assumes that organisms heterozygous at the selected locus have higher fitness than ones
are homozygous at the same locus.
Symmetric overdominance further assumes that fitness is only affected by heterozygous or not;
one allele is equaly fit as any other allele.
On the contrary, homozygotes have higher fitness than heterozygotes under underdominance model.

In this simulator, a pair of distinct parents are randomly chosen to create one offspring
from the entire parental population regardless of their fitness.
One of two parental genes are randomly passed onto the offspring with
possibility of a mutation.
When a mutation occurs, the offspring receives a mutated gene, which is distinguishable from any other extant genes in
the population.
Once the offspring is formed, it undergoes a test to determine the organism survive until maturity.
If determined to die before adulthood, a new offspring is formed from another random pair of offspring.
Otherwise, the offspring is stored in the offspring population.

## Arbitrary number of neutral loci

As mentioned above, this simulator places a selected locus at one end of an autosomal chromosome.
This selected locus is then followed by arbitrary number of neutral loci.
When passed from parents to offspring, recombinations occur between every adjacent pairs of loci
at the site-specific rates.

Mutations also occur at the locus-specific rates. Similar to mutations in the selected locus,
a new allele is distinguishable from any other alleles in the population.

## Recording and constructing a tree of mutations

Whenever a mutation occurs the simulator records a unique ID, time of mutation, allelic state,
and parental gene of this mutation.
In general, allelic state and unique ID of a gene do not match, because each gene is assigned
a distinguishable ID at the beginning of simulation and immediately after turnover events.
However, allelic states are not modified.
This allows detection of turnover events while keeping allelic states of the population intact.
The detection is achieved by collecting the ancestral IDs of all organisms in the current population.
A turnover event occurs when only one ancestral ID remains in the population.
Using the same information, the simulator can also build a tree of mutations.

## Termination of simulations upon k-th turnover event

Directly building on top of the last point, the simulator can terminate a simulation immediately
after all loci undergo k-th turnover events.

## Constant population size

For now, population size is constant.

