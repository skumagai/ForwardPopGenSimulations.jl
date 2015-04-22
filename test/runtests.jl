import SymmetricDominance
using Base.Test

const sd = SymmetricDominance

const g1 = sd.Gene(1, 1)
const g2 = sd.Gene(1, 2)
const g3 = sd.Gene(2, 3)

@test sd.isidbystate(g1, g2) == true
@test sd.isidbystate(g1, g3) != true
@test sd.isidbystate(g2, g3) != true

const ldb = Array(sd.LineageRecord, 0)
append!(ldb, [sd.LineageRecord(0, 1), sd.LineageRecord(0, 2), sd.LineageRecord(0, 3)])
midx, g4 = sd.mutate(g1, 2, ldb, 3)
@test midx == 3
@test g4 == sd.Gene(3, length(ldb))
@test ldb[end] == sd.LineageRecord(1, 3)

const g5 = sd.Gene(3, length(ldb) + 1)
push!(ldb, sd.LineageRecord(4, 5))
@test 1 == sd.getancestor(g5, ldb)

const ps = Array(Int, 2)
for _ = 1:10000
    sd.getparentids!(ps, 10)
    @test ps[1] != ps[2]
end

const o = sd.Organism([sd.Gene(0, 0) sd.Gene(0, 0);
                       sd.Gene(0, 0) sd.Gene(0, 0);
                       sd.Gene(0, 0) sd.Gene(0, 0)])
@test length(o) == 3

empty!(ldb)
append!(ldb, [sd.LineageRecord(0, 0), sd.LineageRecord(0, 0), sd.LineageRecord(0, 0), sd.LineageRecord(0, 0),
              sd.LineageRecord(1, 1), sd.LineageRecord(2, 1), sd.LineageRecord(3, 1), sd.LineageRecord(4, 1),
              sd.LineageRecord(5, 2), sd.LineageRecord(5, 2), sd.LineageRecord(6, 2), sd.LineageRecord(6, 2),
              sd.LineageRecord(7, 2), sd.LineageRecord(7, 2), sd.LineageRecord(8, 2), sd.LineageRecord(8, 2)])

const pops = [sd.Organism([sd.Gene(0, 9) sd.Gene(0, 10); sd.Gene(0, 13) sd.Gene(0, 14)]) sd.Organism([sd.Gene(0, 9) sd.Gene(0, 11); sd.Gene(0, 13) sd.Gene(0, 14)]);
              sd.Organism([sd.Gene(0, 9) sd.Gene(0, 10); sd.Gene(0, 13) sd.Gene(0, 14)]) sd.Organism([sd.Gene(0, 11) sd.Gene(0, 12); sd.Gene(0, 15) sd.Gene(0, 16)])]
@test sd.hascoalesced(pops, ldb, 1)
@test !sd.hascoalesced(pops, ldb, 2)
