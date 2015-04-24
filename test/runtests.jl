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

const ldbs = [[sd.LineageRecord(0, 0), sd.LineageRecord(0, 0), sd.LineageRecord(1, 1), sd.LineageRecord(2, 1),
               sd.LineageRecord(3, 2), sd.LineageRecord(3, 2), sd.LineageRecord(4, 2), sd.LineageRecord(4, 2),
               sd.LineageRecord(5, 3), sd.LineageRecord(5, 3), sd.LineageRecord(6, 3), sd.LineageRecord(6, 3),
               sd.LineageRecord(7, 3), sd.LineageRecord(7, 3), sd.LineageRecord(8, 3), sd.LineageRecord(8, 3)]
              for _ = 1:2]

pops = Array(sd.Organism, 2, 2)

pops[1, 1] = sd.Organism([sd.Gene(1, 9) sd.Gene(2, 10); sd.Gene(3, 9) sd.Gene(4, 13)])
pops[2, 1] = sd.Organism([sd.Gene(5, 11) sd.Gene(6, 12); sd.Gene(7, 10) sd.Gene(8, 14)])

pops[1, 2] = sd.Organism([sd.Gene(9, 11) sd.Gene(10, 15); sd.Gene(11, 13) sd.Gene(12, 14)])
pops[2, 2] = sd.Organism([sd.Gene(13, 12) sd.Gene(14, 16); sd.Gene(15, 15) sd.Gene(16, 16)])

@test sd.hascoalesced(pops, ldbs, 1, 1)
@test !sd.hascoalesced(pops, ldbs, 2, 1)
@test !sd.hascoalesced(pops, ldbs, 1, 2)
@test sd.hascoalesced(pops, ldbs, 2, 2)


sdata3 = Array(Int, 2, 2, 2)
sdata3[1, 1, 1] = 1
sdata3[2, 1, 1] = 5
sdata3[1, 2, 1] = 3
sdata3[2, 2, 1] = 7
sdata3[1, 1, 2] = 2
sdata3[2, 1, 2] = 6
sdata3[1, 2, 2] = 4
sdata3[2, 2, 2] = 8
@test sd.toarray(pops[:,1], :state) == sdata3
sdata3[1, 1, 1] = 9
sdata3[2, 1, 1] = 13
sdata3[1, 2, 1] = 11
sdata3[2, 2, 1] = 15
sdata3[1, 1, 2] = 10
sdata3[2, 1, 2] = 14
sdata3[1, 2, 2] = 12
sdata3[2, 2, 2] = 16
@test sd.toarray(pops[:,2], :state) == sdata3


const acs, gcs, hcs = sd.counts(pops[:,1])
@test acs == [Dict{Int,Int}(1=>1, 2=>1, 5=>1, 6=>1), Dict{Int,Int}(3=>1, 4=>1, 7=>1, 8=>1)]
@test gcs == [Dict{Tuple{Int,Int}, Int}((1,2)=>1, (5,6)=>1), Dict{Tuple{Int,Int}, Int}((3,4)=>1, (7,8)=>1)]
@test hcs == Dict{Tuple{Int,Int}, Int}((1,3)=>1, (2,4)=>1, (5,7)=>1, (6,8)=>1)

const afs, gfs, hfs = sd.spectra(pops[:,1])
@test afs == [Dict{Int,Int}(1=>4), Dict{Int,Int}(1=>4)]
@test gfs == [Dict{Int,Int}(1=>2), Dict{Int,Int}(1=>2)]
@test hfs == Dict{Int,Int}(1=>4)

@test sd.history(ldbs[1], 1) == [0, 1]
@test sd.history(ldbs[1], 2) == [0, 2]
@test sd.history(ldbs[1], 3) == [0, 1, 3]
@test sd.history(ldbs[1], 4) == [0, 2, 4]
@test sd.history(ldbs[1], 5) == [0, 1, 3, 5]
@test sd.history(ldbs[1], 6) == [0, 1, 3, 6]
@test sd.history(ldbs[1], 7) == [0, 2, 4, 7]
@test sd.history(ldbs[1], 8) == [0, 2, 4, 8]
@test sd.history(ldbs[1], 9) == [0, 1, 3, 5 ,9]
@test sd.history(ldbs[1], 10) == [0, 1, 3, 5, 10]
@test sd.history(ldbs[1], 11) == [0, 1, 3, 6, 11]
@test sd.history(ldbs[1], 12) == [0, 1, 3, 6, 12]
@test sd.history(ldbs[1], 13) == [0, 2, 4, 7, 13]
@test sd.history(ldbs[1], 14) == [0, 2, 4, 7, 14]
@test sd.history(ldbs[1], 15) == [0, 2, 4, 8, 15]
@test sd.history(ldbs[1], 16) == [0, 2, 4, 8, 16]

ds = sd.distances(pops[:,1], ldbs)
@test ds == [1 1 2 2;
             1 1 5 4;
             1 1 6 4;
             1 2 5 4;
             1 2 6 4;
             1 5 6 2;
             2 3 4 8;
             2 3 7 2;
             2 3 8 8;
             2 4 7 8;
             2 4 8 2;
             2 7 8 8]
ds = sd.distances(pops[:,2], ldbs)
@test ds == [1 9 10 8;
             1 9 13 2;
             1 9 14 8;
             1 10 13 8;
             1 10 14 2;
             1 13 14 8;
             2 11 12 2;
             2 11 15 4;
             2 11 16 4;
             2 12 15 4;
             2 12 16 4;
             2 15 16 2]

