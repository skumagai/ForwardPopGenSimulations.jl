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

const pops = [sd.Organism([sd.Gene(1, 9) sd.Gene(2, 10); sd.Gene(3, 13) sd.Gene(4, 14)]) sd.Organism([sd.Gene(5, 9) sd.Gene(6, 11); sd.Gene(7, 13) sd.Gene(8, 14)]);
              sd.Organism([sd.Gene(9, 9) sd.Gene(10, 10); sd.Gene(11, 13) sd.Gene(12, 14)]) sd.Organism([sd.Gene(13, 11) sd.Gene(14, 12); sd.Gene(15, 15) sd.Gene(16, 16)])]
@test sd.hascoalesced(pops, ldb, 1)
@test !sd.hascoalesced(pops, ldb, 2)


sdata3 = Array(Int, 2, 2, 2)
sdata3[1, 1, 1] = 1
sdata3[2, 1, 1] = 9
sdata3[1, 2, 1] = 3
sdata3[2, 2, 1] = 11
sdata3[1, 1, 2] = 2
sdata3[2, 1, 2] = 10
sdata3[1, 2, 2] = 4
sdata3[2, 2, 2] = 12
@test sd.toarray(pops[:,1], :state) == sdata3
sdata3[1, 1, 1] = 5
sdata3[2, 1, 1] = 13
sdata3[1, 2, 1] = 7
sdata3[2, 2, 1] = 15
sdata3[1, 1, 2] = 6
sdata3[2, 1, 2] = 14
sdata3[1, 2, 2] = 8
sdata3[2, 2, 2] = 16
@test sd.toarray(pops[:,2], :state) == sdata3


const acs, gcs, hcs = sd.counts(pops[:,1])
@test acs == [Dict{Int,Int}(1=>1, 2=>1, 9=>1, 10=>1), Dict{Int,Int}(3=>1, 4=>1, 11=>1, 12=>1)]
@test gcs == [Dict{Tuple{Int,Int}, Int}((1,2)=>1, (9,10)=>1), Dict{Tuple{Int,Int}, Int}((3,4)=>1, (11, 12)=>1)]
@test hcs == Dict{Tuple{Int,Int}, Int}((1,3)=>1, (2,4)=>1, (9, 11)=>1, (10,12)=>1)

const afs, gfs, hfs = sd.spectra(pops[:,1])
@test afs == [Dict{Int,Int}(1=>4), Dict{Int,Int}(1=>4)]
@test gfs == [Dict{Int,Int}(1=>2), Dict{Int,Int}(1=>2)]
@test hfs == Dict{Int,Int}(1=>4)

@test sd.history(ldb, 1) == [0, 1]
@test sd.history(ldb, 2) == [0, 2]
@test sd.history(ldb, 3) == [0, 3]
@test sd.history(ldb, 4) == [0, 4]
@test sd.history(ldb, 5) == [0, 1, 5]
@test sd.history(ldb, 6) == [0, 2, 6]
@test sd.history(ldb, 7) == [0, 3, 7]
@test sd.history(ldb, 8) == [0, 4, 8]
@test sd.history(ldb, 9) == [0, 1, 5 ,9]
@test sd.history(ldb, 10) == [0, 1, 5, 10]
@test sd.history(ldb, 11) == [0, 2, 6, 11]
@test sd.history(ldb, 12) == [0, 2, 6, 12]
@test sd.history(ldb, 13) == [0, 3, 7, 13]
@test sd.history(ldb, 14) == [0, 3, 7, 14]
@test sd.history(ldb, 15) == [0, 4, 8, 15]
@test sd.history(ldb, 16) == [0, 4, 8, 16]

ds = sd.distances(pops[:,1], ldb)
@test ds == [1 9 10 2;
             2 13 14 2]
ds = sd.distances(pops[:,2], ldb)
@test ds == [1 9 11 6;
             1 9 12 6;
             1 11 12 2;
             2 13 14 2;
             2 13 15 6;
             2 13 16 6;
             2 14 15 6;
             2 14 16 6;
             2 15 16 2];
