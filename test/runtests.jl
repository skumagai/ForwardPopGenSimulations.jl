import SymmetricDominance
using Base.Test

const sd = SymmetricDominance

const gdb = sd.GeneDB()

for _ = 1:2
    sd.insert!(sd.WithNewAllele, gdb, 1, 0)
end
for i = 1:2, _ = 1:2
    sd.insert!(sd.WithNewAllele, gdb, 2, i)
end
for i = 3:6, _ = 1:2
    sd.insert!(sd.WithNewAllele, gdb, 3, i)
end
@test gdb.currentid == 14
@test gdb.currentstate == 14
sd.insert!(sd.WithoutNewAllele, gdb, 4, 1, sd.select(gdb, 1, :state)[1])
@test gdb.currentid == 15
@test gdb.currentstate == 14

@test sd.isidbystate(gdb, 1, 1) == true
@test sd.isidbystate(gdb, 1, 15) == true
@test sd.isidbystate(gdb, 15, 1) == true
@test sd.isidbystate(gdb, 1, 2) != true
@test sd.isidbystate(gdb, 2, 1) != true
@test sd.isidbystate(gdb, 2, 6) != true
@test sd.isidbystate(gdb, 6, 2) != true

@test sd.getancestor(gdb, 1, 1) == 1
@test sd.getancestor(gdb, 2, 1) == 2
@test sd.getancestor(gdb, 3, 1) == 1
@test sd.getancestor(gdb, 4, 1) == 1
@test sd.getancestor(gdb, 5, 1) == 2
@test sd.getancestor(gdb, 6, 1) == 2
@test sd.getancestor(gdb, 7, 1) == 1
@test sd.getancestor(gdb, 8, 1) == 1
@test sd.getancestor(gdb, 9, 1) == 1
@test sd.getancestor(gdb, 10, 1) == 1
@test sd.getancestor(gdb, 11, 1) == 2
@test sd.getancestor(gdb, 12, 1) == 2
@test sd.getancestor(gdb, 13, 1) == 2
@test sd.getancestor(gdb, 14, 1) == 2
@test sd.getancestor(gdb, 15, 1) == 1
@test sd.getancestor(gdb, 3, 2) == 3
@test sd.getancestor(gdb, 4, 2) == 4
@test sd.getancestor(gdb, 5, 2) == 5
@test sd.getancestor(gdb, 6, 2) == 6
@test sd.getancestor(gdb, 7, 2) == 3
@test sd.getancestor(gdb, 8, 2) == 3
@test sd.getancestor(gdb, 9, 2) == 4
@test sd.getancestor(gdb, 10, 2) == 4
@test sd.getancestor(gdb, 11, 2) == 5
@test sd.getancestor(gdb, 12, 2) == 5
@test sd.getancestor(gdb, 13, 2) == 6
@test sd.getancestor(gdb, 14, 2) == 6
@test sd.getancestor(gdb, 15, 2) == 1
@test sd.getancestor(gdb, 7, 3) == 7
@test sd.getancestor(gdb, 8, 3) == 8
@test sd.getancestor(gdb, 9, 3) == 9
@test sd.getancestor(gdb, 10, 3) == 10
@test sd.getancestor(gdb, 11, 3) == 11
@test sd.getancestor(gdb, 12, 3) == 12
@test sd.getancestor(gdb, 13, 3) == 13
@test sd.getancestor(gdb, 14, 3) == 14
@test sd.getancestor(gdb, 15, 3) == 1

const o = sd.Organism(3)
@test sd.nloci(o) == 3

const pop = sd.Population(2, 1)
pop[1, 1, 1] = 7
pop[2, 1, 1] = 8
pop[1, 1, 2] = 9
pop[2, 1, 2] = 10
@test sd.hascoalesced(gdb, pop, 1, 1)
@test !sd.hascoalesced(gdb, pop, 1, 2)

for i = 1:15
    gdb.state[i] += 10
end

sdata3 = Array(Int, 2, 1, 2)
sdata3[1, 1, 1] = 7
sdata3[2, 1, 1] = 8
sdata3[1, 1, 2] = 9
sdata3[2, 1, 2] = 10
@test sd.toarray(gdb, pop, :id) == sdata3
sdata3[1, 1, 1] = 17
sdata3[2, 1, 1] = 18
sdata3[1, 1, 2] = 19
sdata3[2, 1, 2] = 20
@test sd.toarray(gdb, pop, :state) == sdata3
sdata3[1, 1, 1] = 3
sdata3[2, 1, 1] = 3
sdata3[1, 1, 2] = 4
sdata3[2, 1, 2] = 4
@test sd.toarray(gdb, pop, :parent) == sdata3
sdata3[1, 1, 1] = 3
sdata3[2, 1, 1] = 3
sdata3[1, 1, 2] = 3
sdata3[2, 1, 2] = 3
@test sd.toarray(gdb, pop, :epoch) == sdata3

const acs, gcs, hcs = sd.counts(gdb, pop)
@test acs == [Dict{Int,Int}(17=>1, 18=>1, 19=>1, 20=>1)]
@test gcs == [Dict{Tuple{Int,Int}, Int}((17, 19)=>1, (18, 20)=>1)]
@test hcs == Dict{Tuple{Int}, Int}((17,)=>1, (18,)=>1, (19,)=>1, (20,)=>1)

const afs, gfs, hfs = sd.spectra(gdb, pop)
@test afs == [Dict{Int,Int}(1=>4)]
@test gfs == [Dict{Int,Int}(1=>2)]
@test hfs == Dict{Int,Int}(1=>4)

@test sd.history(gdb, 1, 1) == [1]
@test sd.history(gdb, 2, 1) == [2]
@test sd.history(gdb, 3, 1) == [1, 3]
@test sd.history(gdb, 4, 1) == [1, 4]
@test sd.history(gdb, 5, 1) == [2, 5]
@test sd.history(gdb, 6, 1) == [2, 6]
@test sd.history(gdb, 7, 1) == [1, 3, 7]
@test sd.history(gdb, 8, 1) == [1, 3, 8]
@test sd.history(gdb, 9, 1) == [1, 4 ,9]
@test sd.history(gdb, 10, 1) == [1, 4, 10]
@test sd.history(gdb, 11, 1) == [2, 5, 11]
@test sd.history(gdb, 12, 1) == [2, 5, 12]
@test sd.history(gdb, 13, 1) == [2, 6, 13]
@test sd.history(gdb, 14, 1) == [2, 6, 14]
@test sd.history(gdb, 15, 1) == [1, 15]
@test sd.history(gdb, 3, 2) == [3]
@test sd.history(gdb, 4, 2) == [4]
@test sd.history(gdb, 5, 2) == [5]
@test sd.history(gdb, 6, 2) == [6]
@test sd.history(gdb, 7, 2) == [3, 7]
@test sd.history(gdb, 8, 2) == [3, 8]
@test sd.history(gdb, 9, 2) == [4 ,9]
@test sd.history(gdb, 10, 2) == [4, 10]
@test sd.history(gdb, 11, 2) == [5, 11]
@test sd.history(gdb, 12, 2) == [5, 12]
@test sd.history(gdb, 13, 2) == [6, 13]
@test sd.history(gdb, 14, 2) == [6, 14]
@test sd.history(gdb, 15, 2) == [1, 15]
@test sd.history(gdb, 7, 3) == [7]
@test sd.history(gdb, 8, 3) == [8]
@test sd.history(gdb, 9, 3) == [9]
@test sd.history(gdb, 10, 3) == [10]
@test sd.history(gdb, 11, 3) == [11]
@test sd.history(gdb, 12, 3) == [12]
@test sd.history(gdb, 13, 3) == [13]
@test sd.history(gdb, 14, 3) == [14]
@test sd.history(gdb, 15, 3) == [1, 15]

ds = sd.distances(gdb, pop, 1)
@test ds == [1 17 18 2;
             1 17 19 4;
             1 17 20 4;
             1 18 19 4;
             1 18 20 4;
             1 19 20 2]
ds = sd.distances(gdb, pop, 2)
@test ds == [1 17 18 2;
             1 17 19 -1;
             1 17 20 -1;
             1 18 19 -1;
             1 18 20 -1;
             1 19 20 2]
ds = sd.distances(gdb, pop, 3)
@test ds == [1 17 18 -1;
             1 17 19 -1;
             1 17 20 -1;
             1 18 19 -1;
             1 18 20 -1;
             1 19 20 -1]

@test sd.nsegsites(gdb, pop, 1) == [6]
@test sd.nsegsites(gdb, pop, 2) == [-1]
@test sd.nsegsites(gdb, pop, 3) == [-1]

@test sd.tmrca(gdb, pop) == [1]
