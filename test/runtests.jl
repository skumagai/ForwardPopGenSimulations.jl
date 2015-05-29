import SymmetricDominance
using Base.Test

const sd = SymmetricDominance

gdb = sd.GeneDB()

type StateCounter
    state::Int
end

nextstate!(s::StateCounter) = (s.state += 1; s.state)

sc = StateCounter(0)

for _ = 1:2
    sd.insert!(gdb, sd.GeneRecord(1, nextstate!(sc)))
end
for i = 1:2, _ = 1:2
    sd.insert!(gdb, sd.GeneRecord(2, nextstate!(sc), gdb[i]))
end
for i = 3:6, _ = 1:2
    sd.insert!(gdb, sd.GeneRecord(3, nextstate!(sc), gdb[i]))
end
@test gdb.currentid == 14
@test sc.state == 14
sd.insert!(gdb, sd.GeneRecord(4, gdb[1]))
@test gdb.currentid == 15
@test sc.state == 14
@test Set(keys(gdb)) == Set(1:15)
for i = 1:15
    @test haskey(gdb, i) == true
end
@test haskey(gdb, 16) == false
for i = 1:2
    @test isa(gdb[i].event, sd.Transmission) == true
end
for i = 3:14
    @test isa(gdb[i].event, sd.Mutation) == true
end
@test isa(gdb[15].event, sd.Transmission) == true

@test sd.isidbystate(gdb, 1, 1) == true
@test sd.isidbystate(gdb, 1, 15) == true
@test sd.isidbystate(gdb, 15, 1) == true
@test sd.isidbystate(gdb, 15, 15) == true
@test sd.isidbystate(gdb, 1, 2) == false
@test sd.isidbystate(gdb, 2, 1) == false
@test sd.isidbystate(gdb, 2, 6) == false
@test sd.isidbystate(gdb, 6, 2) == false

const o = sd.Organism(3)
@test sd.nloci(o) == 3

gids = [11, 8, 9, 10]
@test sd.mrca(gdb, gids).id == 0
gids = [7:10;]
@test sd.mrca(gdb, gids).id == 1

for i = 1:15
    gdb[i].state += 10
end

pop = sd.Population(2, 1)
pop[1, 1, 1] = 7
pop[2, 1, 1] = 8
pop[1, 1, 2] = 9
pop[2, 1, 2] = 10

sdata3 = Array(Int, 2, 1, 2)
sdata3[1, 1, 1] = 7
sdata3[2, 1, 1] = 8
sdata3[1, 1, 2] = 9
sdata3[2, 1, 2] = 10
@test sd.toarray(gdb, pop, :id) == sdata3
@test sd.toarray(gdb, gids, :id) == gids
sdata3[1, 1, 1] = 17
sdata3[2, 1, 1] = 18
sdata3[1, 1, 2] = 19
sdata3[2, 1, 2] = 20
@test sd.toarray(gdb, pop, :state) == sdata3
@test sd.toarray(gdb, gids, :state) == [17:20;]
sdata3[1, 1, 1] = 3
sdata3[2, 1, 1] = 3
sdata3[1, 1, 2] = 4
sdata3[2, 1, 2] = 4
@test sd.toarray(gdb, pop, :parent) == map(x -> gdb[x], sdata3)
@test sd.toarray(gdb, gids, :parent) == map(x -> gdb[x], [3, 3, 4, 4])
sdata3[1, 1, 1] = 3
sdata3[2, 1, 1] = 3
sdata3[1, 1, 2] = 3
sdata3[2, 1, 2] = 3
@test sd.toarray(gdb, pop, :epoch) == sdata3
@test sd.toarray(gdb, gids, :epoch) == [3, 3, 3, 3]

const acs, gcs, hcs = sd.counts(gdb, pop)
@test acs == [Dict{Int,Int}(17=>1, 18=>1, 19=>1, 20=>1)]
@test gcs == [Dict{Tuple{Int,Int}, Int}((17, 19)=>1, (18, 20)=>1)]
@test hcs == Dict{Tuple{Int}, Int}((17,)=>1, (18,)=>1, (19,)=>1, (20,)=>1)

const afs, gfs, hfs = sd.spectra(gdb, pop)
@test afs == [Dict{Int,Int}(1=>4)]
@test gfs == [Dict{Int,Int}(1=>2)]
@test hfs == Dict{Int,Int}(1=>4)

@test sd.history(gdb, 1) == [1]
@test sd.history(gdb, 2) == [2]
@test sd.history(gdb, 3) == [1, 3]
@test sd.history(gdb, 4) == [1, 4]
@test sd.history(gdb, 5) == [2, 5]
@test sd.history(gdb, 6) == [2, 6]
@test sd.history(gdb, 7) == [1, 3, 7]
@test sd.history(gdb, 8) == [1, 3, 8]
@test sd.history(gdb, 9) == [1, 4 ,9]
@test sd.history(gdb, 10) == [1, 4, 10]
@test sd.history(gdb, 11) == [2, 5, 11]
@test sd.history(gdb, 12) == [2, 5, 12]
@test sd.history(gdb, 13) == [2, 6, 13]
@test sd.history(gdb, 14) == [2, 6, 14]
@test sd.history(gdb, 15) == [1, 15]

ds = sd.distances(gdb, gids)
@test ds == [0 1 2 2;
             1 0 2 2;
             2 2 0 1;
             2 2 1 0]

@test sd.nsegsites(gdb, gids) == 6

@test sd.mrca(gdb, gids).epoch == 1

gdb = sd.GeneDB()
sc = StateCounter(0)
sd.insert!(gdb, sd.GeneRecord(1, nextstate!(sc)))
sd.insert!(gdb, sd.GeneRecord(2, 1, gdb[1]))
sd.insert!(gdb, sd.GeneRecord(3, 1, gdb[2]))
sd.insert!(gdb, sd.GeneRecord(3, 1, gdb[2]))
@test Set(collect(keys(gdb))) == Set([1, 2, 3, 4])
sd.clean!(gdb, [3, 4])
@test Set(collect(keys(gdb))) == Set([2, 3, 4])

gdb = sd.GeneDB()
sd.insert!(gdb, sd.GeneRecord(1, nextstate!(sc)))
sd.insert!(gdb, sd.GeneRecord(2, 1, gdb[1]))
sd.insert!(gdb, sd.GeneRecord(2, 1, gdb[1]))
sd.insert!(gdb, sd.GeneRecord(3, 1, gdb[2]))
sd.insert!(gdb, sd.GeneRecord(3, 1, gdb[2]))
sd.clean!(gdb, [4, 5])
@test Set(collect(keys(gdb))) == Set([2, 4, 5])

gdb = sd.GeneDB()
sd.insert!(gdb, sd.GeneRecord(1, 1))
sd.insert!(gdb, sd.GeneRecord(2, 2, gdb[1]))
sd.insert!(gdb, sd.GeneRecord(3, gdb[2]))
sd.insert!(gdb, sd.GeneRecord(4, gdb[3]))
sd.insert!(gdb, sd.GeneRecord(5, 3, gdb[4]))
sd.insert!(gdb, sd.GeneRecord(5, gdb[3]))
@test Set(collect(keys(gdb))) == Set([1:6;])
@test sd.mrca(gdb, [5, 6]) == gdb[3]
sd.clean!(gdb, [5, 6])
@test Set(collect(keys(gdb))) == Set([3, 5, 6])
