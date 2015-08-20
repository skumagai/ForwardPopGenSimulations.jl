using ForwardPopGenSimulations
using Base.Test

const fpgs = ForwardPopGenSimulations

core = BasicData()
gdb = db(core)

for _ = 1:2
    insert!(gdb, GeneRecord(1, nextstate!(core)))
end
for i = 1:2, _ = 1:2
    insert!(gdb, GeneRecord(2, nextstate!(core), gdb[i]))
end
for i = 3:6, _ = 1:2
    insert!(gdb, GeneRecord(3, nextstate!(core), gdb[i]))
end
@test gdb.currentid == 14
@test core.state == 14
insert!(gdb, GeneRecord(4, gdb[1]))
@test gdb.currentid == 15
@test core.state == 14
@test Set(keys(gdb)) == Set(0:15)
for i = 0:15
    @test haskey(gdb, i) == true
end
@test haskey(gdb, 16) == false
for i = 1:2
    @test isa(gdb[i].event, fpgs.Transmission) == true
end
for i = 3:14
    @test isa(gdb[i].event, fpgs.Mutation) == true
end
@test isa(gdb[15].event, fpgs.Transmission) == true

@test isidbystate(gdb, 1, 1) == true
@test isidbystate(gdb, 1, 15) == true
@test isidbystate(gdb, 15, 1) == true
@test isidbystate(gdb, 15, 15) == true
@test isidbystate(gdb, 1, 2) == false
@test isidbystate(gdb, 2, 1) == false
@test isidbystate(gdb, 2, 6) == false
@test isidbystate(gdb, 6, 2) == false

gids = [11, 8, 9, 10]
@test mrca(gdb, gids).id == 0
gids = [7:10;]
@test mrca(gdb, gids).id == 1

for i = 1:15
    gdb[i].state += 10
end

@test toarray(gdb, gids, :id) == gids
@test toarray(gdb, gids, :state) == [17:20;]
@test toarray(gdb, gids, :parent) == map(x -> gdb[x], [3, 3, 4, 4])
@test toarray(gdb, gids, :epoch) == [3, 3, 3, 3]

@test history(gdb, 1) == [1]
@test history(gdb, 2) == [2]
@test history(gdb, 3) == [1, 3]
@test history(gdb, 4) == [1, 4]
@test history(gdb, 5) == [2, 5]
@test history(gdb, 6) == [2, 6]
@test history(gdb, 7) == [1, 3, 7]
@test history(gdb, 8) == [1, 3, 8]
@test history(gdb, 9) == [1, 4 ,9]
@test history(gdb, 10) == [1, 4, 10]
@test history(gdb, 11) == [2, 5, 11]
@test history(gdb, 12) == [2, 5, 12]
@test history(gdb, 13) == [2, 6, 13]
@test history(gdb, 14) == [2, 6, 14]
@test history(gdb, 15) == [1, 15]

ds = distances(gdb, gids)
@test ds == [0 1 2 2;
             1 0 2 2;
             2 2 0 1;
             2 2 1 0]

@test fpgs.nsegsites(gdb, gids) == 6

@test mrca(gdb, gids).epoch == 1

core = BasicData()
gdb = db(core)
insert!(gdb, GeneRecord(1, nextstate!(core)))
insert!(gdb, GeneRecord(2, 1, gdb[1]))
insert!(gdb, GeneRecord(3, 1, gdb[2]))
insert!(gdb, GeneRecord(3, 1, gdb[2]))
@test Set(collect(keys(gdb))) == Set([0, 1, 2, 3, 4])
clean!(gdb, 3, 4)
@test Set(collect(keys(gdb))) == Set([0, 2, 3, 4])

gdb = GeneDB()
insert!(gdb, GeneRecord(1, nextstate!(core)))
insert!(gdb, GeneRecord(2, 1, gdb[1]))
insert!(gdb, GeneRecord(2, 1, gdb[1]))
insert!(gdb, GeneRecord(3, 1, gdb[2]))
insert!(gdb, GeneRecord(3, 1, gdb[2]))
clean!(gdb, 4, 5)
@test Set(collect(keys(gdb))) == Set([0, 2, 4, 5])

gdb = GeneDB()
insert!(gdb, GeneRecord(1, 1))
insert!(gdb, GeneRecord(2, 2, gdb[1]))
insert!(gdb, GeneRecord(3, gdb[2]))
insert!(gdb, GeneRecord(4, gdb[3]))
insert!(gdb, GeneRecord(5, 3, gdb[4]))
insert!(gdb, GeneRecord(5, gdb[3]))
@test Set(collect(keys(gdb))) == Set([0:6;])
@test mrca(gdb, [5, 6]) == gdb[3]
clean!(gdb, 5, 6)
@test Set(collect(keys(gdb))) == Set([0, 3, 5, 6])

gdb = GeneDB()
insert!(gdb, GeneRecord(1, 1))
insert!(gdb, GeneRecord(2, 2, gdb[1]))
insert!(gdb, GeneRecord(3, gdb[2]))
insert!(gdb, GeneRecord(4, gdb[3]))

insert!(gdb, GeneRecord(1, 11))
insert!(gdb, GeneRecord(2, 2, gdb[5]))
insert!(gdb, GeneRecord(3, gdb[6]))
insert!(gdb, GeneRecord(4, gdb[7]))

insert!(gdb, GeneRecord(5, 3, gdb[4]))
insert!(gdb, GeneRecord(5, gdb[3]))

insert!(gdb, GeneRecord(5, 3, gdb[8]))
insert!(gdb, GeneRecord(5, gdb[7]))
@test Set(collect(keys(gdb))) == Set([0:12;])
@test mrca(gdb, [9, 10]) == gdb[3]
@test mrca(gdb, [11, 12]) == gdb[7]
clean!(gdb, 9, 12)
@test Set(collect(keys(gdb))) == Set([0, 3, 7, 9, 10, 11, 12])

pars = Array{Int}(1)
selectparents!(pars, 10, replace=false)
@test 0 < pars[1] < 11
selectparents!(pars, 10, replace=true)
@test 0 < pars[1] < 11
pars = Array{Int}(2)
vec = Array{Bool}(100)
for i = 1:100
    selectparents!(pars, 2, replace=false)
    vec[i] = pars[1] == pars[2]
end
@test any(vec) == false

for i = 1:100
    selectparents!(pars, 2, replace=true)
    vec[i] = pars[1] == pars[2]
end
@test any(vec) == true

mutarr = Array{Bool}(3, 2)
rates= zeros(3, 2)
selectmutatedsites!(mutarr, rates)
@test all(mutarr .== false) == true
rates[:,:] = 1.0
selectmutatedsites!(mutarr, rates)
@test all(mutarr .== true) == true
rates[:,1] = 0.0
selectmutatedsites!(mutarr, rates)
@test all(mutarr[:,1] .== false) == true
@test all(mutarr[:,2] .== true) == true

@test value(Female) == 1
@test value(Male) == 2
@test value(Mother) == value(Daughter) == value(Female)
@test value(Father) == value(Son) == value(Male)

@test value(Autosome) == 1
@test value(XChromosome) == 2
@test value(YChromosome) == 3
@test value(Mitochondrion) == 4

core = BasicData()
for i = 1:10
    @test nextstate!(core) == i
    @test core.state == i
end
settmax!(core, 10)
@test core.tmax == 10
gdb = db(core)
@test isa(gdb, GeneDB)
@test time(core) == 0
for (i, t) in enumerate(core)
    @test i == t
end
@test time(core) == 10

core = BasicData()
gdb = db(core)
transmit!(gdb, 1, 0, state=1)
transmit!(gdb, 2, 1, state=2, from=1, to=2)
transmit!(gdb, 3, 3)
transmit!(gdb, 4, 4, state=3)
@test haskey(gdb, 2) == true
@test haskey(gdb, 3) == true
@test gdb[2].event == fpgs.Migration(1, 2)
@test gdb[2].state == 1
@test gdb[3].event == fpgs.Mutation()
@test gdb[3].state == 2
@test gdb[4].event == fpgs.Transmission()
@test gdb[4].state == 2
@test gdb[5].event == fpgs.Mutation()
@test gdb[5].state == 3
