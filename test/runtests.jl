import ForwardPopGenSimulations
using Base.Test

const fpgs = ForwardPopGenSimulations

gdb = fpgs.GeneDB()

type StateCounter
    state::Int
end

nextstate!(s::StateCounter) = (s.state += 1; s.state)

sc = StateCounter(0)

for _ = 1:2
    fpgs.insert!(gdb, fpgs.GeneRecord(1, nextstate!(sc)))
end
for i = 1:2, _ = 1:2
    fpgs.insert!(gdb, fpgs.GeneRecord(2, nextstate!(sc), gdb[i]))
end
for i = 3:6, _ = 1:2
    fpgs.insert!(gdb, fpgs.GeneRecord(3, nextstate!(sc), gdb[i]))
end
@test gdb.currentid == 14
@test sc.state == 14
fpgs.insert!(gdb, fpgs.GeneRecord(4, gdb[1]))
@test gdb.currentid == 15
@test sc.state == 14
@test Set(keys(gdb)) == Set(1:15)
for i = 1:15
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

@test fpgs.isidbystate(gdb, 1, 1) == true
@test fpgs.isidbystate(gdb, 1, 15) == true
@test fpgs.isidbystate(gdb, 15, 1) == true
@test fpgs.isidbystate(gdb, 15, 15) == true
@test fpgs.isidbystate(gdb, 1, 2) == false
@test fpgs.isidbystate(gdb, 2, 1) == false
@test fpgs.isidbystate(gdb, 2, 6) == false
@test fpgs.isidbystate(gdb, 6, 2) == false

gids = [11, 8, 9, 10]
@test fpgs.mrca(gdb, gids).id == 0
gids = [7:10;]
@test fpgs.mrca(gdb, gids).id == 1

for i = 1:15
    gdb[i].state += 10
end

@test fpgs.toarray(gdb, gids, :id) == gids
@test fpgs.toarray(gdb, gids, :state) == [17:20;]
@test fpgs.toarray(gdb, gids, :parent) == map(x -> gdb[x], [3, 3, 4, 4])
@test fpgs.toarray(gdb, gids, :epoch) == [3, 3, 3, 3]

@test fpgs.history(gdb, 1) == [1]
@test fpgs.history(gdb, 2) == [2]
@test fpgs.history(gdb, 3) == [1, 3]
@test fpgs.history(gdb, 4) == [1, 4]
@test fpgs.history(gdb, 5) == [2, 5]
@test fpgs.history(gdb, 6) == [2, 6]
@test fpgs.history(gdb, 7) == [1, 3, 7]
@test fpgs.history(gdb, 8) == [1, 3, 8]
@test fpgs.history(gdb, 9) == [1, 4 ,9]
@test fpgs.history(gdb, 10) == [1, 4, 10]
@test fpgs.history(gdb, 11) == [2, 5, 11]
@test fpgs.history(gdb, 12) == [2, 5, 12]
@test fpgs.history(gdb, 13) == [2, 6, 13]
@test fpgs.history(gdb, 14) == [2, 6, 14]
@test fpgs.history(gdb, 15) == [1, 15]

ds = fpgs.distances(gdb, gids)
@test ds == [0 1 2 2;
             1 0 2 2;
             2 2 0 1;
             2 2 1 0]

@test fpgs.nsegsites(gdb, gids) == 6

@test fpgs.mrca(gdb, gids).epoch == 1

gdb = fpgs.GeneDB()
sc = StateCounter(0)
fpgs.insert!(gdb, fpgs.GeneRecord(1, nextstate!(sc)))
fpgs.insert!(gdb, fpgs.GeneRecord(2, 1, gdb[1]))
fpgs.insert!(gdb, fpgs.GeneRecord(3, 1, gdb[2]))
fpgs.insert!(gdb, fpgs.GeneRecord(3, 1, gdb[2]))
@test Set(collect(keys(gdb))) == Set([1, 2, 3, 4])
fpgs.clean!(gdb, 3, 4)
@test Set(collect(keys(gdb))) == Set([2, 3, 4])

gdb = fpgs.GeneDB()
fpgs.insert!(gdb, fpgs.GeneRecord(1, nextstate!(sc)))
fpgs.insert!(gdb, fpgs.GeneRecord(2, 1, gdb[1]))
fpgs.insert!(gdb, fpgs.GeneRecord(2, 1, gdb[1]))
fpgs.insert!(gdb, fpgs.GeneRecord(3, 1, gdb[2]))
fpgs.insert!(gdb, fpgs.GeneRecord(3, 1, gdb[2]))
fpgs.clean!(gdb, 4, 5)
@test Set(collect(keys(gdb))) == Set([2, 4, 5])

gdb = fpgs.GeneDB()
fpgs.insert!(gdb, fpgs.GeneRecord(1, 1))
fpgs.insert!(gdb, fpgs.GeneRecord(2, 2, gdb[1]))
fpgs.insert!(gdb, fpgs.GeneRecord(3, gdb[2]))
fpgs.insert!(gdb, fpgs.GeneRecord(4, gdb[3]))
fpgs.insert!(gdb, fpgs.GeneRecord(5, 3, gdb[4]))
fpgs.insert!(gdb, fpgs.GeneRecord(5, gdb[3]))
@test Set(collect(keys(gdb))) == Set([1:6;])
@test fpgs.mrca(gdb, [5, 6]) == gdb[3]
fpgs.clean!(gdb, 5, 6)
@test Set(collect(keys(gdb))) == Set([3, 5, 6])

gdb = fpgs.GeneDB()
fpgs.insert!(gdb, fpgs.GeneRecord(1, 1))
fpgs.insert!(gdb, fpgs.GeneRecord(2, 2, gdb[1]))
fpgs.insert!(gdb, fpgs.GeneRecord(3, gdb[2]))
fpgs.insert!(gdb, fpgs.GeneRecord(4, gdb[3]))

fpgs.insert!(gdb, fpgs.GeneRecord(1, 11))
fpgs.insert!(gdb, fpgs.GeneRecord(2, 2, gdb[5]))
fpgs.insert!(gdb, fpgs.GeneRecord(3, gdb[6]))
fpgs.insert!(gdb, fpgs.GeneRecord(4, gdb[7]))

fpgs.insert!(gdb, fpgs.GeneRecord(5, 3, gdb[4]))
fpgs.insert!(gdb, fpgs.GeneRecord(5, gdb[3]))

fpgs.insert!(gdb, fpgs.GeneRecord(5, 3, gdb[8]))
fpgs.insert!(gdb, fpgs.GeneRecord(5, gdb[7]))
@test Set(collect(keys(gdb))) == Set([1:12;])
@test fpgs.mrca(gdb, [9, 10]) == gdb[3]
@test fpgs.mrca(gdb, [11, 12]) == gdb[7]
fpgs.clean!(gdb, 9, 12)
@test Set(collect(keys(gdb))) == Set([3, 7, 9, 10, 11, 12])

pars = Array{Int}(1)
fpgs.selectparents!(pars, 10, replace=false)
@test 0 < pars[1] < 11
fpgs.selectparents!(pars, 10, replace=true)
@test 0 < pars[1] < 11
pars = Array{Int}(2)
vec = Array{Bool}(100)
for i = 1:100
    fpgs.selectparents!(pars, 2, replace=false)
    vec[i] = pars[1] == pars[2]
end
@test any(vec) == false

for i = 1:100
    fpgs.selectparents!(pars, 2, replace=true)
    vec[i] = pars[1] == pars[2]
end
@test any(vec) == true

mutarr = Array{Bool}(3, 2)
rates= zeros(3, 2)
fpgs.selectmutatedsites!(mutarr, rates)
@test all(mutarr .== false) == true
rates[:,:] = 1.0
fpgs.selectmutatedsites!(mutarr, rates)
@test all(mutarr .== true) == true
rates[:,1] = 0.0
fpgs.selectmutatedsites!(mutarr, rates)
@test all(mutarr[:,1] .== false) == true
@test all(mutarr[:,2] .== true) == true
