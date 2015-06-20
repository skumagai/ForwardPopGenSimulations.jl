export BasicData,
       settmax!,
       db,
       nextstate!,
       clean!,
       transmit!

type BasicData
    db::GeneDB
    state::Int
    t::Int
    tmax::Int
    BasicData() = new(GeneDB(), 0, 0)
end

@inline settmax!(d::BasicData, tmax) = d.tmax = tmax
@inline db(d::BasicData) = d.db
@inline Base.time(d::BasicData) = d.t
@inline nextstate!(d::BasicData) = d.state += 1
Base.start(d::BasicData) = d.t = 0
Base.next(d::BasicData, state) = (d.t = state + 1; (d.t, d.t))
Base.done(d::BasicData, state) = state == d.tmax

@inline clean!(d::BasicData, cmin::Int, cmax::Int) = cleandb!(db(d), cmin, cmax)
@inline transmit!(d::BasicData, t::Int, pid::Int) = transmitgene!(db(d), t, pid)
@inline transmit!(d::BasicData, t::Int, s::Int, pid::Int) = transmitgene!(db(d), t, s, pid)
