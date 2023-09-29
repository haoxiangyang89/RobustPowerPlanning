# define the type of variables
export
    costDataType,GenericData,ParsedData

struct costDataType
    model :: Int64
    upCost :: Float64
    downCost :: Float64
    n :: Int64
    params :: Array{Float64,1}
end

struct uncertainData
    DPmax :: Array{Float64,1}
    DQmax :: Array{Float64,1}
    DPmin :: Array{Float64,1}
    DQmin :: Array{Float64,1}
    DP0 :: Array{Float64,1}
    DQ0 :: Array{Float64,1}

    RESPmax :: Float64
    RESPmin :: Float64
    RESP0 :: Array{Float64,1}
end

mutable struct fixedData
    # static network data
    baseMVA :: Float64
    bType :: Dict{Int64,Any}

    IDList :: Array{Int64,1}
    genIDList :: Array{Int64,1}
    brList :: Array{Any,1}
    brRev :: Dict{Any,Any}

    Loc :: Dict{Int64,Any}
    LocRev :: Dict{Int64,Any}
    Vmax :: Dict{Int64,Any}
    Vmin :: Dict{Int64,Any}
    Pmax :: Dict{Int64,Any}
    Pmin :: Dict{Int64,Any}
    Qmax :: Dict{Int64,Any}
    Qmin :: Dict{Int64,Any}
    gs :: Dict{Int64,Any}
    bs :: Dict{Int64,Any}
    Vmag :: Dict{Int64,Any}
    Vang :: Dict{Int64,Any}
    Pd :: Dict{Int64,Any}
    Qd :: Dict{Int64,Any}
    Pg :: Dict{Int64,Any}
    Qg :: Dict{Int64,Any}
    RU :: Dict{Int64,Any}
    RD :: Dict{Int64,Any}

    g :: Dict{Tuple{Int64,Int64,Int64},Any}
    b :: Dict{Tuple{Int64,Int64,Int64},Any}
    Rdict :: Dict{Tuple{Int64,Int64,Int64},Any}
    Xdict :: Dict{Tuple{Int64,Int64,Int64},Any}
    bc :: Dict{Tuple{Int64,Int64,Int64},Any}
    θmax :: Dict{Tuple{Int64,Int64,Int64},Any}
    θmin :: Dict{Tuple{Int64,Int64,Int64},Any}
    rateA :: Dict{Tuple{Int64,Int64,Int64},Any}
    τ1 :: Dict{Tuple{Int64,Int64,Int64},Any}
    τ2 :: Dict{Tuple{Int64,Int64,Int64},Any}
    σ :: Dict{Tuple{Int64,Int64,Int64},Any}

    cp :: Dict{Any,Any}
    cq :: Dict{Any,Any}
    cz :: Any

    busInd :: Dict{Any,Any}
    branchDict1 :: Dict{Any,Any}
    branchDict2 :: Dict{Any,Any}
    connectPair :: Array{Any,1}
    connectDict :: Dict{Any,Any}
    kpDict :: Dict{Any,Any}

    Δt :: Float64

end


mutable struct batteryData
    # battery information: charging/discharging factor, capacity
    IDList :: Array{Any,1}
    Loc :: Dict{Int64,Any}
    bInv :: Dict{Int64,Any}
    ηα :: Dict{Int64,Any}
    ηβ :: Dict{Int64,Any}
    cap :: Dict{Int64,Any}
    cost :: Dict{Int64,Any}
    uCap :: Dict{Int64,Any}
end

struct demandData
  # demand data: active demand and reactive demand
  T :: Int64
  pd :: Dict{Any,Any}
  qd :: Dict{Any,Any}
end

struct windData
    # demand data: active demand and reactive demand
    hList :: Array{Int64,1}
    ph :: Dict{Any,Any}
    cost :: Dict{Any,Any}
  end