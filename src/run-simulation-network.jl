using Distributed
using Dates, TimeZones, Logging, DelimitedFiles, CSV
using DataFrames, Random, Statistics
using DifferentialEquations, Graphs, NetworkDynamics


####################
# set up parameters
####################

### prepare for file generation ###
# get current time
filenamebase = Dates.format(now(tz"Asia/Tokyo"), "yymmdd-HHMMSS")
# helper method for logging
getstrnow() = Dates.format(now(tz"Asia/Tokyo"), "Y-mm-dd HH:MM:SS.s")

# change directory if needed
const wd = "src"
if splitdir(pwd())[end] != wd cd(wd) end
# -> now, we are inside "src" directory!

# make sure that output directory exists
const outdir = "output"
mkpath(outdir)

### logging ###
logpath = joinpath(outdir, filenamebase * "_log.txt")
@info "stdout and stderr will be redirected to $logpath"
logio = open(logpath, "a")

redirect_stdio(; stdout=logio, stderr=logio)

""" flush IOStream for logging from time to time """
function flushlogio(io, c, interval=2.0)
    @info "flushlogio task started"
    while true
        if isready(c)
            # if anything is put in the channel `c`, break the loop
            break
        end
        # flush every interval seconds
        flush(io)
        sleep(interval)  # make room for other tasks by sleeping
    end
    return "flushlogio terminated"
end

# start flush task: make sure at least two threads are available!
flushchannel = Channel(1)
flushtask = Threads.@spawn flushlogio(logio, flushchannel)  # never use @async !!!

### parameters ###
# fixed dynamical parameters
r, D, ξ = (0.05, 0.005, 0.5)

# read command line input
if length(ARGS) != 5
    error("5 args must be passed: $(length(ARGS)) were given")
end
nwparamset::String = ARGS[1]
K::Float64 = parse(Float64, ARGS[2])
firstseed::Int64, numsample::Int64 = parse.(Int, ARGS[3:4])
numworker::Int64 = parse(Int, ARGS[5])

# generate (dynamical) seeds
seeds = collect(firstseed : firstseed + numsample - 1)

# get network parameters `df` and define network generator `nwgen()`
include("nwparams.jl")  # load predefined sets of network parameters

df_p, nwgen = if nwparamset == "cbg256" cbg256()
elseif nwparamset == "cbg512" cbg512()
elseif nwparamset == "ssf256" ssf256()
elseif nwparamset == "ssf512" ssf512()
elseif nwparamset == "ba256" ba256()
elseif nwparamset == "ba512" ba512()
elseif nwparamset == "er256" er256()
elseif nwparamset == "er512" er512()
elseif nwparamset == "rrg256" rrg256()
elseif nwparamset == "rrg512" rrg512()
else 
    throw(ArgumentError("unknown nwparamset: $nwparamset")) 
end

# add other parameters
df_p.K .= K
df_p.r .= r
df_p.D .= D
df_p.ξ .= ξ

# save parameter values
CSV.write(joinpath(outdir, filenamebase * "_params.csv"), df_p)
@info "param file generated"

##################################
# prepare for running simulations
##################################

### add workers ###
addprocs(numworker)
@info "nprocs: $(nprocs())"

### preparation in all processes ###
# -----
@everywhere begin

### load modules in workers ###
using Statistics
using DifferentialEquations, Graphs, NetworkDynamics

### define edge and vertext dynamics ###
function diffusionedge_g!(e_dst, v_src, v_dst, p, t)
    # e_dst, v_src, v_dst are arrays, so we use the broadcasting operator .
    @. e_dst = v_src - v_dst
    nothing
end

nd_diffusion_edge = EdgeModel(; 
    g=AntiSymmetric(diffusionedge_g!), outsym=[:flow]
)

bistable(x, r=$r) = -x * (x - r) * (x - 1)

function bistablediffvertex_f!(dv, v, esum, (K_deg,), t)
    # dv, v and esum are arrays, so we use the broadcasting operator .
    @. dv = bistable(v) + K_deg * esum
    #@. dv = -v * (v - r) * (v - 1) + K_deg * esum
    nothing
end

nd_bistablediff_vertex = VertexModel(; 
    f=bistablediffvertex_f!, g=StateMask(1:1), dim=1, psym=[:K_deg],
)

function noisecoef!(du, u, p, t, σ=sqrt(2 * $D))
    du .= σ
end

### prepare to define EnsembleProblem ###
""" extract escape times from simulation result """
function sol2et(sol, ξ=$ξ)
    et = zeros(size(sol.u[1]))
    for (t, u) in zip(sol.t, sol.u)
        escaped_before = et .> 0  # if escaped before this time step
        escape_now = u .>= ξ
        et[@. (!escaped_before) && escape_now] .= t
    end
    @assert all(et .> 0)
    return et
end

# give a particular seed for each trajectory
prob_func(prob, i, repeat, seeds=$seeds) = remake(prob; seed=seeds[i])

# output reduction: record the mean escape time & retcode
output_func(sol, i) = ((sol.retcode, mean(sol2et(sol))), false)

### define CallBacks for recording and termination ###
# condition for escape time measurement
function escape_cond(out, u, t, int, ξ=$ξ)
    # up-crossing corresponds to an escape
    @. out = u - ξ
end

# function to record escape orders in p: do nothing
function escape_affect!(int, event_index) end

### define methods for termination ###
# function of terminate condition
# - if all nodes are above threshold, terminate the simulation
terminate_cond(u, t, int, ξ=$ξ) = minimum(u) >= ξ

# callback to terminate the simulation
terminate_cb = DiscreteCallback(
    terminate_cond, integrator -> terminate!(integrator); 
    save_positions=(true, false)  # only save the data just before termination
)

end
# -----

""" define and solve EnsembleProblem on the given Graph """
function solve_ensembleprob(nw; 
    K=K, r=r, D=D, ξ=ξ, numsample=numsample,
    tmax=10000.0, maxiters=1e11, abstol=5e-4, reltol=1e-3,
)
    N = nv(nw)
    aet = Vector{Float64}(undef, numsample)

    # generate CallbackSet on every worker
    @everywhere begin
        # Callback to monitor escapes
        escape_cb = VectorContinuousCallback(
            escape_cond, escape_affect!, nothing, $N;
            # only save the data just before an escape
            save_positions=(true, false)
        )
        # set of Callbacks
        cb_set = CallbackSet(escape_cb, terminate_cb)
    end

    ### NetworkDynamics ###
    nd = Network(nw, nd_bistablediff_vertex, nd_diffusion_edge)
    p = NWParameter(nd)
    p.v[:, :K_deg] .= K ./ degree(nw)  # set coupling weights

    ### define the base SDEProblem ###
    u0 = zeros(N)
    sde_prob = SDEProblem(nd, noisecoef!, u0, (0.0, tmax), pflat(p))

    ### define and solve EnsembleProblem
    ensemble_prob = EnsembleProblem(sde_prob, 
        prob_func=prob_func, output_func=output_func,
    )
    sim = solve(ensemble_prob, SOSRA(), EnsembleDistributed();
        abstol=abstol, reltol=reltol, maxiters=maxiters, 
        save_everystep=false, callback=cb_set, 
        trajectories=numsample, pmap_batch_size=1,
    )

    ### check if all nodes escaped & record the results ###
    for (i, out) in enumerate(sim)
        if out[1] === ReturnCode.T(2)  # ReturnCode.Terminated = 2
            aet[i] = out[2]
        else
            @warn "Calculation was not properly terminated!" K=K N=N sample=i retcode=out[1]
            aet[i] = -1  # meaningless value
        end
    end

    return aet
end

degree2κ(degree) = mean(degree.^2) / mean(degree)^2

##################
# run simulations
##################

### prepare output csv ###
respath = joinpath(outdir, filenamebase * "_aet.csv")
# write header
open(respath, "w") do io
    writedlm(io, hcat(["id" "κ";], permutedims(["$s" for s in seeds])), ',')
end

### run simulations ###
# start calculations
for row in eachrow(df_p)
    @info "$(getstrnow())\t start working on $(row.id)"

    nw = nwgen(row)  # generate network
    aet = solve_ensembleprob(nw)  # define & solve EnsembleProblem
    # calculate degree heterogeneity
    csvrow = [[row.id]; [degree2κ(degree(nw))]; aet]
    # append result
    open(respath, "a") do io
        writedlm(io, permutedims(csvrow), ',')
    end
end

### closing ###
put!(flushchannel, :done)  # kill flush task
println(fetch(flushtask))  # wait for the task to close
@info "$(getstrnow()) completed"
close(logio)

