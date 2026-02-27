using DataFrames
using Graphs


### predefined sets of network parameters ###

""" extract the largest connected components as a Graph object """
function main_component(g)
    c = connected_components(g)
    _, i = findmax(length.(c))
    g[c[i]]
end


### Complete Bipartite Graph ###
gen_cbg(row::DataFrameRow) = complete_bipartite_graph(row.N - row.n2, row.n2)

""" complete_bipartite_graph() with N = 256 """
function cbg256()
    df = DataFrame([:id => 0,  # to be updated
        :N => 256,
        :n2 => [1, 2, 4, 9, 20, 40, 80],
    ])
    df.id = ["cbg256-$i" for i in 1:size(df, 1)]
    return (df, gen_cbg)
end

""" complete_bipartite_graph() with N =512 """
function cbg512()
    df = DataFrame([:id => 0,  # to be updated
        :N => 512,
        :n2 => [1, 2, 4, 9, 20, 40, 80, 160],
    ])
    df.id = ["cbg512-$i" for i in 1:size(df, 1)]
    return (df, gen_cbg)
end


### Static Scale-Free graph ###
gen_ssf(row::DataFrameRow) = main_component(static_scale_free(
    row.N, row.m * row.N, row.γ; seed=row.nwseed
))

""" static_scale_free() with N nodes """
function ssfbase(N::Int)
    t_m = (3, 6)
    t_γ = (2.0, 2.4, 3.4)
    t_seed = (1, 2, 3)

    # cartesian product of three parameters
    v_m = Vector{Int64}(); v_γ = Vector{Float64}()
    v_seed = Vector{Int64}()
    for m in t_m, γ in t_γ, s in t_seed
        push!(v_m, m); push!(v_γ, γ); push!(v_seed, s)
    end

    df = DataFrame([:id => 0,  # to be updated
        :N => N,
        :m => v_m, :γ => v_γ, :nwseed => v_seed,
    ])
    df.id = ["ssf$N-$i" for i in 1:size(df, 1)]
    return (df, gen_ssf)
end

""" static_scale_free() with N = 256 """
function ssf256() return ssfbase(256) end

""" static_scale_free() with N = 512 """
function ssf512() return ssfbase(512) end


### Barabasi-Albert model ###
gen_ba(row::DataFrameRow) = main_component(barabasi_albert(
    row.N, row.k; seed=row.nwseed
))

""" barabasi_albert() with N nodes """
function babase(N::Int)
    t_k = (3, 6)
    t_seed = (1, 2, 3)

    # cartesian product of two parameters
    v_k = Vector{Int64}()
    v_seed = Vector{Int64}()
    for k in t_k, s in t_seed
        push!(v_k, k); push!(v_seed, s)
    end

    df = DataFrame([:id => 0,  # to be updated
        :N => N,
        :k => v_k, :nwseed => v_seed,
    ])
    df.id = ["ba$N-$i" for i in 1:size(df, 1)]
    return (df, gen_ba)
end

""" barabasi_albert() with N = 256 """
function ba256() return babase(256) end

""" barabasi_albert() with N = 512 """
function ba512() return babase(512) end


### Erdos-Renyi network ###
gen_er(row::DataFrameRow) = main_component(
    erdos_renyi(row.N, Int(row.N * row.m / 2); seed=row.nwseed)
)

""" erdos_renyi() with N nodes """
function erbase(N::Int)
    t_m = (4, 16)
    t_seed = (1, 2, 3)

    # cartesian product of two parameters
    v_m = Vector{Int64}()
    v_seed = Vector{Int64}()
    for m in t_m, s in t_seed
        push!(v_m, m); push!(v_seed, s)
    end

    df = DataFrame([:id => 0,  # to be updated
        :N => N,
        :m => v_m, :nwseed => v_seed,
    ])
    df.id = ["er$N-$i" for i in 1:size(df, 1)]
    return (df, gen_er)
end

""" erdos_renyi() with N = 256 """
function er256() return erbase(256) end

""" erdos_renyi() with N = 512 """
function er512() return erbase(512) end


### Random Regular Graph ###
""" random_regular_graph() with N nodes """
function gen_rrg(row::DataFrameRow)
    g = random_regular_graph(row.N, row.k; seed=row.nwseed)
    @assert is_connected(g) "RRG is not connected with N=$(row.N), k=$(row.k), seed=$(row.nwseed)"
    return g
end

function rrgbase(N::Int)
    t_k = (4, 32)
    t_seed = (1, 2, 3)

    # cartesian product of two parameters
    v_k = Vector{Int64}()
    v_seed = Vector{Int64}()
    for k in t_k, s in t_seed
        push!(v_k, k); push!(v_seed, s)
    end

    df = DataFrame([:id => 0,  # to be updated
        :N => N,
        :k => v_k, :nwseed => v_seed,
    ])
    df.id = ["rrg$N-$i" for i in 1:size(df, 1)]
    return (df, gen_rrg)
end

""" random_regular_graph() with N = 256 """
function rrg256() return rrgbase(256) end

""" random_regular_graph() with N = 512 """
function rrg512() return rrgbase(512) end

