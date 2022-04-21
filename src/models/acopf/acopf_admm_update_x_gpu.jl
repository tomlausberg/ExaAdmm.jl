function acopf_admm_update_x_gen(
    env::AdmmEnv{Float64,CuArray{Float64,1},CuArray{Int,1},CuArray{Float64,2}},
    mod::ModelAcopf{Float64,CuArray{Float64,1},CuArray{Int,1},CuArray{Float64,2}},
    gen_solution::EmptyGeneratorSolution{Float64,CuArray{Float64,1}}
)
    sol, info = mod.solution, mod.info
    time_gen = generator_kernel_two_level(mod, mod.baseMVA, sol.u_curr, sol.v_curr, sol.z_curr, sol.l_curr, sol.rho)
    info.user.time_generators += time_gen.time
    info.time_x_update += time_gen.time
end

function acopf_admm_update_x_line(
    env::AdmmEnv{Float64,CuArray{Float64,1},CuArray{Int,1},CuArray{Float64,2}},
    mod::ModelAcopf{Float64,CuArray{Float64,1},CuArray{Int,1},CuArray{Float64,2}}
)
    par, sol, info = env.params, mod.solution, mod.info
    shmem_size = env.params.shmem_size

#=
    tmp = mod.nline
    mod.nline = 1
    par.shift_lines=3990
    @printf("GPU x_line_update\n")
=#
    if env.use_linelimit
        time_br = CUDA.@timed @cuda threads=32 blocks=mod.nline shmem=shmem_size auglag_linelimit_two_level_alternative(
                                            mod.n, mod.nline, mod.line_start,
                                            info.inner, par.max_auglag, par.mu_max, par.scale,
                                            sol.u_curr, sol.v_curr, sol.z_curr, sol.l_curr, sol.rho,
                                            par.shift_lines, mod.membuf, mod.YffR, mod.YffI, mod.YftR, mod.YftI,
                                            mod.YttR, mod.YttI, mod.YtfR, mod.YtfI,
                                            mod.FrVmBound, mod.ToVmBound, mod.FrVaBound, mod.ToVaBound)
    else
        time_br = CUDA.@timed @cuda threads=32 blocks=mode.nline shmem=shmem_size polar_kernel_two_level_alternative(mod.n, mod.nline, mod.line_start, par.scale,
                                                        sol.u_curr, sol.v_curr, sol.z_curr, sol.l_curr, sol.rho,
                                                        par.shift_lines, mod.membuf, mod.YffR, mod.YffI, mod.YftR, mod.YftI,
                                                        mod.YttR, mod.YttI, mod.YtfR, mod.YtfI, mod.FrVmBound, mod.ToVmBound)
    end
#=
    @printf("GPU x_line_update DONE\n")
    mod.nline = tmp
=#
    info.time_x_update += time_br.time
    info.user.time_branches += time_br.time
    return
end

function admm_update_x(
    env::AdmmEnv{Float64,CuArray{Float64,1},CuArray{Int,1},CuArray{Float64,2}},
    mod::ModelAcopf{Float64,CuArray{Float64,1},CuArray{Int,1},CuArray{Float64,2}}
)
    acopf_admm_update_x_gen(env, mod, mod.gen_solution)
    acopf_admm_update_x_line(env, mod)
    return
end