
module RunUcmp
export run_ucmp

using ExaAdmm, CSV, DataFrames
# Functions for collecting solutions
function collect_ucmpmodel_solution(ucmpmodel::ExaAdmm.UCMPModel)
    T = ucmpmodel.mpmodel.len_horizon
    ngen = ucmpmodel.mpmodel.models[1].grid_data.ngen
    nbus = ucmpmodel.mpmodel.models[1].grid_data.nbus
    nline = ucmpmodel.mpmodel.models[1].grid_data.nline
    pg = zeros(T, ngen)
    qg = zeros(T, ngen)
    for t in 1:T
        for i in 1:ngen
            pg[t,i] = ucmpmodel.mpmodel.models[t].solution.u_curr[2*i-1]
            qg[t,i] = ucmpmodel.mpmodel.models[t].solution.u_curr[2*i]
        end
    end
    return pg, qg
end

function extract_generation(ucmpmodel::ExaAdmm.UCMPModel)
    T = ucmpmodel.mpmodel.len_horizon
    ngen = ucmpmodel.mpmodel.models[1].grid_data.ngen
    pg = zeros(T, ngen)
    qg = zeros(T, ngen)
    for t in 1:T
        for i in 1:ngen
            pg[t,i] = ucmpmodel.mpmodel.models[t].solution.u_curr[2*i-1]
            qg[t,i] = ucmpmodel.mpmodel.models[t].solution.u_curr[2*i]
        end
    end
    return pg, qg
end

function extract_unit_commitment(ucmpmodel::ExaAdmm.UCMPModel)::Matrix{Int64}
    # Converts vector of vectors to Matrix of Int64s
    return reduce(hcat,ucmpmodel.mpmodel.on_status)
    
end

function extract_voltages(ucmpmodel::ExaAdmm.UCMPModel)
    mod = ucmpmodel.mpmodel
    T = ucmpmodel.mpmodel.len_horizon
    nbus = mod.models[1].grid_data.nbus
    line_start = mod.models[1].line_start

    voltage_magnitude = zeros(T, nbus)
    voltage_angle = zeros(T, nbus)

    for t in 1:T
        # Extract voltage magnitudes and angles from u_curr
        u_curr = mod.models[t].solution.u_curr
        FrStart = mod.models[t].grid_data.FrStart
        FrIdx = mod.models[t].grid_data.FrIdx
        ToStart = mod.models[t].grid_data.ToStart
        ToIdx = mod.models[t].grid_data.ToIdx
        for b=1:nbus
            count = 0
            for k=FrStart[b]:FrStart[b+1]-1
                l = FrIdx[k]
                voltage_magnitude[t,b] += sqrt(u_curr[line_start + 8*(l-1)+4])
                voltage_angle[t,b] += u_curr[line_start + 8*(l-1)+6]
                count += 1
            end
            for k=ToStart[b]:ToStart[b+1]-1
                l = ToIdx[k]
                voltage_magnitude[t,b] += sqrt(u_curr[line_start + 8*(l-1)+5])
                voltage_angle[t,b] += u_curr[line_start + 8*(l-1)+7]
                count += 1
            end

            voltage_magnitude[t,b] /= count
            voltage_angle[t,b] /= count
        end
    end
    return voltage_magnitude, voltage_angle
end

function run_ucmp(data_dir, case_name)
    case_file = joinpath(data_dir, case_name * ".m")
    load_file = joinpath(data_dir, "load")
    gen_prefix = joinpath(data_dir, "gen")

    # get number of timesteps from load file
    p_load_file = joinpath(data_dir, "load.Pd")
    load_data = DataFrame(CSV.File(p_load_file))
    # CSV.read(p_load_file)
    timesteps = size(load_data, 2)
    scaling_factor = 1
    println("Number of timesteps: ", timesteps)
    env, mod = ExaAdmm.solve_ucmp(case_file, load_file, gen_prefix; ramp_ratio=0.03, rho_pq=1e3, rho_va=1e4, rho_uc=1e4, outer_iterlim=10, inner_iterlim=10, start_period=1, end_period=timesteps, scale=1e-4, tight_factor=0.99, use_gpu=false, warm_start=false)

    # Active and reactive power generation
    pg, qg = extract_generation(mod)
    # Voltage magnitudes and angles
    vm, va = extract_voltages(mod)
    # Unit commitment keys (Int64)
    uc = extract_unit_commitment(mod)


    pg_df = DataFrame(pg* scaling_factor, :auto) 
    qg_df = DataFrame(qg* scaling_factor, :auto)
    vm_df = DataFrame(vm, :auto)
    va_df = DataFrame(va, :auto)
    uc_df = DataFrame(uc, :auto)

    # Save pg, qg to csv in data_dir
    CSV.write(joinpath(data_dir, "pg.csv"), pg_df)
    CSV.write(joinpath(data_dir, "qg.csv"), qg_df)
    CSV.write(joinpath(data_dir, "vm.csv"), vm_df)
    CSV.write(joinpath(data_dir, "va.csv"), va_df)
    CSV.write(joinpath(data_dir, "uc.csv"), uc_df)
end

end