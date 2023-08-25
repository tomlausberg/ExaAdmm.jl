function ucmp_auglag_generator_kernel(
    t::Int, n::Int, ngen::Int, gen_start::Int,
    major_iter::Int, max_auglag::Int, xi_max::Float64, scale::Float64,
    u::CuDeviceArray{Float64,1}, v::CuDeviceArray{Float64,1}, z::CuDeviceArray{Float64,1},
    l::CuDeviceArray{Float64,1}, rho::CuDeviceArray{Float64,1},
    r_u::CuDeviceArray{Float64,1}, r_v::CuDeviceArray{Float64,1}, r_z::CuDeviceArray{Float64,1},
    r_l::CuDeviceArray{Float64,1}, r_rho::CuDeviceArray{Float64,1}, r_s::CuDeviceArray{Float64,1},
    vr_u::CuDeviceArray{Float64,1}, vr_v::CuDeviceArray{Float64,1}, vr_z::CuDeviceArray{Float64,1},
    vr_l::CuDeviceArray{Float64,1}, vr_rho::CuDeviceArray{Float64,1}, vr_s::CuDeviceArray{Float64,1},
    uc_u::CuDeviceArray{Float64,2}, uc_v::CuDeviceArray{Float64,2}, uc_z::CuDeviceArray{Float64,2},
    uc_l::CuDeviceArray{Float64,2}, uc_rho::CuDeviceArray{Float64,2}, uc_s::CuDeviceArray{Float64,2},
    param::CuDeviceArray{Float64,2},
    pgmin::CuDeviceArray{Float64,1}, pgmax::CuDeviceArray{Float64,1},
    qgmin::CuDeviceArray{Float64,1}, qgmax::CuDeviceArray{Float64,1},
    ramp_limit::CuDeviceArray{Float64,1},
    _c2::CuDeviceArray{Float64,1}, _c1::CuDeviceArray{Float64,1}, _c0::CuDeviceArray{Float64,1}, baseMVA::Float64
)
    tx = threadIdx().x
    I = blockIdx().x

    x = CuDynamicSharedArray(Float64, n)
    xl = CuDynamicSharedArray(Float64, n, n*sizeof(Float64))
    xu = CuDynamicSharedArray(Float64, n, (2*n)*sizeof(Float64))

    @inbounds begin
        pg_idx = gen_start + 2*(I-1)
        qg_idx = gen_start + 2*(I-1) + 1
        v_idx = 3*(t-1) + 1
        w_idx = 3*(t-1) + 2
        y_idx = 3*t

        c2 = _c2[I]; c1 = _c1[I]; c0 = _c0[I]

		xl[1] = pgmin[I]
		xl[3] = pgmin[I]
		xl[2] = qgmin[I]
        xl[4] = 0.0
        xl[5] = 0.0
        xl[6] = 0.0
        xl[7] = 0.0
		xl[8] = -4*ramp_limit[I]
		xl[9] = -4*ramp_limit[I]
		xl[10] = -abs(pgmax[I])-abs(pgmin[I])
		xl[11] = -abs(pgmax[I])-abs(pgmin[I])
		xl[12] = -abs(qgmax[I])-abs(qgmin[I])
		xl[13] = -abs(qgmax[I])-abs(qgmin[I])
#=
        xl[1] = xl[3] = pgmin[I]
        xl[2] = qgmin[I]
        xl[8] = xl[9] = -4*ramp_limit[I]
        xl[10] = xl[11] = -abs(pgmax[I])-abs(pgmin[I])
        xl[12] = xl[13] = -abs(qgmax[I])-abs(qgmin[I])
=#
        # COMMENTED: THESE BOUNDS ARE TOO TIGHT FOR S VARIABLES. LEAD TO DIVERGENCE FOR 9-BUS CASE.
        # xl[8] = xl[9] = -2*ramp_limit[I]
        # xl[10] = xl[11] = -pgmax[I]
        # xl[12] = xl[13] = -qgmax[I]
		xu[1] = pgmax[I]
		xu[3] = pgmax[I]
		xu[2] = qgmax[I]
		xu[4] = 1
		xu[5] = 1
		xu[6] = 1
		xu[7] = 1
		xu[8] = 0
		xu[9] = 0
		xu[10] = 0
		xu[11] = 0
		xu[12] = 0
		xu[13] = 0
#=
        xu[1] = xu[3] = pgmax[I]
        xu[2] = qgmax[I]
        xu[4] = xu[5] = xu[6] = xu[7] = 1
        xu[8] = xu[9] = xu[10] = xu[11] = xu[12] = xu[13] = 0
=#
        CUDA.sync_threads()

        x[1] = min(xu[1], max(xl[1], u[pg_idx]))
        x[2] = min(xu[2], max(xl[2], u[qg_idx]))
        x[3] = min(xu[3], max(xl[3], r_u[I]))
        x[4] = uc_u[I, v_idx]
        x[5] = uc_u[I, w_idx]
        x[6] = uc_u[I, y_idx]
        t > 1 ? x[7] = vr_u[I] : x[7] = 0.
        x[8] = min(xu[8], max(xl[8], r_s[2*I-1]))
        x[9] = min(xu[9], max(xl[9], r_s[2*I]))
        x[10] = min(xu[10], max(xl[10], uc_s[I, 4*t-3]))
        x[11] = min(xu[11], max(xl[11], uc_s[I, 4*t-2]))
        x[12] = min(xu[12], max(xl[12], uc_s[I, 4*t-1]))
        x[13] = min(xu[13], max(xl[13], uc_s[I, 4*t]))

        param[1,I] = l[pg_idx]
        param[2,I] = l[qg_idx]
        param[3,I] = r_l[I]
        param[4,I] = uc_l[I, v_idx]
        param[5,I] = uc_l[I, w_idx]
        param[6,I] = uc_l[I, y_idx]
        param[7,I] = vr_l[I]
        # param[8,I] = r_l[4*I-1]
        # param[9,I] = r_l[4*I]
        # param[10,I] = uc_l[I, 7*t-3]
        # param[11,I] = uc_l[I, 7*t-2]
        # param[12,I] = uc_l[I, 7*t-1]
        # param[13,I] = uc_l[I, 7*t]
        param[14,I] = rho[pg_idx]
        param[15,I] = rho[qg_idx]
        param[16,I] = r_rho[I]
        param[17,I] = uc_rho[I, v_idx]
        param[18,I] = uc_rho[I, w_idx]
        param[19,I] = uc_rho[I, y_idx]
        param[20,I] = vr_rho[I]
        # param[21,I] = r_rho[4*I-1]
        # param[22,I] = r_rho[4*I]
        # param[23,I] = uc_rho[7*I-3]
        # param[24,I] = uc_rho[7*I-2]
        # param[25,I] = uc_rho[7*I-1]
        # param[26,I] = uc_rho[7*I]
        param[27,I] = v[pg_idx] - z[pg_idx]
        param[28,I] = v[qg_idx] - z[qg_idx]
        t > 1 ? param[29,I] = r_v[I] - r_z[I] : param[29,I] = 0.
        param[30,I] = uc_v[I, v_idx] - uc_z[I, v_idx]
        param[31,I] = uc_v[I, w_idx] - uc_z[I, w_idx]
        param[32,I] = uc_v[I, y_idx] - uc_z[I, y_idx]
        t > 1 ? param[33,I] = uc_v[I, v_idx-3] - vr_z[I] : param[33,I] = 0.
        param[34,I] = ramp_limit[I]
        param[35,I] = ramp_limit[I]
        param[36,I] = ramp_limit[I]
        param[37,I] = ramp_limit[I]
        param[38,I] = pgmax[I]
        param[39,I] = pgmin[I]
        param[40,I] = qgmax[I]
        param[41,I] = qgmin[I]

        # Initialization of Augmented Lagrangian Method Parameters
        if major_iter <= 1
            # mu (Lagrangian term parameter)
            param[10,I] = 10.0
            param[11,I] = 10.0
            param[12,I] = 10.0
            param[13,I] = 10.0
            # xi (augmented term parameter)
            param[23,I] = 10.0
            param[24,I] = 10.0
            param[25,I] = 10.0
            param[26,I] = 10.0
            if t > 1
                # mu
                param[8,I] = 10.0
                param[9,I] = 10.0
                # xi
                param[21,I] = 10.0
                param[22,I] = 10.0
            end
        # else
        #     xi = param[43,I]
        end

        CUDA.sync_threads()


        t > 1 ? _st = 1 : _st = 3
        etas = CuDynamicSharedArray(Float64, 6, (3*n)*sizeof(Float64))
        omegas = CuDynamicSharedArray(Float64, 6, (3*n+6)*sizeof(Float64))
        cvios = CuDynamicSharedArray(Float64, 6, (3*n+12)*sizeof(Float64))
        for i in _st:6
            etas[i] = 1 / param[20+i,I]^0.1
            omegas[i] = 1 / param[20+i,I]
        end

        CUDA.sync_threads()

        it = 0
        terminate = false

        while !terminate
            it += 1

            status, minor_iter = ucmp_tron_generator_kernel(t, n, 500, 200, 1e-6, scale, x, xl, xu, param, c2, c1, c0, baseMVA)
            # Check the termination condition.
            if tx == 1
                cvios[1] = 0.0
                cvios[2] = 0.0
                cvios[3] = x[1] - pgmax[I]*x[4] - x[10]
                cvios[4] = x[1] - pgmin[I]*x[4] + x[11]
                cvios[5] = x[2] - qgmax[I]*x[4] - x[12]
                cvios[6] = x[2] - qgmin[I]*x[4] + x[13]
                if t > 1
                    cvios[1] = x[1] - x[3] - param[34,I]*x[7] - param[35,I]*x[5] - x[8]
                    cvios[2] = x[1] - x[3] + param[36,I]*x[4] + param[37,I]*x[6] + x[9]
                end
            end

            CUDA.sync_threads()

            terminate = true
            for i in _st:6
                cvios_i = cvios[i]
                etas_i = etas[i]
                CUDA.sync_threads()

                if abs(cvios_i) <= etas_i
                    if abs(cvios_i) > 1e-6
                        terminate = false
                        if tx == 1
                            param[7+i,I] += param[20+i,I]*cvios[i]
                            etas[i] = etas[i] / param[20+i,I]^0.9
                            omegas[i] = omegas[i] / param[20+i,I]
                        end
                    end
                else
                    terminate = false
                    if tx == 1
                        param[20+i,I] = min(xi_max, param[20+i,I]*10)
                        etas[i] = 1 / param[20+i,I]^0.1
                        omegas[i] = 1 / param[20+i,I]
                    end
                end

				CUDA.sync_threads()
            end
            if it >= max_auglag
                terminate = true
            end
            CUDA.sync_threads()
        end

        u[pg_idx] = x[1]
        u[qg_idx] = x[2]
        uc_u[I, v_idx] = x[4]
        uc_u[I, w_idx] = x[5]
        uc_u[I, y_idx] = x[6]
        uc_s[I, 4*t-3] = x[10]
        uc_s[I, 4*t-2] = x[11]
        uc_s[I, 4*t-1] = x[12]
        uc_s[I, 4*t] = x[13]

        if t > 1
            r_u[I] = x[3]
            vr_u[I] = x[7]
            r_s[2*I-1] = x[8]
            r_s[2*I] = x[9]
        end

        CUDA.sync_threads()
    end
    return
end
