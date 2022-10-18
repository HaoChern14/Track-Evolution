using LinearAlgebra
using DifferentialEquations
using StaticArrays
using DelimitedFiles


mu0 = readdlm("track_data/initial_scale.txt")[1]
target_scales = readdlm("track_data/target_scales.txt")


up_scales = filter(x -> x > mu0, target_scales)
down_scales = filter(x -> x < mu0, target_scales)


# constants

const nf = 5
const CF = 4/3
const CA = 3
const TF = 1/2
const zeta3 = 1.20205690315959

const beta0 = (11*CA)/3 - (4*nf*TF)/3
const beta1 = (34*CA^2)/3 - (20*CA*nf*TF)/3 - 4*CF*nf*TF
const beta2 = (2857*CA^3)/54 + ((-1415*CA^2)/27 - (205*CA*CF)/9 + 2*CF^2)*nf*TF + ((158*CA)/27 + (44*CF)/9)*nf^2*TF^2
const beta3 = 149753/6 + (1093*nf^3)/729 + nf*(-1078361/162 - (6508*zeta3)/27) + 3564*zeta3 + nf^2*(50065/162 + (6472*zeta3)/81)

const Beta1 = beta1/beta0
const Beta2 = beta2/beta0
const Beta3 = beta3/beta0

# running coupling
const lambda_QCD = 213/1000
ll0 = 2*log(mu0/lambda_QCD)

mu_to_t(mu) = 2*log(mu) - 2*log(mu0)
t_to_mu(t) = mu0 * exp(t/2) 

function a_s(t) 
    # t = log(mu^2) - log(mu0^2) = 2*log(mu) - 2*log(mu0) = 2*log(mu/lambda_QCD) - 2*log(mu0/lambda_QCD)
    ll = t + ll0 # ll = 2*log(mu/lambda_QCD) 
    LL = beta0 * ll
    return (1 + Beta3/(2*LL^3) + Beta2/LL^2 - (3*Beta1*Beta2*log(ll))/LL^3 - (Beta1*log(ll))/LL + (Beta1^2*(-1 - log(ll) + log(ll)^2))/LL^2 + (Beta1^3*(-1/2 + 2*log(ll) + (5*log(ll)^2)/2 - log(ll)^3))/LL^3)/LL
end

# import kernels

LO_order = 100;
NLO_order = 40;
LO_path_name = string("kernels/LO_kernels/order_", LO_order, "_mat/");
NLO_TT_path = string("kernels/NLO_TT_kernels/order_", NLO_order, "_mat/");
NLO_TTT_path = string("kernels/NLO_TTT_kernels/");

qg_LO_mat = [readdlm(string(LO_path_name, "qg_array_", n, ".txt")) for n = 1:LO_order];
gg_LO_mat = [readdlm(string(LO_path_name, "gg_array_", n, ".txt")) for n = 1:LO_order];
qqb_LO_mat = [readdlm(string(LO_path_name, "qqb_array_", n, ".txt")) for n = 1:LO_order];

qg_NLO_mat = [readdlm(string(NLO_TT_path, "qg_no_nf_array_", n, ".txt")) + nf*readdlm(string(NLO_TT_path, "qg_nf_array_", n, ".txt")) for n = 1:NLO_order];
gg_NLO_mat = [readdlm(string(NLO_TT_path, "gg_no_nf_array_", n, ".txt")) + nf*readdlm(string(NLO_TT_path, "gg_nf_array_", n, ".txt")) for n = 1:NLO_order];
qqb_NLO_mat = [readdlm(string(NLO_TT_path, "qq_no_nf_array_", n, ".txt")) + nf*readdlm(string(NLO_TT_path, "qq_nf_array_", n, ".txt")) for n = 1:NLO_order];

NLO_Tq_coeff=(CF*nf*TF*(-206/27 - (8*log(2))/3) + CA*CF*(1715/54 - 45*zeta3) + 
  CF^2*(3/2 - 4*pi^2 + 24*zeta3))-2*(-CF^2*pi^2);

NLO_Tg_coeff= (nf*TF*(-4*CF + CA*(-332/27 + (16*pi^2)/9 - (8*log(2))/3)) + 
  CA^2*(1069/27 - (44*pi^2)/9 - 21*zeta3)) - 2*(-11*CA + 4*nf*TF)*CA*pi^2/9;

function NLO_TTT_import_array(tag_name, order)
    path_name = string(NLO_TTT_path, tag_name, "_mat/order_", order, "_mat/") 
    array_size = 2 * order + 1
    return [reshape(readdlm(string(path_name, tag_name, "_array_", n, ".txt")), array_size, array_size, array_size) for n = 1:order]
end

ggg_NLO_mat = NLO_TTT_import_array("ggg", NLO_order);
qgg_NLO_mat = NLO_TTT_import_array("qgg", NLO_order);
gqq_NLO_mat = NLO_TTT_import_array("gqq", NLO_order);
qqq_NLO_mat = NLO_TTT_import_array("qqq", NLO_order);
qQQb_NLO_mat = NLO_TTT_import_array("qQQb", NLO_order);

function array_multiply(arr, vec1, vec2, vec3, vec_size)
    sum(vec3[i]*dot(conj(vec1), arr[:, :, i], vec2) for i = 1:vec_size)
end

function matrix_multiply(matrix, vec1, vec2)
    dot(conj(vec1), matrix, vec2)
end

# import initial conditions
Tq_initial_value_re = readdlm("track_data/Tq_initial_coeff_re.txt")
Tg_initial_value_re = readdlm("track_data/Tg_initial_coeff_re.txt")
Tq_initial_value_im = readdlm("track_data/Tq_initial_coeff_im.txt")
Tg_initial_value_im = readdlm("track_data/Tg_initial_coeff_im.txt")

Tq_initial_value = Tq_initial_value_re + im * Tq_initial_value_im
Tg_initial_value = Tg_initial_value_re + im * Tg_initial_value_im
initial_value = vcat(Tq_initial_value, Tg_initial_value)


# evolution equations to higher scale
function evolution_up!(du, u, p, t)
    q_vector = vcat(reverse([conj(u[i]) for i = 1:LO_order]), [1], [u[i] for i = 1:LO_order])
    g_vector = vcat(reverse([conj(u[i]) for i = (LO_order + 1):(2 * LO_order)]), [1], [u[i] for i = (LO_order + 1):(2 * LO_order)])
    q_vector_trunc = SVector{2 * NLO_order + 1}(vcat(reverse([conj(u[i]) for i = 1:NLO_order]), [1], [u[i] for i = 1:NLO_order]))
    g_vector_trunc = SVector{2 * NLO_order + 1}(vcat(reverse([conj(u[i]) for i = (LO_order + 1):(LO_order + NLO_order)]), [1], [u[i] for i = (LO_order + 1):(LO_order + NLO_order)]))
    for i = 1:NLO_order 
        # T_q evolution
        du[i] = a_s(t) * (3 * CF * u[i] + dot(conj(q_vector), qg_LO_mat[i], g_vector)) + a_s(t)^2 * (
            NLO_Tq_coeff * u[i] + dot(conj(q_vector_trunc), qg_NLO_mat[i], g_vector_trunc)
            + array_multiply(qgg_NLO_mat[i], q_vector_trunc, g_vector_trunc, g_vector_trunc, 2*NLO_order+1)
            + array_multiply(qqq_NLO_mat[i], q_vector_trunc, q_vector_trunc, q_vector_trunc, 2*NLO_order+1)
            + (nf - 1) * array_multiply(qQQb_NLO_mat[i], q_vector_trunc, q_vector_trunc, q_vector_trunc, 2*NLO_order+1)
        )
        # T_g evolution
        du[i + LO_order] = a_s(t) * (beta0 * u[i + LO_order] + dot(conj(g_vector), gg_LO_mat[i], g_vector) + nf * dot(conj(q_vector), qqb_LO_mat[i], q_vector)) + a_s(t)^2 * (
            NLO_Tg_coeff * u[i + LO_order] + dot(conj(g_vector_trunc), gg_NLO_mat[i], g_vector_trunc) + nf * dot(conj(q_vector_trunc), qqb_NLO_mat[i], q_vector_trunc)
            + array_multiply(ggg_NLO_mat[i], g_vector_trunc, g_vector_trunc, g_vector_trunc, 2*NLO_order+1)
            + nf * array_multiply(gqq_NLO_mat[i], g_vector_trunc, q_vector_trunc, q_vector_trunc, 2*NLO_order+1)
        )
    end
    
    for i = (NLO_order+1):LO_order
        # T_q evolution
        du[i] = a_s(t) * (3 * CF * u[i] + dot(conj(q_vector), qg_LO_mat[i], g_vector))
        # T_g evolution
        du[i + LO_order] = a_s(t) * (beta0 * u[i + LO_order] + dot(conj(g_vector), gg_LO_mat[i], g_vector) + nf * dot(conj(q_vector), qqb_LO_mat[i], q_vector))
    end
end

# evolution equations to lower scale
function evolution_down!(du, u, p, t)
    q_vector = vcat(reverse([conj(u[i]) for i = 1:LO_order]), [1], [u[i] for i = 1:LO_order])
    g_vector = vcat(reverse([conj(u[i]) for i = (LO_order + 1):(2 * LO_order)]), [1], [u[i] for i = (LO_order + 1):(2 * LO_order)])
    q_vector_trunc = SVector{2 * NLO_order + 1}(vcat(reverse([conj(u[i]) for i = 1:NLO_order]), [1], [u[i] for i = 1:NLO_order]))
    g_vector_trunc = SVector{2 * NLO_order + 1}(vcat(reverse([conj(u[i]) for i = (LO_order + 1):(LO_order + NLO_order)]), [1], [u[i] for i = (LO_order + 1):(LO_order + NLO_order)]))
    for i = 1:NLO_order 
        # T_q evolution
        du[i] = -( a_s(-t) * (3 * CF * u[i] + dot(conj(q_vector), qg_LO_mat[i], g_vector)) + a_s(-t)^2 * (
            NLO_Tq_coeff * u[i] + dot(conj(q_vector_trunc), qg_NLO_mat[i], g_vector_trunc)
            + array_multiply(qgg_NLO_mat[i], q_vector_trunc, g_vector_trunc, g_vector_trunc, 2*NLO_order+1)
            + array_multiply(qqq_NLO_mat[i], q_vector_trunc, q_vector_trunc, q_vector_trunc, 2*NLO_order+1)
            + (nf - 1) * array_multiply(qQQb_NLO_mat[i], q_vector_trunc, q_vector_trunc, q_vector_trunc, 2*NLO_order+1)
        ) )
        # T_g evolution
        du[i + LO_order] = -( a_s(-t) * (beta0 * u[i + LO_order] + dot(conj(g_vector), gg_LO_mat[i], g_vector) + nf * dot(conj(q_vector), qqb_LO_mat[i], q_vector)) + a_s(-t)^2 * (
            NLO_Tg_coeff * u[i + LO_order] + dot(conj(g_vector_trunc), gg_NLO_mat[i], g_vector_trunc) + nf * dot(conj(q_vector_trunc), qqb_NLO_mat[i], q_vector_trunc)
            + array_multiply(ggg_NLO_mat[i], g_vector_trunc, g_vector_trunc, g_vector_trunc, 2*NLO_order+1)
            + nf * array_multiply(gqq_NLO_mat[i], g_vector_trunc, q_vector_trunc, q_vector_trunc, 2*NLO_order+1)
        ) )
    end
    
    for i = (NLO_order+1):LO_order
        # T_q evolution
        du[i] = - a_s(-t) * (3 * CF * u[i] + dot(conj(q_vector), qg_LO_mat[i], g_vector))
        # T_g evolution
        du[i + LO_order] = - a_s(-t) * (beta0 * u[i + LO_order] + dot(conj(g_vector), gg_LO_mat[i], g_vector) + nf * dot(conj(q_vector), qqb_LO_mat[i], q_vector))
    end
end


# track_data output
# export fourier coefficients
function NLO_coefficients_export(mu, vec)
    outfile = string("track_data/output/NLO_fourier_coeff_mu_", mu, ".txt")
    open(outfile, "w") do os
        writedlm(os, vec)
    end
end
# export discretized curves
# export x-axis
x = 0.0:0.001:1.0
open(string("track_data/output/x_axis.txt"), "w") do os
    writedlm(os, x)
end
# export track values
function track_value_export(mu, Tq_value_vec, Tg_value_vec)
    outfile = string("track_data/output/Tq_value_mu_", mu, ".txt")
    open(outfile, "w") do os
        writedlm(os, Tq_value_vec)
    end
    outfile = string("track_data/output/Tg_value_mu_", mu, ".txt")
    open(outfile, "w") do os
        writedlm(os, Tg_value_vec)
    end
end





# evolve to higher energy scale
if !isempty(up_scales)
    upper_bound = ceil(mu_to_t(maximum(up_scales)))

    # solving evolution equation
    tspan = (0.0, upper_bound)
    prob = ODEProblem(evolution_up!, initial_value, tspan)
    sol = solve(prob, reltol=1e-10, abstol=1e-10)


    # export desired fourier coefficients
    broadcast(mu->NLO_coefficients_export(mu, sol(mu_to_t(mu))), up_scales)


    # export curves
    Tq_res_t(t, x) = 1 + 2 * real(dot(sol(t)[1:LO_order], [exp(- im * 2 * n * pi * x) for n = 1:LO_order]))
    Tg_res_t(t, x) = 1 + 2 * real(dot(sol(t)[(LO_order + 1):(2 * LO_order)], [exp(- im * 2 * n * pi * x) for n = 1:LO_order]))
    Tq_res_mu(mu, x) = Tq_res_t(mu_to_t(mu), x)
    Tg_res_mu(mu, x) = Tg_res_t(mu_to_t(mu), x)

    broadcast(mu -> track_value_export(mu, [Tq_res_mu(mu, 0.001 * i) for i = 0:1000], [Tg_res_mu(mu, 0.001 * i) for i = 0:1000]), up_scales)

end


# evolve to lower energy scale
if !isempty(down_scales)
    lower_bound = ceil(- mu_to_t(minimum(down_scales)))
    
    # solving evolution equation
    tspan = (0.0, lower_bound)
    prob = ODEProblem(evolution_down!, initial_value, tspan)
    sol = solve(prob, reltol=1e-10, abstol=1e-10)

    # export desired fourier coefficients
    broadcast(mu->NLO_coefficients_export(mu, sol(-mu_to_t(mu))), down_scales)

    # export curves
    Tq_res_t(t, x) = 1 + 2 * real(dot(sol(t)[1:LO_order], [exp(- im * 2 * n * pi * x) for n = 1:LO_order]))
    Tg_res_t(t, x) = 1 + 2 * real(dot(sol(t)[(LO_order + 1):(2 * LO_order)], [exp(- im * 2 * n * pi * x) for n = 1:LO_order]))
    Tq_res_mu(mu, x) = Tq_res_t(-mu_to_t(mu), x)
    Tg_res_mu(mu, x) = Tg_res_t(-mu_to_t(mu), x)

    broadcast(mu -> track_value_export(mu, [Tq_res_mu(mu, 0.001 * i) for i = 0:1000], [Tg_res_mu(mu, 0.001 * i) for i = 0:1000]), down_scales)

end