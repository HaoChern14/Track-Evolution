using LinearAlgebra
using DifferentialEquations
using StaticArrays
using DelimitedFiles


# wavelet parameters
KK = 5;
n_intervals = 2^(KK-1);
MM = 3;

wavelet_index_mapping(n, l) = l * n_intervals + n
M_index_mapping(m, m1, m2) = MM^2 * m + MM * m1 + m2 + 1
M_index_mapping(m, m1, m2, m3) = MM^3 * m + MM^2 * m1 + MM * m2 + m3  + 1

# import initial conditions
initial_condition_q = readdlm("data/Tq_initial_coeff.dat")
initial_condition_g = readdlm("data/Tg_initial_coeff.dat")
initial_value = vcat(initial_condition_q, initial_condition_g);

mu0 = readdlm("data/initial_scale.txt")[1]
target_scales = readdlm("data/target_scales.txt")

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

gg_LO_array=Array{Any}(undef, MM^3);
qg_LO_array=Array{Any}(undef, MM^3);
gq_LO_array=Array{Any}(undef, MM^3);
qq_LO_array=Array{Any}(undef, MM^3);

LO_path = "kernels/LO_arrays/";
for m in 0:(MM-1)
    for m1 in 0:(MM-1)
        for m2 in 0:(MM-1)
            gg_file_name = string(LO_path, "gg_", m, "_", m1, "_", m2, ".dat")
            gg_LO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(gg_file_name), n_intervals, n_intervals, n_intervals)
            qg_file_name = string(LO_path, "qg_", m, "_", m1, "_", m2, ".dat")
            qg_LO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(qg_file_name), n_intervals, n_intervals, n_intervals)
            gq_file_name = string(LO_path, "gq_", m, "_", m1, "_", m2, ".dat")
            gq_LO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(gq_file_name), n_intervals, n_intervals, n_intervals)
            qq_file_name = string(LO_path, "qq_", m, "_", m1, "_", m2, ".dat")
            qq_LO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(qq_file_name), n_intervals, n_intervals, n_intervals)
        end
    end
end

gg_NLO_array=Array{Any}(undef, MM^3);
qg_NLO_array=Array{Any}(undef, MM^3);
gq_NLO_array=Array{Any}(undef, MM^3);
qq_NLO_array=Array{Any}(undef, MM^3);

NLO_path = "kernels/NLO_arrays/";
for m in 0:(MM-1)
    for m1 in 0:(MM-1)
        for m2 in 0:(MM-1)
            gg_nf_file_name = string(NLO_path, "gg_nf_NLO_TT_", m, "_", m1, "_", m2, ".dat")
            gg_no_nf_file_name = string(NLO_path, "gg_no_nf_NLO_TT_", m, "_", m1, "_", m2, ".dat") 
            gg_NLO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(gg_no_nf_file_name) + nf * readdlm(gg_nf_file_name), n_intervals, n_intervals, n_intervals)
            
            qg_no_nf_file_name = string(NLO_path, "qg_NLO_TT_", m, "_", m1, "_", m2, ".dat") 
            qg_NLO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(qg_no_nf_file_name), n_intervals, n_intervals, n_intervals)
            
            gq_no_nf_file_name = string(NLO_path, "gq_NLO_TT_", m, "_", m1, "_", m2, ".dat") 
            gq_NLO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(gq_no_nf_file_name), n_intervals, n_intervals, n_intervals)

            qq_nf_file_name = string(NLO_path, "qq_nf_NLO_TT_", m, "_", m1, "_", m2, ".dat")
            qq_no_nf_file_name = string(NLO_path, "qq_no_nf_NLO_TT_", m, "_", m1, "_", m2, ".dat") 
            qq_NLO_array[M_index_mapping(m, m1, m2)] = reshape(readdlm(qq_no_nf_file_name) + nf * readdlm(qq_nf_file_name), n_intervals, n_intervals, n_intervals)
        end
    end
end

ddu_array=Array{Any}(undef, MM^4);
dud_array=Array{Any}(undef, MM^4);
udd_array=Array{Any}(undef, MM^4);

ggg_array=Array{Any}(undef, MM^4);

qgg_array=Array{Any}(undef, MM^4);
gqg_array=Array{Any}(undef, MM^4);
ggq_array=Array{Any}(undef, MM^4);

gqq_array=Array{Any}(undef, MM^4);
qgq_array=Array{Any}(undef, MM^4);
qqg_array=Array{Any}(undef, MM^4);

qbqq_array=Array{Any}(undef, MM^4);
qqbq_array=Array{Any}(undef, MM^4);
qqqb_array=Array{Any}(undef, MM^4);

for m in 0:(MM-1)
    for m1 in 0:(MM-1)
        for m2 in 0:(MM-1)
            for m3 in 0:(MM-1)
                ddu_file_name = string(NLO_path, "ddu_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                ddu_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(ddu_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                dud_file_name = string(NLO_path, "dud_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                dud_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(dud_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                udd_file_name = string(NLO_path, "udd_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                udd_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(udd_file_name), n_intervals, n_intervals, n_intervals, n_intervals)

                ggg_file_name = string(NLO_path, "ggg_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                ggg_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(ggg_file_name), n_intervals, n_intervals, n_intervals, n_intervals)

                qgg_file_name = string(NLO_path, "qgg_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                qgg_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(qgg_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                gqg_file_name = string(NLO_path, "gqg_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                gqg_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(gqg_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                ggq_file_name = string(NLO_path, "ggq_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                ggq_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(ggq_file_name), n_intervals, n_intervals, n_intervals, n_intervals)

                gqq_file_name = string(NLO_path, "gqq_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                gqq_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(gqq_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                qgq_file_name = string(NLO_path, "qgq_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                qgq_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(qgq_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                qqg_file_name = string(NLO_path, "qqg_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                qqg_array[M_index_mapping(m, m1, m2, m3)] = reshape(readdlm(qqg_file_name), n_intervals, n_intervals, n_intervals, n_intervals)


                # in kernel generation, forget to multiply a factor of 2
                qbqq_file_name = string(NLO_path, "qbqq_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                qbqq_array[M_index_mapping(m, m1, m2, m3)] = reshape(2 * readdlm(qbqq_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                qqbq_file_name = string(NLO_path, "qqbq_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                qqbq_array[M_index_mapping(m, m1, m2, m3)] = reshape(2 * readdlm(qqbq_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
                qqqb_file_name = string(NLO_path, "qqqb_NLO_TTT_", m, "_", m1, "_", m2, "_", m3, ".dat")
                qqqb_array[M_index_mapping(m, m1, m2, m3)] = reshape(2 * readdlm(qqqb_file_name), n_intervals, n_intervals, n_intervals, n_intervals)
            end
        end
    end
end


NLO_Tq_coeff=CF*nf*TF*(-532/27 + (92*log(2))/9 + (16*log(2)^2)/3) + CA*CF*(3337/54 - (274*log(2))/9 - (44*log(2)^2)/3 - 49*zeta3) + CF^2*(3/2 - 2*pi^2 + 24*zeta3);
NLO_Tg_coeff= nf*TF*(-4*CF + CA*(-658/27 + (8*pi^2)/9 + (92*log(2))/9 + (16*log(2)^2)/3)) + CA^2*(1880/27 - (22*pi^2)/9 - (274*log(2))/9 - (44*log(2)^2)/3 - 25*zeta3);


gg_LO_mat(n, m, m1, m2) = gg_LO_array[M_index_mapping(m, m1, m2)][:, :, n]
qg_LO_mat(n, m, m1, m2) = qg_LO_array[M_index_mapping(m, m1, m2)][:, :, n]
gq_LO_mat(n, m, m1, m2) = gq_LO_array[M_index_mapping(m, m1, m2)][:, :, n]
qq_LO_mat(n, m, m1, m2) = qq_LO_array[M_index_mapping(m, m1, m2)][:, :, n]

gg_NLO_mat(n, m, m1, m2) = gg_NLO_array[M_index_mapping(m, m1, m2)][:, :, n]
qg_NLO_mat(n, m, m1, m2) = qg_NLO_array[M_index_mapping(m, m1, m2)][:, :, n]
gq_NLO_mat(n, m, m1, m2) = gq_NLO_array[M_index_mapping(m, m1, m2)][:, :, n]
qq_NLO_mat(n, m, m1, m2) = qq_NLO_array[M_index_mapping(m, m1, m2)][:, :, n]

ddu_mat(n, m, m1, m2, m3) = ddu_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
dud_mat(n, m, m1, m2, m3) = dud_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
udd_mat(n, m, m1, m2, m3) = udd_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]

ggg_mat(n, m, m1, m2, m3) = ggg_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]

qgg_mat(n, m, m1, m2, m3) = qgg_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
gqg_mat(n, m, m1, m2, m3) = gqg_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
ggq_mat(n, m, m1, m2, m3) = ggq_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]

gqq_mat(n, m, m1, m2, m3) = gqq_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
qgq_mat(n, m, m1, m2, m3) = qgq_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
qqg_mat(n, m, m1, m2, m3) = qqg_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]

qbqq_mat(n, m, m1, m2, m3) = qbqq_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
qqbq_mat(n, m, m1, m2, m3) = qqbq_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]
qqqb_mat(n, m, m1, m2, m3) = qqqb_array[M_index_mapping(m, m1, m2, m3)][:, :, :, n]


function array_multiply(array::Array{Float64, 3}, vec_n1, vec_n2, vec_n3)
    sum(vec_n1[i] * dot(vec_n3, array[:, :, i], vec_n2) for i = 1:n_intervals )
end

array_multiply(array::Array{Float64, 3}, single_vec) = array_multiply(array, single_vec, single_vec, single_vec)


q_coeff_vec = Array{Any}(undef, MM);
g_coeff_vec = Array{Any}(undef, MM);



# evolution equation (to higher scale)
function evolution_up!(du, u, p, t)
    for m = 0:(MM-1)
        q_coeff_vec[m + 1] =  SVector{n_intervals}([u[i] for i = (1 + m * n_intervals):((m + 1) * n_intervals)])
        g_coeff_vec[m + 1] =  SVector{n_intervals}([u[i] for i = (1 + (m + 3) * n_intervals):((m + 4) * n_intervals)])
    end
    for n = 1:n_intervals
        for m = 0:(MM-1)
            # T_q evolution
            du[wavelet_index_mapping(n, m)] =  a_s(t) * (3 * CF * u[wavelet_index_mapping(n, m)] + sum(dot(g_coeff_vec[m2 + 1], qg_LO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1]) + dot(q_coeff_vec[m2 + 1], gq_LO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1)) ) + a_s(t)^2 * (
                NLO_Tq_coeff * u[wavelet_index_mapping(n, m)] + sum(dot(g_coeff_vec[m2 + 1], qg_NLO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1]) + dot(q_coeff_vec[m2 + 1], gq_NLO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1))  
                # one quark + two gluons part
                + sum( array_multiply(qgg_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                + sum( array_multiply(gqg_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                + sum( array_multiply(ggq_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # three quarks part
                # + sum( array_multiply(qbqq_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + sum( array_multiply(qqbq_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + sum( array_multiply(qqqb_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + 2 * (nf - 1) * sum( array_multiply(udd_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + 2 * (nf - 1) * sum( array_multiply(dud_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + 2 * (nf - 1) * sum( array_multiply(ddu_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))

                # simplification for quarks being the same
                + sum( array_multiply(
                    qbqq_mat(n, m, m1, m2, m3) + qqbq_mat(n, m, m1, m2, m3) + qqqb_mat(n, m, m1, m2, m3) + 2 * (nf - 1) * (udd_mat(n, m, m1, m2, m3) + dud_mat(n, m, m1, m2, m3) + ddu_mat(n, m, m1, m2, m3))
                    , q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
            )
            # T_g evolution
            du[wavelet_index_mapping(n, m) + MM * n_intervals] =  a_s(t) * (beta0 * u[wavelet_index_mapping(n, m) + MM * n_intervals] + sum(dot(g_coeff_vec[m2 + 1], gg_LO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) + 2 * nf * dot(q_coeff_vec[m2 + 1], qq_LO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1])  for m1=0:(MM - 1), m2=0:(MM - 1))) + a_s(t)^2 * (
                NLO_Tg_coeff * u[wavelet_index_mapping(n, m) + MM * n_intervals] +  sum(dot(g_coeff_vec[m2 + 1], gg_NLO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) + 2 * nf * dot(q_coeff_vec[m2 + 1], qq_NLO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1])  for m1=0:(MM - 1), m2=0:(MM - 1))  
                # three gluons part
                + sum( array_multiply(ggg_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # one gluon + two quarks part
                + 2 * nf *(
                    sum( array_multiply(gqq_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                    + sum( array_multiply(qgq_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1)) 
                    + sum( array_multiply(qqg_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                )

            )
        end
    end
end

# evolution equation (to lower scale)
function evolution_down!(du, u, p, t)
    for m = 0:(MM-1)
        q_coeff_vec[m + 1] =  SVector{n_intervals}([u[i] for i = (1 + m * n_intervals):((m + 1) * n_intervals)])
        g_coeff_vec[m + 1] =  SVector{n_intervals}([u[i] for i = (1 + (m + 3) * n_intervals):((m + 4) * n_intervals)])
    end
    for n = 1:n_intervals
        for m = 0:(MM-1)
            # T_q evolution
            du[wavelet_index_mapping(n, m)] = -(  a_s(-t) * (3 * CF * u[wavelet_index_mapping(n, m)] + sum(dot(g_coeff_vec[m2 + 1], qg_LO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1]) + dot(q_coeff_vec[m2 + 1], gq_LO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1)) ) + a_s(-t)^2 * (
                NLO_Tq_coeff * u[wavelet_index_mapping(n, m)] + sum(dot(g_coeff_vec[m2 + 1], qg_NLO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1]) + dot(q_coeff_vec[m2 + 1], gq_NLO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1))  
                # one quark + two gluons part
                + sum( array_multiply(qgg_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                + sum( array_multiply(gqg_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                + sum( array_multiply(ggq_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # three quarks part
                # + sum( array_multiply(qbqq_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + sum( array_multiply(qqbq_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + sum( array_multiply(qqqb_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + 2 * (nf - 1) * sum( array_multiply(udd_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + 2 * (nf - 1) * sum( array_multiply(dud_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # + 2 * (nf - 1) * sum( array_multiply(ddu_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))

                # simplification for quarks being the same
                + sum( array_multiply(
                    qbqq_mat(n, m, m1, m2, m3) + qqbq_mat(n, m, m1, m2, m3) + qqqb_mat(n, m, m1, m2, m3) + 2 * (nf - 1) * (udd_mat(n, m, m1, m2, m3) + dud_mat(n, m, m1, m2, m3) + ddu_mat(n, m, m1, m2, m3))
                    , q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
            ) )
            # T_g evolution
            du[wavelet_index_mapping(n, m) + MM * n_intervals] = -(  a_s(-t) * (beta0 * u[wavelet_index_mapping(n, m) + MM * n_intervals] + sum(dot(g_coeff_vec[m2 + 1], gg_LO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) + 2 * nf * dot(q_coeff_vec[m2 + 1], qq_LO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1])  for m1=0:(MM - 1), m2=0:(MM - 1))) + a_s(-t)^2 * (
                NLO_Tg_coeff * u[wavelet_index_mapping(n, m) + MM * n_intervals] +  sum(dot(g_coeff_vec[m2 + 1], gg_NLO_mat(n, m, m1, m2), g_coeff_vec[m1 + 1]) + 2 * nf * dot(q_coeff_vec[m2 + 1], qq_NLO_mat(n, m, m1, m2), q_coeff_vec[m1 + 1])  for m1=0:(MM - 1), m2=0:(MM - 1))  
                # three gluons part
                + sum( array_multiply(ggg_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                # one gluon + two quarks part
                + 2 * nf *(
                    sum( array_multiply(gqq_mat(n, m, m1, m2, m3), g_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                    + sum( array_multiply(qgq_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], g_coeff_vec[m2 + 1], q_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1)) 
                    + sum( array_multiply(qqg_mat(n, m, m1, m2, m3), q_coeff_vec[m1 + 1], q_coeff_vec[m2 + 1], g_coeff_vec[m3 + 1]) for m1=0:(MM - 1), m2=0:(MM - 1), m3=0:(MM - 1))
                )

            ) )
        end
    end
end

# export wavelet coefficients
function NLO_coefficients_export(mu, vec)
    outfile = string("data/output/NLO_wavelet_coeff_mu_", mu, ".txt")
    open(outfile, "w") do os
        writedlm(os, vec)
    end
end

# evolve to higher energy scale
if !isempty(up_scales)
    upper_bound = ceil(mu_to_t(maximum(up_scales)))
    

    # solving evolution equation
    tspan = (0.0, upper_bound)
    prob = ODEProblem(evolution_up!, initial_value, tspan)
    sol = solve(prob, reltol=1e-10, abstol=1e-10)

    # export desired wavelet coefficients
    broadcast(mu->NLO_coefficients_export(mu, sol(mu_to_t(mu))), up_scales)
end


# evolve to lower energy scale
if !isempty(down_scales)
    lower_bound = ceil(- mu_to_t(minimum(down_scales)))
    
    # solving evolution equation
    tspan = (0.0, lower_bound)
    prob = ODEProblem(evolution_down!, initial_value, tspan)
    sol = solve(prob, reltol=1e-10, abstol=1e-10)

    # export desired wavelet coefficients
    broadcast(mu->NLO_coefficients_export(mu, sol(-mu_to_t(mu))), down_scales)
end