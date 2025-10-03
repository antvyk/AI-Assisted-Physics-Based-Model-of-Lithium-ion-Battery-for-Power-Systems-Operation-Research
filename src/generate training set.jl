
using DataFrames
using XLSX
using Dierckx
using DelimitedFiles

println("Start")

directory = 
directory_for_parameters = 
file_with_SPM_parameters = directory_for_parameters*"Parameters for SPM - LG Chen2020Kane2022.xlsx"
SPM_Parameters = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Parameters"))
OCP_coefficients_p_df = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Positive Electrode"))
OCP_coefficients_n_df = DataFrame(XLSX.readtable(file_with_SPM_parameters, "Negative Electrode"))
OCP_coefficients_p = Matrix(OCP_coefficients_p_df)
OCP_coefficients_n = Matrix(OCP_coefficients_n_df)
spl_OCP_p = Spline1D(OCP_coefficients_p[:,1], OCP_coefficients_p[:,2], k=1)
spl_OCP_n = Spline1D(OCP_coefficients_n[:,1], OCP_coefficients_n[:,2], k=1)


function OCP_p_spline(theta)
    evaluate(spl_OCP_p, theta)
end

function OCP_n_spline(theta)
    evaluate(spl_OCP_n, theta)
end

######################################################################################

const R_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Particle radius","Positive electrode"][1]
const R_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Particle radius","Negative electrode"][1]
const El_thickness_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode thickness","Positive electrode"][1]
const El_thickness_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode thickness","Negative electrode"][1]
const El_length_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode length","Positive electrode"][1]
const El_length_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode length","Negative electrode"][1]
const El_width_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode width","Positive electrode"][1]
const El_width_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Electrode width","Negative electrode"][1]
const eps_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Active material volume fraction","Positive electrode"][1]
const eps_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Active material volume fraction","Negative electrode"][1]
const k_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Reaction rate","Positive electrode"][1]
const k_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Reaction rate","Negative electrode"][1]
const alpha_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Charge transfer coefficient","Positive electrode"][1]
const alpha_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Charge transfer coefficient","Negative electrode"][1]
const D_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Diffusion coefficient","Positive electrode"][1]
const D_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Diffusion coefficient","Negative electrode"][1]
const c_max_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Maximum concentration","Positive electrode"][1]
const c_max_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Maximum concentration","Negative electrode"][1]
const c_el = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Li concentration in  Electrolyte","Electrolyte"][1]
const theta_soc0_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 0% SoC","Positive electrode"][1]
const theta_soc0_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 0% SoC","Negative electrode"][1]
const theta_soc1_p = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 100% SoC","Positive electrode"][1]
const theta_soc1_n = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Stoichiometry at 100% SoC","Negative electrode"][1]
const V_min = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Minimum Voltage","Nameplate data"][1]
const V_max = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Maximum Voltage","Nameplate data"][1]
const Q = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Cell Charge Capacity","Nameplate data"][1]
const Q_E = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Nominal energy capacity","Nameplate data"][1]


###   SEI parameters
const R_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI resistivity","SEI"][1]
const δ00_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Initial SEI thickness","SEI"][1]
const c_bulk_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Concentration of EC in bulk electrolyte","SEI"][1]
const k_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="Kinetic rate constant for SEI reaction","SEI"][1]
const α_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI charge transfer coefficient","SEI"][1]
const ρ_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI density","SEI"][1]
const M_sei = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI molar weight","SEI"][1]
const D_EC = SPM_Parameters[SPM_Parameters[:,"Parameter"].=="SEI layer diffusivity","SEI"][1]
###

#Thermodynamics Constants
const F = 96485.309
const R = 8.31451
const T = 293




const T_h = 1
const dt = 0.1
const N_dr = 40 #40
const N_time = trunc(Int,(3600/dt))*T_h

###################################################################################
const dr_p = R_p/N_dr
const dr_n = R_n/N_dr
PDE_p = zeros(N_dr+1, N_dr+1)
PDE_p[1,1] = dt*D_p/dr_p^2*(dr_p^2/dt/D_p-2)
PDE_p[1,2] = dt*D_p/dr_p^2*2
for i in 2:N_dr
    PDE_p[i,i-1] = dt*D_p/dr_p^2
    PDE_p[i,i] = dt*D_p/dr_p^2*(dr_p^2/dt/D_p-2/i-2)
    PDE_p[i,i+1] = dt*D_p/dr_p^2*(2/i+1)
end
PDE_p[N_dr+1,N_dr] = 2*dt*D_p/dr_p^2
PDE_p[N_dr+1,N_dr+1] = 1-2*dt*D_p/dr_p^2


PDE_n = zeros(N_dr+1, N_dr+1)
PDE_n[1,1] = dt*D_n/dr_n^2*(dr_n^2/dt/D_n-2)
PDE_n[1,2] = dt*D_n/dr_n^2*2
for i in 2:N_dr
    PDE_n[i,i-1] = dt*D_n/dr_n^2
    PDE_n[i,i] = dt*D_n/dr_n^2*(dr_n^2/dt/D_n-2/i-2)
    PDE_n[i,i+1] = dt*D_n/dr_n^2*(2/i+1)
end

PDE_n[N_dr+1,N_dr] = 2*dt*D_n/dr_n^2
PDE_n[N_dr+1,N_dr+1] = 1-2*dt*D_n/dr_n^2
###################################################################################


######################################################################################

function Final_SoC_Caploss(SoC_init,PowerInput,SoC_max_init)
    Capacity_loss_error = 0
    SoC_final = 0
    Capacity_loss  = 0
    δ_sei_prev = (1-SoC_max_init)*Q*3600*M_sei/ρ_sei/(2*F)+δ00_sei
    LLI_prev = 0
    SoC_prev = SoC_init
    SoC = 0
    Final_SoC_Caploss = 0
    LLI = LLI_prev
    C_EC_prev = c_bulk_EC
    c_initial_p = (theta_soc0_p*(1-SoC_init)+theta_soc1_p*SoC_init)* c_max_p
    c_initial_n = (theta_soc0_n*(1-SoC_init)+theta_soc1_n*SoC_init)*c_max_n
    P = PowerInput
    SPM_c_p_prev = c_initial_p * ones(N_dr+1,1)
    SPM_c_n_prev = c_initial_n * ones(N_dr+1,1)
    theta_p = c_initial_p/c_max_p
    theta_n = c_initial_n/c_max_n
    I = P/(OCP_p_spline(theta_p) - OCP_n_spline(theta_n))
  
    
    for k in 2:N_time
        J_p = -I/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p)
        J_n = I/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)
        SPM_c_p = PDE_p*SPM_c_p_prev
        SPM_c_p_N_dr_1 = SPM_c_p[N_dr+1]
        SPM_c_p[N_dr+1] = SPM_c_p_N_dr_1 + (-J_p/D_p)*(2*dt*D_p/dr_p+2*dt*D_p/R_p)
        theta_p = SPM_c_p[N_dr+1]/c_max_p
        SPM_c_n = PDE_n*SPM_c_n_prev
        SPM_c_n_N_dr_1 = SPM_c_n[N_dr+1]
        SPM_c_n[N_dr+1] = SPM_c_n_N_dr_1 + (-J_n/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
        theta_n =  SPM_c_n[N_dr+1]/c_max_n
      

        if theta_p>=theta_soc1_p && theta_p<=theta_soc0_p && theta_n<=theta_soc1_n && theta_n>=theta_soc0_n
            phi_p = 2*R*T/F*asinh(J_p*F/(k_p*theta_p^alpha_p*c_el^(1-alpha_p)*c_max_p*(1-theta_p)^(1-alpha_p))/2) + OCP_p_spline(theta_p)
            phi_n = 2*R*T/F*asinh(J_n*F/(k_n*theta_n^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-theta_n)^(1-alpha_n))/2) + OCP_n_spline(theta_n)
            I = P/(phi_p-phi_n)

            for kk in 1:3

                J_p = -I/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p)
                SPM_c_p[N_dr+1] = SPM_c_p_N_dr_1 + (-J_p/D_p)*(2*dt*D_p/dr_p+2*dt*D_p/R_p)
                theta_p = SPM_c_p[N_dr+1]/c_max_p
                phi_p = 2*R*T/F*asinh(J_p*F/(k_p*theta_p^alpha_p*c_el^(1-alpha_p)*c_max_p*(1-theta_p)^(1-alpha_p))/2) + OCP_p_spline(theta_p)
                J_n = I/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)
                SPM_c_n[N_dr+1] = SPM_c_n_N_dr_1 + (-J_n/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
                theta_n = SPM_c_n[N_dr+1]/c_max_n
                phi_n = 2*R*T/F*asinh(J_n*F/(k_n*theta_n^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-theta_n)^(1-alpha_n))/2) + OCP_n_spline(theta_n)
                I = P/(phi_p-phi_n)
            end
          
            eta_sei = phi_n - J_n*F*δ_sei_prev*R_sei
            j_sei = -F*k_sei*C_EC_prev*exp(-α_sei*F*eta_sei /R/T)
            δ_sei = δ_sei_prev - dt*j_sei*M_sei/ρ_sei/(2*F)
            C_EC = δ_sei/D_EC*(j_sei/F)+c_bulk_EC
            J_n = I/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)-j_sei/F
            SPM_c_n[N_dr+1] = SPM_c_n_N_dr_1 + (-J_n/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
            theta_n = SPM_c_n[N_dr+1]/c_max_n
            phi_n = 2*R*T/F*asinh(J_n*F/(k_n*theta_n^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-theta_n)^(1-alpha_n))/2) + OCP_n_spline(theta_n)+J_n*F*δ_sei*R_sei
            I = P/(phi_p-phi_n)

            for kk in 1:3
                J_p = -I/(3*eps_p*El_length_p*El_width_p*El_thickness_p*F/R_p)
                SPM_c_p[N_dr+1] = SPM_c_p_N_dr_1 + (-J_p/D_p)*(2*dt*D_p/dr_p+2*dt*D_p/R_p)
                theta_p = SPM_c_p[N_dr+1]/c_max_p
                phi_p = 2*R*T/F*asinh(J_p*F/(k_p*theta_p^alpha_p*c_el^(1-alpha_p)*c_max_p*(1-theta_p)^(1-alpha_p))/2) + OCP_p_spline(theta_p)
                eta_sei = phi_n - J_n*F*δ_sei*R_sei
                j_sei = -F*k_sei*C_EC*exp(-α_sei*F*eta_sei /R/T)
                δ_sei = δ_sei_prev - dt*j_sei*M_sei/ρ_sei/(2*F)
                C_EC = δ_sei/D_EC*(j_sei/F)+c_bulk_EC                   
                J_n = I/(3*eps_n*El_length_n*El_width_n*El_thickness_n*F/R_n)-j_sei/F
                SPM_c_n[N_dr+1] = SPM_c_n_N_dr_1 + (-J_n/D_n)*(2*dt*D_n/dr_n+2*dt*D_n/R_n)
                theta_n = SPM_c_n[N_dr+1]/c_max_n
                phi_n = 2*R*T/F*asinh(J_n*F/(k_n*theta_n^alpha_n*c_el^(1-alpha_n)*c_max_n*(1-theta_n)^(1-alpha_n))/2) + OCP_n_spline(theta_n)+J_n*F*δ_sei*R_sei
                I = P/(phi_p-phi_n)
            end

            eta_sei = phi_n - J_n*F*δ_sei*R_sei
            j_sei = -F*k_sei*C_EC*exp(-α_sei*F*eta_sei /R/T)
            δ_sei = δ_sei_prev - dt*j_sei*M_sei/ρ_sei/(2*F)
            C_EC = δ_sei/D_EC*(j_sei/F)+c_bulk_EC                   
            
           
            SoC = SoC_prev -I*dt/3600/Q  #SoC through current
            LLI = LLI_prev - (j_sei*(3*eps_n*El_length_n*El_width_n*El_thickness_n/R_n)*dt/3600)/Q #LLI

        else
            SoC_final = 3
            Capacity_loss  = 0
            break
        end

        δ_sei_prev = δ_sei
        C_EC_prev = C_EC
        SPM_c_p_prev = SPM_c_p
        SPM_c_n_prev = SPM_c_n
        SoC_prev = SoC 
        LLI_prev = LLI
    end
    
    SoC_final = SoC
    Capacity_loss = LLI*1e4
    if P/I > V_max ||  P/I < V_min
        #     Capacity_loss =  0
        Capacity_loss_error = 1
    end
    if theta_p<theta_soc1_p || theta_p>theta_soc0_p || theta_n>theta_soc1_n || theta_n<theta_soc0_n
        #Capacity_loss =  0
        Capacity_loss_error = 1
        SoC_final = 3
    end
    if SoC_init>SoC_max_init
        Capacity_loss_error = 1
    #    Capacity_loss =  0
    end
    return SoC_final, Capacity_loss, Capacity_loss_error
end



############################

# SoC_init = 0.09
# PowerInput = -3
# SoC_max_init = 1
#
# SOC_final, Capacity_loss = Final_SoC_Caploss(SoC_init,PowerInput,SoC_max_init)
# @show SOC_final
# @show Capacity_loss

#@show Final_SoC_Caploss(0.2,-8,0.99)

N_steps = 10

function SPM_SEI(N_steps)
    SoC_init_array = LinRange(0, 0.1, N_steps)
    PowerInput_array = 1*LinRange(0, Q_E, 18*N_steps)
    SoC_max_init_array = LinRange(0.8, 1, 2*N_steps)
    SoC_final_array = zeros(size(SoC_init_array)[1]*size(PowerInput_array)[1]*size(SoC_max_init_array)[1],6)

    Threads.@threads for i in 1:size(SoC_init_array)[1]
        for k in 1:size(PowerInput_array)[1]
            for j in 1:size(SoC_max_init_array)[1]
                index_ik = size(PowerInput_array)[1]*size(SoC_max_init_array)[1]*(i-1)+(k-1)*size(SoC_max_init_array)[1]+j
                SoC_init = SoC_init_array[i]
                PowerInput = PowerInput_array[k]
                SoC_max_init = SoC_max_init_array[j]
                SoC_final_array[index_ik,1] = SoC_init
                SoC_final_array[index_ik,2] = PowerInput
                SoC_final_array[index_ik,3] = SoC_max_init
                SoC_final_array[index_ik,4:6] .= Final_SoC_Caploss(SoC_init_array[i],PowerInput_array[k],SoC_max_init_array[j])
            
            end
        end
    end
    return SoC_final_array
end

#time_sum(N_steps) = @time SPM_SEI(N_steps)
#time_sum(N_steps)
@time  Output_Data = SPM_SEI(N_steps)


