using DataFrames
using DelimitedFiles
using Gurobi
using JuMP
using StatsPlots
using XLSX

N_hours = 24

SoC_max = 0.9


day1 = 130
day_last = 365

SoC_max_tracker =  0.9899778736360052


function build_1day(N_hours,day1)
    
   
    Neurons_LLI = 14
    Neurons_Status = 12
    directory = 
    directory_for_parameters =
    file_with_system_parameters = directory_for_parameters*"System_level_data.xlsx"
    LIBESS_Parameters = DataFrame(XLSX.readtable(file_with_system_parameters, "SystemLevel"))

    #System level parameters
    Ch_max = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Maximum charging power","Value"][1]
    Dis_max = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Maximum discharging power","Value"][1]
    E_nom = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Energy capacity","Value"][1]
    SoC_init = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Initial SoC","Value"][1]/100
    SoC_final = LIBESS_Parameters[LIBESS_Parameters[:,"Parameter"].=="Final SoC","Value"][1]/100

    P_max = 18.2*0.75
    P_min = 18.2*0.1
    SoC_operation_max = 0.90
    SoC_operation_min = 0.1

    #Neural Network parameters
    eps_w = 1e-3
   
    directory_for_NN_LLI = "D:\\OneDrive - University of Calgary\\Work\\NN_LIBESS\\data_generation_SPM_SEI\\Data_2\\14N2L_700E_0.001LR_LLI\\"
    w1_LLI = readdlm(directory_for_NN_LLI*"w1.csv", ',', Float64)
    w2_LLI = readdlm(directory_for_NN_LLI*"w2.csv", ',', Float64)
    w3_LLI = readdlm(directory_for_NN_LLI*"w3.csv", ',', Float64)
    w1_LLI[abs.(w1_LLI) .< eps_w] .= 0
    w2_LLI[abs.(w2_LLI) .< eps_w] .= 0
    w3_LLI[abs.(w3_LLI) .< eps_w] .= 0
    b1_LLI = readdlm(directory_for_NN_LLI*"b1.csv", ',', Float64)
    b2_LLI = readdlm(directory_for_NN_LLI*"b2.csv", ',', Float64)
    b3_LLI = readdlm(directory_for_NN_LLI*"b3.csv", ',', Float64)
    M1_min_LLI = 1*readdlm(directory_for_NN_LLI*"M1_min.csv", ',', Float64).-0
    M1_max_LLI = 1*readdlm(directory_for_NN_LLI*"M1_max.csv", ',', Float64).+ 0
    M2_min_LLI = 1*readdlm(directory_for_NN_LLI*"M2_min.csv", ',', Float64).- 0
    M2_max_LLI = 1*readdlm(directory_for_NN_LLI*"M2_max.csv", ',', Float64).+ 0

    directory_for_NN_Status = "D:\\OneDrive - University of Calgary\\Work\\NN_LIBESS\\data_generation_SPM_SEI\\Data_2\\12N2L_700E_0.001LR_Status\\"
    w1_Status = readdlm(directory_for_NN_Status*"w1.csv", ',', Float64)
    w2_Status = readdlm(directory_for_NN_Status*"w2.csv", ',', Float64)
    w3_Status = readdlm(directory_for_NN_Status*"w3.csv", ',', Float64)
    w1_Status[abs.(w1_Status) .< eps_w] .= 0
    w2_Status[abs.(w2_Status) .< eps_w] .= 0
    w3_Status[abs.(w3_Status) .< eps_w] .= 0
    b1_Status = readdlm(directory_for_NN_Status*"b1.csv", ',', Float64)
    b2_Status = readdlm(directory_for_NN_Status*"b2.csv", ',', Float64)
    b3_Status = readdlm(directory_for_NN_Status*"b3.csv", ',', Float64)
    b1_Status[abs.(b1_Status) .< eps_w] .= 0
    b2_Status[abs.(b2_Status) .< eps_w] .= 0
    b3_Status[abs.(b3_Status) .< eps_w] .= 0
    M1_min_Status = 1*readdlm(directory_for_NN_Status*"M1_min.csv", ',', Float64).-0
    M1_max_Status = 1*readdlm(directory_for_NN_Status*"M1_max.csv", ',', Float64).+ 0
    M2_min_Status = 1*readdlm(directory_for_NN_Status*"M2_min.csv", ',', Float64).- 0
    M2_max_Status = 1*readdlm(directory_for_NN_Status*"M2_max.csv", ',', Float64).+ 0
    M1_min_Status[abs.(M1_min_Status) .< eps_w] .= 0
    M1_max_Status[abs.(M1_max_Status) .< eps_w] .= 0
    M1_max_Status[abs.(M2_min_Status) .< eps_w] .= 0
    M2_max_Status[abs.(M2_max_Status) .< eps_w] .= 0


   
    #####################################
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "LogFile", directory*"YEAR_1_NN_LLI_ONLY_Start_DAY__"*string(day1)*"____log_file.txt")
    set_optimizer_attribute(model, "MIPGap", 0.05)
    set_optimizer_attribute(model, "TimeLimit", 900)
    set_optimizer_attribute(model,"NumericFocus", 0)
    set_optimizer_attribute(model,"Threads", 12)
    set_optimizer_attribute(model,"Method", 2)  #Barrier 2 
    set_optimizer_attribute(model,"Heuristics", 0.02) 
    set_optimizer_attribute(model,"Cuts", 3) 
    set_optimizer_attribute(model,"MIPFocus", 3)   
    
    @variable(model, SoC[1:N_hours]>=0)
    @variable(model, P[1:N_hours])
    for i in 1:N_hours
        MOI.set(model, Gurobi.VariableAttribute("Start"), P[i], 0)
    end
       
    @variable(model, Dis[1:N_hours]>=0)
    @variable(model,Ch[1:N_hours]>=0)
    @variable(model, u[1:N_hours],Bin)
    @variable(model, uch[1:N_hours],Bin)
    @variable(model, udis[1:N_hours],Bin)
    @variable(model, SoH>=0)
    @variable(model,LLI[1:N_hours]>=0)
    @variable(model,Status[1:N_hours])


    @variable(model, Z1_LLI[1:N_hours,1:Neurons_LLI]>=0)
    @variable(model, y1_LLI[1:N_hours,1:Neurons_LLI],Bin)
    @variable(model, Z2_LLI[1:N_hours,1:Neurons_LLI]>=0)
    @variable(model, y2_LLI[1:N_hours,1:Neurons_LLI],Bin)


    @variable(model, Z1_Status[1:N_hours,1:Neurons_Status]>=0)
    @variable(model, y1_Status[1:N_hours,1:Neurons_Status],Bin)
    @variable(model, Z2_Status[1:N_hours,1:Neurons_Status]>=0)
    @variable(model, y2_Status[1:N_hours,1:Neurons_Status],Bin)

    ############### Fixing some binary variables
    h_set = [1:N_hours...]

    # Tuned for 12N2L_700E_0.001LR_Status
    y1_Status_fixed_0 = [1:12...]
    y1_Status_fixed_1 = []

    for i1 in h_set
        for i2 in y1_Status_fixed_0
            fix(y1_Status[i1,i2], 0)
        end
        for i2 in y1_Status_fixed_1
            fix(y1_Status[i1,i2], 1)
        end
    end

    y2_Status_fixed_0 = [1:5...,7:12...]
    y2_Status_fixed_1 = [6]

    for i1 in h_set
        for i2 in y2_Status_fixed_0
            fix(y2_Status[i1,i2], 0)
        end
        for i2 in y2_Status_fixed_1
            fix(y2_Status[i1,i2], 1)
        end
    end
    #############################################################

    # Tuned for 6N2L_300E_0.001LR_LLI
    # for i1 in h_set
    #     fix(y2_LLI[i1,5], 0)
    # end


    # Tuned for 6N2L_275_700E_0005LR_SoC_no3_SoC_restricted
    # for i1 in h_set
    #     fix(y1_SoC[i1,3], 1)
    #     fix(y2_SoC[i1,2], 1)
    # end


    ######################  Extra Cuts 
    # Tuned for 6N2L_300E_0.001LR_LLI

    # @constraint(model,Cut_LLI_1[h = 1:N_hours], y1_LLI[h,1]+y1_LLI[h,2]<=1)
    # @constraint(model,Cut_LLI_2[h = 1:N_hours], y1_LLI[h,3]+y1_LLI[h,4]>=1)
    # @constraint(model,Cut_LLI_3[h = 1:N_hours], y1_LLI[h,3]+y1_LLI[h,5]>=1)
    # @constraint(model,Cut_LLI_4[h = 1:N_hours], y2_LLI[h,2]+y2_LLI[h,3]>=1)

    ###############################################


    
   
    
    @constraint(model,Constr_State_of_Charge1[h = [N_hours]],  SoC[h] >= SoC_final) 
    @constraint(model,Constr_State_of_Charge2_0[h = [1]],  SoC[h] == SoC_init + 0.95*Ch[h]/18.2-Dis[h]/18.2)
    @constraint(model,Constr_State_of_Charge2[h = 2:N_hours],  SoC[h] == SoC[h-1]+0.95*Ch[h]/18.2-Dis[h]/18.2)
    @constraint(model,Constr_State_of_Charge3[h = 1:N_hours],  SoC[h] <= SoC_operation_max)
    @constraint(model,Constr_State_of_Charge4[h = 1:N_hours],  SoC[h] >= SoC_operation_min)
    @constraint(model,Constr_State_of_Charge5[h = 1:N_hours],  SoC[h] <= SoH)
    
   


    @constraint(model,Constr_for_Z1_1_LLI[h = [1], k =1:Neurons_LLI], Z1_LLI[h,k] >= w1_LLI[1,k]*SoC_init+w1_LLI[2,k]*P[h]+w1_LLI[3,k]*SoH+b1_LLI[k])
    @constraint(model,Constr_for_Z1_2_LLI[h = [1], k =1:Neurons_LLI], Z1_LLI[h,k] <= w1_LLI[1,k]*SoC_init+w1_LLI[2,k]*P[h]+w1_LLI[3,k]*SoH+b1_LLI[k]-M1_min_LLI[k]*(1-y1_LLI[h,k]))
    @constraint(model,Constr_for_Z1_3_LLI[h = 2:N_hours, k =1:Neurons_LLI], Z1_LLI[h,k] >= w1_LLI[1,k]*SoC[h-1]+w1_LLI[2,k]*P[h]+w1_LLI[3,k]*SoH+b1_LLI[k])
    @constraint(model,Constr_for_Z1_4_LLI[h = 2:N_hours, k =1:Neurons_LLI], Z1_LLI[h,k] <= w1_LLI[1,k]*SoC[h-1]+w1_LLI[2,k]*P[h]+w1_LLI[3,k]*SoH+b1_LLI[k]-M1_min_LLI[k]*(1-y1_LLI[h,k]))
    @constraint(model,Constr_for_Z1_5_LLI[h = 1:N_hours, k =1:Neurons_LLI], Z1_LLI[h,k] <= M1_max_LLI[k]*y1_LLI[h,k])
    @constraint(model,Constr_for_Z2_1_LLI[h = 1:N_hours, k =1:Neurons_LLI], Z2_LLI[h,k] >= sum(w2_LLI[k1,k]*Z1_LLI[h,k1] for k1 in 1:Neurons_LLI)+b2_LLI[k])
    @constraint(model,Constr_for_Z2_2_LLI[h = 1:N_hours, k =1:Neurons_LLI], Z2_LLI[h,k] <= sum(w2_LLI[k1,k]*Z1_LLI[h,k1] for k1 in 1:Neurons_LLI)+b2_LLI[k]-M2_min_LLI[k]*(1-y2_LLI[h,k]))
    @constraint(model,Constr_for_Z2_3_LLI[h = 1:N_hours, k =1:Neurons_LLI], Z2_LLI[h,k] <= M2_max_LLI[k]*y2_LLI[h,k])

    @constraint(model,Constr_LLI_1[h = 1:N_hours],  LLI[h] == sum(w3_LLI[k1,1]*Z2_LLI[h,k1] for k1 in 1:Neurons_LLI)+b3_LLI[1])
    @constraint(model,Constr_LLI_2[h = 1:N_hours],  LLI[h] <= 1.2)
    @constraint(model,Constr_LLI_3[h = 1:N_hours], LLI[h] >= 0.001)


    @constraint(model,Constr_for_Z1_1_Status[h = [1], k =1:Neurons_Status], Z1_Status[h,k] >= w1_Status[1,k]*SoC_init+w1_Status[2,k]*P[h]+w1_Status[3,k]*SoH+b1_Status[k])
    @constraint(model,Constr_for_Z1_2_Status[h = [1], k =1:Neurons_Status], Z1_Status[h,k] <= w1_Status[1,k]*SoC_init+w1_Status[2,k]*P[h]+w1_Status[3,k]*SoH+b1_Status[k]-M1_min_Status[k]*(1-y1_Status[h,k]))
    @constraint(model,Constr_for_Z1_3_Status[h = 2:N_hours, k =1:Neurons_Status], Z1_Status[h,k] >= w1_Status[1,k]*SoC[h-1]+w1_Status[2,k]*P[h]+w1_Status[3,k]*SoH+b1_Status[k])
    @constraint(model,Constr_for_Z1_4_Status[h = 2:N_hours, k =1:Neurons_Status], Z1_Status[h,k] <= w1_Status[1,k]*SoC[h-1]+w1_Status[2,k]*P[h]+w1_Status[3,k]*SoH+b1_Status[k]-M1_min_Status[k]*(1-y1_Status[h,k]))
    @constraint(model,Constr_for_Z1_5_Status[h = 1:N_hours, k =1:Neurons_Status], Z1_Status[h,k] <= M1_max_Status[k]*y1_Status[h,k])
    @constraint(model,Constr_for_Z2_6_Status[h = 1:N_hours, k =1:Neurons_Status], Z2_Status[h,k] >= sum(w2_Status[k1,k]*Z1_Status[h,k1] for k1 in 1:Neurons_Status)+b2_Status[k])
    @constraint(model,Constr_for_Z2_7_Status[h = 1:N_hours, k =1:Neurons_Status], Z2_Status[h,k] <= sum(w2_Status[k1,k]*Z1_Status[h,k1] for k1 in 1:Neurons_Status)+b2_Status[k]-M2_min_Status[k]*(1-y2_Status[h,k]))
    @constraint(model,Constr_for_Z2_8_Status[h = 1:N_hours, k =1:Neurons_Status], Z2_Status[h,k] <= M2_max_Status[k]*y2_Status[h,k])

    @constraint(model,Constr_Status_1[h = 1:N_hours],  Status[h] == sum(w3_Status[k1,1]*Z2_Status[h,k1] for k1 in 1:Neurons_Status)+b3_Status[1])
    @constraint(model,Constr_Status_2[h = 1:N_hours],  Status[h] <= 0.5)

  

    #Constarints to remove outliers and impact of simultaneous charging/discharging
    @constraint(model,Constr_Power_Additional_1[h = 1:N_hours],  Ch[h] <= uch[h]*P_max)
    @constraint(model,Constr_Power_Additional_2[h = 1:N_hours],  Ch[h] >= uch[h]*P_min)
    @constraint(model,Constr_Power_Additional_3[h = 1:N_hours],  Dis[h] <= udis[h]*P_max)
    @constraint(model,Constr_Power_Additional_4[h = 1:N_hours],  Dis[h] >= udis[h]*P_min)
    @constraint(model,Constr_Power_Additional_5[h = 1:N_hours],  P[h] == Dis[h]-Ch[h])
    @constraint(model,Constr_Power_Additional_6[h = 1:N_hours],  Dis[h] <= (1-u[h])*P_max)
    @constraint(model,Constr_Power_Additional_7[h = 1:N_hours],  Ch[h] <= u[h]*P_max)
    

  
    return model, P, Dis, LLI, h_set, SoC, SoH
    
end    


directory = 
directory_for_market = 
directory_for_parameters = 
file_with_system_parameters 
file_with_market_data = 
Data_market = DataFrame(XLSX.readtable(file_with_market_data, "2022"))
LIBESS_Parameters = DataFrame(XLSX.readtable(file_with_system_parameters, "SystemLevel"))



Results_in_array_Dispatch = zeros(N_hours,366)
Results_in_array_SoC = zeros(N_hours,366)
Results_in_array_Obj = zeros(366,9)
Results_in_array_Dispatch[1:N_hours,1] .= LinRange(0, 23, N_hours)
Results_in_array_SoC[1:N_hours,1] .= LinRange(0, 23, N_hours)

model, P, Dis, LLI, h_set, SoC, SoH = build_1day(N_hours,day1)


day_loop = 0
for day in day1:day_last
    global SoC_max_tracker
    global day_loop
    day_loop = day
    fix(SoH, SoC_max_tracker; force=true)
    @objective(model, Max, sum(P[h]*Data_market[h+24*(day-1),"Pool Price"] - LLI[h]*1420-Dis[h]*2  for h in 1:N_hours))

    println("________________________________________")
    println("DAY: ---> ",day)
    println("________________________________________")

    optimize!(model)


    Results_in_array_Dispatch[1:N_hours,day+1] = JuMP.value.(P)
    Results_in_array_SoC[1:N_hours,day+1] = JuMP.value.(SoC)

    SoC_max_tracker = SoC_max_tracker - 1e-4 * JuMP.value(sum(LLI[h]  for h in 1:N_hours))
    
    Results_in_array_Obj[day,1] = day
    Results_in_array_Obj[day,2] = JuMP.objective_value(model)
    Results_in_array_Obj[day,3] = JuMP.value(sum(P[h]*Data_market[h+24*(day-1),"Pool Price"] for h in 1:N_hours))
    Results_in_array_Obj[day,4] = JuMP.value(sum(Dis[h]*2  for h in 1:N_hours))
    Results_in_array_Obj[day,5] = JuMP.value(sum(LLI[h]*1420  for h in 1:N_hours))
    Results_in_array_Obj[day,6] = JuMP.value(sum(LLI[h]*1e-4  for h in 1:N_hours))
    Results_in_array_Obj[day,7] = SoC_max_tracker
    Results_in_array_Obj[day,8] = JuMP.relative_gap(model)*100
    Results_in_array_Obj[day,9] = solve_time(model)

    println()
    println("Objective value: ---> ", JuMP.objective_value(model))
    println("Market Revenue: ---> ", JuMP.value(sum(P[h]*Data_market[h+24*(day-1),"Pool Price"] for h in 1:N_hours)))
    println("Cost of degradation: ---> ", JuMP.value(sum(LLI[h]*1420  for h in 1:N_hours)))
    println("Current SoH: ---> ", SoC_max_tracker)
    println("Time to get a solution [s]: ---> ", solve_time(model))
    
    if SoC_max_tracker <= SoC_max-0.01
        break
    end
end



Results_in_df_disp = DataFrame(Results_in_array_Dispatch,:auto)
Results_in_df_SoC = DataFrame(Results_in_array_SoC,:auto)
Results_in_df_Obj = DataFrame(Results_in_array_Obj,[:Day,:Profit,:Revenue,:Operation_Cost,:Degradation_Cost,:Degradation_amount,:SoC_max_tracker ,:MIP_Gap, :Solution_Time])


Output_name = "__DAYS___"*string(day1)*"-"*string(day_loop)


Results_output_file_name_disp = "YEAR_1_DISPATCH_NN_LLI_ONLY_year2022"*Output_name
if isfile(directory*Results_output_file_name_disp*".xlsx") == true
   rm(directory*Results_output_file_name_disp*".xlsx", force=true)
end
XLSX.writetable(directory*Results_output_file_name_disp*".xlsx", RESULTS=(collect(DataFrames.eachcol(Results_in_df_disp)), DataFrames.names(Results_in_df_disp) ))

Results_output_file_name_SoC = "YEAR_1_SoC_NN_LLI_ONLY_year2022"*Output_name
if isfile(directory*Results_output_file_name_SoC*".xlsx") == true
   rm(directory*Results_output_file_name_SoC*".xlsx", force=true)
end
XLSX.writetable(directory*Results_output_file_name_SoC*".xlsx", RESULTS=(collect(DataFrames.eachcol(Results_in_df_SoC)), DataFrames.names(Results_in_df_SoC) ))


Results_output_file_name_obj = "YEAR_1_OBJECTIVE_NN_LLI_ONLY_year2022"*Output_name
if isfile(directory*Results_output_file_name_obj*".xlsx") == true
   rm(directory*Results_output_file_name_obj*".xlsx", force=true)
end
XLSX.writetable(directory*Results_output_file_name_obj*".xlsx", RESULTS=(collect(DataFrames.eachcol(Results_in_df_Obj)), DataFrames.names(Results_in_df_Obj) ))
