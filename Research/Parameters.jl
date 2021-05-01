# Parameter Settings
const N_homes = Int64(1000)
const connection_density = Int64(4)
const K = Int64(4000) 
const sim_period = Float64(26)
const lockout_period = Float64(15/60)
const DR_period = Float64(20/60)
const H = Int64(20)

# 0 = Constant outside temperature, 1 = Variable outside temperature (TMY)
const disturbance_type = Int64(1) 
# 0 = No Demand Response, 1 = Active Demand Response
const DR_flag = Int64(1)
# 0 = Binary Control, 1 = Continuous Control
const control_type = Int64(0)

# Boundary constraint parameters.
const outside_temp = 95
const init_temp_mean = 68.0         
const init_temp_std = 0
const init_temp_dot_mean = 0.0
const init_temp_dot_std = 0
const T_sp_mean = 70.0
const T_sp_std = 0.25
const db_mean = 2.0
const db_std = 0.025

# House thermal parameters.
const std_percent = 0.05              
const C_A_mean = 1080.0
const C_A_std = C_A_mean * std_percent
const C_M_mean = 4280.0
const C_M_std = C_M_mean * std_percent
const U_A_mean = 520.0
const U_A_std = U_A_mean * std_percent
const H_M_mean = 7050.0
const H_M_std = H_M_mean * std_percent

# Time conditions
const ti = 266.5
const tf = ti + sim_period
const Δt = abs(tf - ti) / K
const t_dr_i = abs((tf + ti) / 2) - DR_period / 2
const t_dr_f = abs((tf + ti) / 2) + DR_period / 2
