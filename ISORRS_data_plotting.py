# VARIABLE PLOTTING TOOL written by H. S. Joyce, 2022. 

# This module is designed to be able to graphically show data from the 
# input and output files produced by the ISORRS model

# BELOW is data for what information each line of both types of file contains
# so that the code appropriate information can be extracted for what you wish 
# to visualise. The input file column and output file column are for the 
# input data type and output data type respectively

# NEED TO UPD8 LINE NUMBER BC WIDTHS

# LINE NUMBER ------ INPUT FILE  ----------------- OUTPUT FILE
#-----------------------------------------------------------------
# 1 ---------------- 'setup' ------------------- 'results'
# 2 --------------- 'planet = ' ----- 'TPS - electron total particle source'
# 3 ---------------'timestep = ' ------ 'TPS' -  ion total particle source'
# 4 -------------- 'spatial step = ' ------- 'TMS - total mass source'
# 5 ---------------- 'iterations = ' ------------ 'setup'
# 6 ------------- 'outer limit = (RJ)' ---------- 'planet'
# 7 --------------- 'l shell' ------------------- 'timestep'
# 8 --------------- 'FAC on?' ------------------ 'spatial step'
# 9 ------------- 'FAC Strength = ' ------------ 'iterations = '
# 10 --------------- 'CF on?' --------------- 'outer limit = (RJ)'
# 11 --------------- 'Lat Width =' ---------- 'l shell'
# 12 --------------- 'Long Width =' --------- 'FAC on?'
# 12 ------------ 'Colatitude = ' --------------- 'FAC on?'
# 13 --- 'Initial Ionospheric Temperature = ' -- 'FAC Strength = '
# 14 ---------- 'Number of Ions = ' -------------- 'CF on?'
# 15 -------- 'Number of Neutrals = ' -------- 'Lat Width = '
# 16 -------- ---- Ion Species 1: ' ---------- 'Long Width = '
# 16 ------------- Ion Species 1: ' ------------ 'Colatitude'
# 17 ------------- Ion Species 2: ' ----- 'Initial Ionospheric Temperature'
# 18 ----------- Neutral Species 1: ' --------- 'Number of Ions = '
# 19 ----------- Neutral Species 2: ' ------- 'Number of Neutrals = '
# 20 ----------- Neutral Species 3: ' ---------- 'Ion Species 1  '
# 21 -------------- 'Vectors' ------------------ 'Ion Species 2: '
# 22 ---------- 'Current Density' -------------- 'Neutral Species 1: '
# 23 -------- current density values ----------- 'Neutral Species 2: '
# 24 -------- 'Cross Sectional Area' ----------- 'Neutral Species 3: '
# 25 ------------ area values ---------------------- 'Vectors'
# 26 -------------- 'Grid' --------------------- 'Current Density'
# 27 ----------- grid values ----------------- current density values 
# 28 ------ 'Gravitational Acceleration' ------ 'Cross Sectional Area'
# 29 --- gravitational acceleration values --------- area values
# 30 ------- 'Centrifugal Acceleration' -------------- 'Grid'
# 31 ---- centrifugal acceleration values ---------- grid values
# 32 ------------ 'Electrons' -------------- 'Gravitational Acceleration'
# 33 ------------ 'rho' for e ------------ gravitational acceleration values
# 34 ----------- rho e values --------------- 'Centrifugal Acceleration'
# 35 ------------- 'n' for e ------------- centrifugal acceleration values
# 36 ------------- n e values ---------------------- 'E Field'
# 37 ------------- 'u' for e -------------------- E Field values
# 38 ------------- u e values -------------------- rho e values
# 39 ------------- 'P' for e --------------------- 'rho' for e
# 40 ------------- P e values -------------------- rho e values
# 41 ------------- 'T' for e ----------------------- 'n' for e
# 42 ------------- T e values ---------------------- n e values
# 43 ----------- 'kappa' for e --------------------- 'u' for e
# 44 ----------- kappa e values -------------------- u e values
# 45 --------------- 'Ion 1' ----------------------- 'P' for e
# 46 ------------ 'rho' for ion 1 ------------------ P e values
# 47 ----------- rho ion 1 values ------------------ 'T' for e
# 48 ------------- 'n' for ion 1 ------------------- T e values
# 49 ------------ n ion 1 values ----------------- 'kappa' for e
# 50 ------------- 'u' for ion 1 ----------------- kappa e values
# 51 ------------- u ion 1 values -------------------- 'Ion 1'
# 52 ------------- 'P' for ion 1 ----------------- 'rho' for ion 1
# 53 ------------- P ion 1 values ---------------- rho ion 1 values
# 54 ------------ 'T' for ion 1 ------------------- 'n' for ion 1
# 55 ------------ T ion 1 values ------------------ n ion 1 values
# 56 ---------- 'kappa' for ion 1 ----------------- 'u' for ion 1
# 57 ---------- kappa ion 1 values ---------------- u ion 1 values
# 58 --------------- 'Ion 2' ---------------------- 'P' for ion 1
# 59 ----------- 'rho' for ion 2 ------------------ P ion 1 values
# 60 ----------- rho ion 2 values ----------------- 'T' for ion 1
# 61 ------------ 'n' for ion 2 ------------------- T ion 1 values
# 62 ------------ n ion 2 values ---------------- 'kappa' for ion 1
# 63 ------------ 'u' for ion 2 ----------------- kappa ion 1 values
# 64 ------------ u ion 2 values --------------------- 'Ion 2'
# 65 ------------ 'P' for ion 2 ------------------ 'rho' for ion 2
# 66 ------------ P ion 2 values ----------------- rho ion 2 values
# 67 ------------ 'T' for ion 2 ------------------- 'n' for ion 2
# 68 ------------ T ion 2 values ------------------ n ion 2 values
# 69 ---------- 'kappa' for ion 2 ----------------- 'u' for ion 2
# 70 ---------- kappa ion 2 values ---------------- u ion 2 values
# 71 ---------------------------------------------- 'P' for ion 2
# 72 ---------------------------------------------- P ion 2 values
# 73 ---------------------------------------------- 'T' for ion 2
# 74 ---------------------------------------------- T ion 2 values
# 75 -------------------------------------------- 'kappa' for ion 2
# 76 -------------------------------------------- kappa ion 2 values


# import needed modules
import matplotlib.pyplot as pl
import numpy as np

# establish a 'counter' to move through files line by line
# line_count = 1

# # open file
# run_single = open("/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/ISORRS_input_jupiter_dusk_high_temp_flux_look.txt","r")
# #run_single = open("/Users/hannah/Documents/HEC_outputs/outputs/ISORRS_output_jupiter_test.txt","r")
# # loop to go through the file and pull out specific lines
# for line in run_single:
#     # move through initial un-needed lines
#     if line_count < 27:
#         line_count += 1
#     elif line_count == 27:
#         # remove any excess spacing
#         grid = line.strip()
#         # establish comma as delimmiter between values
#         grid = line.split(",")
#         # explicitly set as an array with numbers as floats
#         grid = np.array(grid[:-1]).astype(float)
#         # convert grid into km
#         grid = (grid/1000)
#         line_count +=1
#         # move to the next line in the file
#     elif line_count == 51:
#         u_H = line.strip()
#         u_H = line.split(",")
#         u_H = np.array(u_H[:-1]).astype(float)
#         #n_H3 = (T_H/1e-7)
#         line_count += 1
#     else:
#         line_count += 1
#         # keep moving through file to end
        
        
# line_count2 =1
# # open file
# run_single_C = open("/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/ISORRS_input_jupiter_auroral_carley_dusk.txt","r")
# #run_single = open("/Users/hannah/Documents/HEC_outputs/outputs/ISORRS_output_jupiter_test.txt","r")
# # loop to go through the file and pull out specific lines
# for line in run_single_C:
#     # move through initial un-needed lines
#     if line_count2 < 27:
#         line_count2 += 1
#     elif line_count2 == 27:
#         # remove any excess spacing
#         grid2 = line.strip()
#         # establish comma as delimmiter between values
#         grid2 = line.split(",")
#         # explicitly set as an array with numbers as floats
#         grid2 = np.array(grid2[:-1]).astype(float)
#         # convert grid into km
#         grid2 = (grid2/1000)
#         line_count2 +=1
#         # move to the next line in the file
#     elif line_count2 == 62:
#         n2_H3 = line.strip()
#         n2_H3 = line.split(",")
#         n2_H3 = np.array(n2_H3[:-1]).astype(float)
#         #n_H3 = (T_H/1e-7)
#         line_count2 += 1
#     else:
#         line_count2 += 1
#         # keep moving through file to end  
        
# dens_jim = [2e-06,5e-09,5e-09,5e-09,33040000000,1382200000,3745000000,8330000000,141110000000,86880000000,55490000000,40270000000,30220000000,22922000000,17200000000,12619000000,9134000000,6665000000,5084000000,4011000000,3314000000,2849000000,2371000000,1982300000,1542500000,1203100000,877700000,627700000,422700000,263500000]
# alts_jim = [714,763.4,822.6,889.8,966.4,1053.5,1152.7,1264.9,1390,1524.8,1665.2,1810.4,1960.2,2114,2269,2427,2585,2745,2905,3067,3228,3391,3554,3718,3884,4050,4216,4384,4553,4723]

# plotting       
# fig = pl.figure(figsize=(10,8))
# pl.subplots_adjust(wspace=0, hspace=0) 
# # save location
# #folder = '/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/HEC_trials/'
# folder = '/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/'
# # establish axes
# ax1 = fig.add_subplot(1,1,1)
# ax1.set_yscale('log')
# # plot data
# ax1.plot(grid,u_H,linestyle='-',color='orange')
# # ax1.plot(grid2,n2_H3,linestyle='-',color='blue')
# # ax1.plot(alts_jim[9:],dens_jim[9:],linestyle='-',color='black')
# # can scale axes if want to zoom in on certain section 
# #ax1.set_xlim(1400, 60000)
# #ax1.xaxis.set_ticks(np.arange(1000, 50000,3500))
# # label axes
# ax1.set_ylabel('T')
# ax1.set_xlabel('Distance Along Field Line (km)')
# ax1.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
# #ax1.legend(['Scale Heights', 'Carleys Function', 'JIM'],loc='upper right',bbox_to_anchor=(0.99,0.95))
# # save file as 'name'
# pl.savefig(folder+'u_H_highT.png')
# pl.show()


line_count=1


run_1 = open("/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/ISORRS_input_jupiter_dusk_density_asym_1.txt","r")
for line in run_1:
    #print(line_count)
    if line_count < 27:
        line_count +=1
        continue
    elif line_count == 27:
        grid  = line.strip()
        grid = line.split(",")
        grid = np.array(grid[:-1]).astype(float)
        grid = (grid/1000)
        line_count +=1
    elif line_count == 36:
        n_e = line.strip()
        #print(n_e)
        n_e = line.split(",")
        n_e = np.array(n_e[:-1]).astype(float)
        line_count +=1
    elif line_count == 38:
        u_e = line.strip()
        u_e = line.split(",")
        u_e = np.array(u_e[:-1]).astype(float)
        line_count += 1
    elif line_count == 49:
        n_H_plus = line.strip()
        #print(n_H_plus)
        # use comma as delimiter
        n_H_plus = line.split(",")
        n_H_plus = np.array(n_H_plus[:-1]).astype(float)
        line_count +=1
    elif line_count == 51:
        u_H_plus = line.strip()
        u_H_plus = line.split(",")
        u_H_plus = np.array(u_H_plus[:-1]).astype(float)
        line_count +=1
    elif line_count == 62:
        n_H3_plus = line.strip()
        # use comma as delimiter
        n_H3_plus = line.split(",")
        n_H3_plus = np.array(n_H3_plus[:-1]).astype(float)
        line_count +=1
    elif line_count == 64:
        u_H3_plus = line.strip()
        u_H3_plus = line.split(",")
        u_H3_plus = np.array(u_H3_plus[:-1]).astype(float)
        line_count +=1
    else:
        line_count +=1
    
        
#print(grid)

line_count_2 = 1

run_2 = open("/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/ISORRS_input_jupiter_dusk_density_asym_2.txt","r")
for line in run_2:
    #print(line_count)
    if line_count_2 < 27:
        line_count_2 +=1
        continue
    elif line_count_2 == 27:
        grid_2  = line.strip()
        grid_2 = line.split(",")
        grid_2 = np.array(grid_2[:-1]).astype(float) 
        line_count_2 +=1
    elif line_count_2 == 36:
        n_e_2 = line.strip()
        n_e_2 = line.split(",")
        n_e_2 = np.array(n_e_2[:-1]).astype(float)
        line_count_2 +=1
    elif line_count_2 == 38:
        u_e_2 = line.strip()
        u_e_2 = line.split(",")
        u_e_2 = np.array(u_e_2[:-1]).astype(float)
        line_count_2 += 1
    elif line_count_2 == 49:
        n_H_plus_2 = line.strip()
        # use comma as delimiter
        n_H_plus_2 = line.split(",")
        n_H_plus_2 = np.array(n_H_plus_2[:-1]).astype(float)
        line_count_2 +=1
    elif line_count_2 == 51:
        u_H_plus_2 = line.strip()
        u_H_plus_2 = line.split(",")
        u_H_plus_2 = np.array(u_H_plus_2[:-1]).astype(float)
        line_count_2 += 1
    elif line_count_2 == 62:
        n_H3_plus_2 = line.strip()
        # use comma as delimiter
        n_H3_plus_2 = line.split(",")
        n_H3_plus_2 = np.array(n_H3_plus_2[:-1]).astype(float)
        line_count_2 +=1
    elif line_count_2 == 64:
        u_H3_plus_2 = line.strip()
        u_H3_plus_2 = line.split(",")
        u_H3_plus_2 = np.array(u_H3_plus_2[:-1]).astype(float)
        line_count_2 += 1
    else:
        line_count_2 +=1
        
line_count_3 = 1
        
run_3 = open("/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/ISORRS_input_jupiter_dusk_density_asym_3.txt","r")
for line in run_3:
    #print(line_count)
    if line_count_3 < 27:
        line_count_3 +=1
        continue
    elif line_count_3 == 27:
        grid_3  = line.strip()
        grid_3 = line.split(",")
        grid_3 = np.array(grid_3[:-1]).astype(float)
        line_count_3 +=1
    elif line_count_3 == 36:
        n_e_3 = line.strip()
        n_e_3 = line.split(",")
        n_e_3 = np.array(n_e_3[:-1]).astype(float)
        line_count_3 +=1
    elif line_count_3 == 38:
        u_e_3 = line.strip()
        u_e_3 = line.split(",")
        u_e_3 = np.array(u_e_3[:-1]).astype(float)
        line_count_3 += 1
    elif line_count_3 == 49:
        n_H_plus_3 = line.strip()
        # use comma as delimiter
        n_H_plus_3 = line.split(",")
        n_H_plus_3 = np.array(n_H_plus_3[:-1]).astype(float)
        line_count_3 +=1
    elif line_count_3 == 51:
        u_H_plus_3 = line.strip()
        u_H_plus_3 = line.split(",")
        u_H_plus_3 = np.array(u_H_plus_3[:-1]).astype(float)
        line_count_3 += 1
    elif line_count_3 == 62:
        n_H3_plus_3 = line.strip()
        # use comma as delimiter
        n_H3_plus_3 = line.split(",")
        n_H3_plus_3 = np.array(n_H3_plus_3[:-1]).astype(float)
        line_count_3 +=1
    elif line_count_3 == 64:
        u_H3_plus_3 = line.strip()
        u_H3_plus_3 = line.split(",")
        u_H3_plus_3 = np.array(u_H3_plus_3[:-1]).astype(float)
        line_count_3 += 1
    else:
        line_count_3 +=1

line_count_4 = 1
        
run_4 = open("/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/ISORRS_input_jupiter_dusk_density_asym_4.txt","r")
for line in run_4:
    #print(line_count)
    if line_count_4 < 27:
        line_count_4 +=1
        continue
    elif line_count_4 == 27:
        grid_4  = line.strip()
        grid_4 = line.split(",")
        grid_4 = np.array(grid_4[:-1]).astype(float)
        line_count_4 +=1
    elif line_count_4 == 36:
        n_e_4 = line.strip()
        n_e_4 = line.split(",")
        n_e_4 = np.array(n_e_4[:-1]).astype(float)
        line_count_4 +=1
    elif line_count_4 == 38:
        u_e_4 = line.strip()
        u_e_4 = line.split(",")
        u_e_4 = np.array(u_e_4[:-1]).astype(float)
        line_count_4 += 1
    elif line_count_4 == 49:
        n_H_plus_4 = line.strip()
        # use comma as delimiter
        n_H_plus_4 = line.split(",")
        n_H_plus_4 = np.array(n_H_plus_4[:-1]).astype(float)
        line_count_4 +=1
    elif line_count_4 == 51:
        u_H_plus_4 = line.strip()
        u_H_plus_4 = line.split(",")
        u_H_plus_4 = np.array(u_H_plus_4[:-1]).astype(float)
        line_count_4 += 1
    elif line_count_4 == 62:
        n_H3_plus_4 = line.strip()
        # use comma as delimiter
        n_H3_plus_4 = line.split(",")
        n_H3_plus_4 = np.array(n_H3_plus_4[:-1]).astype(float)
        line_count_4 +=1
    elif line_count_4 == 64:
        u_H3_plus_4 = line.strip()
        u_H3_plus_4 = line.split(",")
        u_H3_plus_4 = np.array(u_H3_plus_4[:-1]).astype(float)
        line_count_4 += 1
    else:
        line_count_4 +=1
        
  
line_count_5 = 1
    
run_5 = open("/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/ISORRS_input_jupiter_dusk_density_asym_5.txt","r")
for line in run_5:
    #print(line_count)
    if line_count_5 < 27:
        line_count_5 +=1
        continue
    elif line_count_5 == 27:
        grid_5  = line.strip()
        grid_5 = line.split(",")
        grid_5 = np.array(grid_5[:-1]).astype(float)
        line_count_5 +=1
    elif line_count_5 == 36:
        n_e_5 = line.strip()
        n_e_5 = line.split(",")
        n_e_5 = np.array(n_e_5[:-1]).astype(float)
        line_count_5 +=1
    elif line_count_5 == 38:
        u_e_5 = line.strip()
        u_e_5 = line.split(",")
        u_e_5 = np.array(u_e_5[:-1]).astype(float)
        line_count_5 += 1
    elif line_count_5 == 49:
        n_H_plus_5 = line.strip()
        # use comma as delimiter
        n_H_plus_5 = line.split(",")
        n_H_plus_5 = np.array(n_H_plus_5[:-1]).astype(float)
        line_count_5 +=1
    elif line_count_5 == 51:
        u_H_plus_5 = line.strip()
        u_H_plus_5 = line.split(",")
        u_H_plus_5 = np.array(u_H_plus_5[:-1]).astype(float)
        line_count_5 += 1
    elif line_count_5 == 62:
        n_H3_plus_5 = line.strip()
        # use comma as delimiter
        n_H3_plus_5 = line.split(",")
        n_H3_plus_5 = np.array(n_H3_plus_5[:-1]).astype(float)
        line_count_5 +=1
    elif line_count_5 == 64:
        u_H3_plus_5 = line.strip()
        u_H3_plus_5 = line.split(",")
        u_H3_plus_5 = np.array(u_H3_plus_5[:-1]).astype(float)
        line_count_5 += 1
    else:
        line_count_5 +=1
        
 
flux_e = np.multiply(n_e,u_e)
flux_e_2 = np.multiply(n_e_2,u_e_2)
flux_e_3 = np.multiply(n_e_3,u_e_3)
flux_e_4 = np.multiply(n_e_4,u_e_4)
flux_e_5 = np.multiply(n_e_5,u_e_5)
print(flux_e)

flux_H_plus = np.multiply(n_H_plus,u_H_plus)
flux_H_plus_2 = np.multiply(n_H_plus_2,u_H_plus_2)
flux_H_plus_3 = np.multiply(n_H_plus_3,u_H_plus_3)
flux_H_plus_4 = np.multiply(n_H_plus_4,u_H_plus_4)
flux_H_plus_5 = np.multiply(n_H_plus_5,u_H_plus_5)

flux_H3_plus = np.multiply(n_H3_plus,u_H3_plus)
flux_H3_plus_2 = np.multiply(n_H3_plus_2,u_H3_plus_2)
flux_H3_plus_3 = np.multiply(n_H3_plus_3,u_H3_plus_3)
flux_H3_plus_4 = np.multiply(n_H3_plus_4,u_H3_plus_4)
flux_H3_plus_5 = np.multiply(n_H3_plus_5,u_H3_plus_5)
#from scipy.signal import savgol_filter

fig = pl.figure(figsize=(7,7))
pl.subplots_adjust(wspace=0, hspace=0) 
folder = '/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/plots/'
ax1 = fig.add_subplot(3,1,1)
ax1.set_yscale('log')
#ax1.set_yscale('log')
# yhat = savgol_filter(n_H3_plus, 51, 3)
# yhat2 = savgol_filter(n_H3_plus_2, 51, 3)
# yhat3 = savgol_filter(n_H3_plus_3, 51, 3)
# yhat4 = savgol_filter(n_H3_plus_4, 51, 3)
# yhat5 = savgol_filter(n_H3_plus_5, 51, 3) # window size 51, polynomial order 3
ax1.plot(grid,flux_e,linestyle='-',color='orange')
ax1.plot(grid,flux_e_2,linestyle='-',color='blue')
ax1.plot(grid,flux_e_3,linestyle='-',color='green')
ax1.plot(grid,flux_e_4,linestyle='-',color='red')
ax1.plot(grid,flux_e_5,linestyle='-',color='black')
ax1.set_xlim(1400, 30000)
#formatter = ax1.get_major_formatter()
#ax1.set_minor_formatter(formatter)
#ax1.xaxis.set_ticks(np.arange(1000, 50000,3500))
ax1.set_ylabel('Electron Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
#pl.xlim([0,z_ext[-1]/1000])
#ax1.locator_params(axis='y', nbins=6)
#ax1.set_xlabel('Distance Along Field Line \n (km)')
#ax1.xaxis.set_ticklabels([])
ax1.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
ax1.xaxis.set_ticklabels([])
#Ωax1.yaxis.set_major_locator(pl.MaxNLocator(6))

ax2 = fig.add_subplot(3,1,2)
ax2.set_yscale('log')
ax2.plot(grid,flux_H_plus,linestyle='-',color='orange')
ax2.plot(grid,flux_H_plus_2,linestyle='-',color='blue')
ax2.plot(grid,flux_H_plus_3,linestyle='-',color='green')
ax2.plot(grid,flux_H_plus_4,linestyle='-',color='red')
ax2.plot(grid,flux_H_plus_5,linestyle='-',color='black')
ax2.set_xlim(1400, 30000)
#ax2.xaxis.set_ticks(np.arange(1000, 50000,3500))
ax2.set_ylabel('$H^+$ Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
ax2.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
ax2.xaxis.set_ticklabels([])
#pl.grid('on') 

ax3 = fig.add_subplot(3,1,3)
ax3.set_yscale('log')
ax3.plot(grid,flux_H3_plus,linestyle='-',color='orange')
ax3.plot(grid,flux_H3_plus_2,linestyle='-',color='blue')
ax3.plot(grid,flux_H3_plus_3,linestyle='-',color='green')
ax3.plot(grid,flux_H3_plus_4,linestyle='-',color='red')
ax3.plot(grid,flux_H3_plus_5,linestyle='-',color='black')
ax3.set_xlim(1400, 30000)
#ax3.xaxis.set_ticks(np.arange(1400, 75000,7000))
ax3.set_ylabel('$H_3^+$ Flux \n $(m^{-2}$ $s^{-1})$ x A $(m^2)$')
ax3.set_xlabel('Distance Along Field Line \n (km)')
ax3.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
pl.savefig(folder+'flux_compare_zoomed_all.png')
pl.show()


# plotting       
fig = pl.figure(figsize=(5,3))
pl.subplots_adjust(wspace=0, hspace=0) 
# save location
#folder = '/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/HEC_trials/'
folder = '/Users/hannah/OneDrive - Lancaster University/ISORRS/2023/march/asym_runs/plots'
# establish axes
ax1 = fig.add_subplot(1,1,1)
ax1.set_yscale('log')
# plot data
ax1.plot(grid,n_H_plus,linestyle='-',color='orange')
ax1.plot(grid,n_H_plus_2,linestyle='-',color='blue')
ax1.plot(grid,n_H_plus_3,linestyle='-',color='green')
ax1.plot(grid,n_H_plus_4,linestyle='-',color='red')
ax1.plot(grid,n_H_plus_5,linestyle='-',color='black')
# can scale axes if want to zoom in on certain section 
ax1.set_xlim(1400, 30000)
#ax1.xaxis.set_ticks(np.arange(1000, 50000,3500))
# label axes
ax1.set_ylabel('$H^+$ Number Density $(m^{-3})$')
ax1.set_xlabel('Distance Along Field Line (km)')
ax1.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
#ax1.legend(['Scale Heights', 'Carleys Function', 'JIM'],loc='upper right',bbox_to_anchor=(0.99,0.95))
# save file as 'name'
pl.savefig(folder+'n_H_plus_dens.png')
pl.show()
