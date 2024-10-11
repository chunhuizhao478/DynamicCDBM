python3 examples/buried_fault/pythonscripts/preprocess_optimized.py parameters.json

mpirun -np 8 ./dynamic_cdbm-opt -i examples/buried_fault/static_solve/static_solve_combined.i
mpirun -np 8 ./dynamic_cdbm-opt -i examples/buried_fault/dynamic_solve/dynamic_solve_combined.i

: '
Model Input Parameters Conversion
Python/MOOSE file parameter -> User Interface Input parameter

##Geometry Parameters##
xmin -> xmin (m)
xmax -> xmax (m)
ymin -> ymin (m)
ymax -> ymax (m)
zmin -> zmin (m)
zmax -> zmax (m)
lc -> mesh size at far boundary (m)
lc_fault -> mesh size at near fault (m)

##Material Parameters##
lambda_o -> Lamé Constant $\lambda$ (Pa)
shear_modulus_o -> Shear Modulus $\mu$ (Pa)
rho -> Density $\rho$ (kg/m^3)

##Continuum Damage-Breakage Model Parameters##
xi_0 -> Strain invariant ratio at onset of damage $\xi_0$
xi_d -> Strain invariant ratio at onset of breakage $\xi_d$
Cd_constant -> Damage accumulation rate $C_d$ (1/s)
CdCb_multiplier -> Breakage accumulation rate multiplier Cm ($C_B$ (1/s) = Cm $C_d$)
CBH_constant -> Breakage healing rate $C_{BH}$ (1/s)
C_1 -> Damage healing rate $C_1$ (1/s)
C_2 -> Damage healing rate $C_2$ (1/s)
beta_width -> Width of transitional region $\beta$ 
C_g -> Compliance of fluidity of the fine grain material $C_g$ (1/(Pa s))
m1 -> Coefficient of power law index $m_1$
m2 -> Coefficient of power law index $m_2$
chi -> Ratio of two energy state $\chi$

##Initial Damage Parameters##
nucl_center_x -> Nucleation center x coordinate (m)
nucl_center_y -> Nucleation center y coordinate (m)
nucl_center_z -> Nucleation center z coordinate (m)
len_fault_strike -> Length of initial damage strip along strike direction (m)
len_fault_dip -> Length of initial damage strip along dip direction (m)
len_fault_normal -> Length of initial damage strip along normal direction (m)
nucl_distance -> Length of nucleation patch along strike/dip direction (m)
nucl_thickness -> Lenght of nucleation patch along normal direction (m) 
nucl_damage -> Value of initial damage strip
e_damage -> Value of additional nucleation damage
e_sigma -> Value of initial damage expontial decay
duration -> Simulation time (s)

##Boundary Condition##
normal_traction_x -> Normal boundary traction along x direction (Pa)
normal_traction_z -> Normal boundary traction along z direction (Pa)
normal_traction_y -> Normal boundary traction along y direction (Pa)
shear_traction 55e6 -> Shear boundary traction along x-z direction (Pa)

##Simulation Parameters##
dt -> Simulation time step (s)
end_time -> Total simulation time (s)
time_step_interval -> Result saved time interval
'