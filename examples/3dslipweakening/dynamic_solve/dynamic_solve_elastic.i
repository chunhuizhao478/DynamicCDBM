#parameters

##mesh parameters
###need to check gmsh file for changing the parameters here
bottom_nodes_coord =' -60000 -60000 -60000;
                      60000 -60000 -60000;
                      60000 60000  -60000;
                     -60000 60000  -60000'

xmin_fault = -15000 #xmin of fault
xmax_fault = 15000 #xmax of fault
zmin_fault = -15000 #zmin of fault
# zmax_fault = 0 #zmax of fault
elem_size = 100 #!!! element size near the fault, need to be consistent with the mesh file
##-------------------------##

##material properties##
density = 2670 #density
lambda_o = 3.204e10 #first lame constant
shear_modulus_o = 3.204e10 #second lame constant
# Cs = '${fparse shear_modulus_o / density }' #shear wave speed
# Cp = '${fparse (lambda_o + 2 * shear_modulus_o) / density }' #pressure wave speed
##-------------------------##

##Slip weakening parameters##
Dc = 0.4 #characteristic length (m)
q = 0.4 #damping ratio
mu_s = 0.677 #static friction coefficient
mu_d = 0.525 #dynamic friction coefficient

use_cohesion = true #use cohesion
cohesion_expression = 'if(z >= -1000, 20e3 * z + 20e6, 0)'
##-------------------------##

##CDB model parameters##
xi_0 = -0.8 #strain invariants ratio: onset of damage evolution
xi_d = -0.9 #strain invariants ratio: onset of breakage healing

###constant Cd
Cd_constant = 0 #coefficient gives positive damage evolution
###

###rate dependent Cd
###note: if use_strain_rate_dependent_Cd is true, Cd_constant will be ignored
use_strain_rate_dependent_Cd = false #use strain rate dependent Cd
m_exponent = 0.8 #strain rate dependent parameters
strain_rate_hat = 1e-4 #strain rate dependent parameters
cd_hat = 1.0 #strain rate dependent parameters
###

CdCb_multiplier = 100 #multiplier between Cd and Cb
CBH_constant = 1e4 #coefficient of healing for breakage evolution
C_1 = 300 #coefficient of healing for damage evolution
C_2 = 0.05 #coefficient of healing for damage evolution
beta_width = 0.05 #coefficient gives width of transitional region
C_g = 1e-10 #material parameter: compliance or fluidity of the fine grain granular material
m1 = 10 #coefficient of power law indexes
m2 = 1 #coefficient of power law indexes
chi = 0.8 #energy ratio
##-------------------------##

##initial damage parameters
sigma = 5e2
peak_val = 0.7
len_of_fault_strike = 30000
len_of_fault_dip = 15000
fault_center = '0 0 -7500'
##-------------------------##

##initial stress parameters##
peak_shear_value = 81.6e6 #initial shear stress perturbation peak value
nucl_center_x = 0 #nucleation center x coordinate
nucl_center_z = -7500 #nucleation center y coordinate
nucl_size = 3000 #nucleation size
##-------------------------##

##model parameters##
dt = 0.0025 #time step size

end_time = 6 #end time for simulation

# num_steps = 40 #end_time or num_steps only one of them is needed
exodus_time_step_interval = 10 #time step interval for output
csv_time_step_interval = 40 #time step interval for csv output
checkpoint_time_step_interval = 80 #time step interval for checkpoint output
checkpoint_num_files = 2 #number of files for checkpoint output
##-------------------------##

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = '../mesh/tpv2053d_100m.msh'
  []
  [./new_block_1]
    type = ParsedSubdomainMeshGenerator
    input = msh
    combinatorial_geometry = 'x >= ${xmin_fault} & x <= ${xmax_fault} & z >= ${zmin_fault} & y < 0'
    block_id = 100
  []
  [./new_block_2]
    type = ParsedSubdomainMeshGenerator
    input = new_block_1
    combinatorial_geometry = 'x >= ${xmin_fault} & x <= ${xmax_fault} & z >= ${zmin_fault} & y > 0'
    block_id = 200
  []       
  [./split_1]
    type = BreakMeshByBlockGenerator
    input = new_block_2
    split_interface = true
    block_pairs = '100 200'
  []      
  [./sidesets]
    input = split_1
    type = SideSetsFromNormalsGenerator
    normals = '-1 0 0
                1 0 0
                0 -1 0
                0 1 0
                0 0 -1
                0 0 1'
    new_boundary = 'left right bottom top back front'
  [] 
  [./extranodeset1]
      type = ExtraNodesetGenerator
      coord = ${bottom_nodes_coord}
      new_boundary = corner_ptr
      input = sidesets
  []
[]

[GlobalParams]

  ##------------slip weakening------------##
  displacements = 'disp_x disp_y disp_z'
  
  #damping ratio
  q = ${q}

  ##----continuum damage breakage model----##
  #initial lambda value (first lame constant) [Pa]
  lambda_o = ${lambda_o}
  
  #initial shear modulus value (second lame constant) [Pa]
  shear_modulus_o = ${shear_modulus_o}

  #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
  xi_0 = ${xi_0}

  #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
  xi_d = ${xi_d}

  #<strain invariants ratio: maximum allowable value>: set boundary
  #Xu_etal_P15-2D
  #may need a bit space, use 1.5 as boundary
  xi_max = 1.8

  #<strain invariants ratio: minimum allowable value>: set boundary
  #Xu_etal_P15-2D
  xi_min = -1.8

  #if option 2, use Cd_constant
  Cd_constant = ${Cd_constant}

  #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
  CdCb_multiplier = ${CdCb_multiplier}

  #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
  # CBCBH_multiplier = 0.0
  CBH_constant = ${CBH_constant}

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_1 = ${C_1}

  #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
  C_2 = ${C_2}

  #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  beta_width = ${beta_width}

  #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  C_g = ${C_g}

  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
  m1 = ${m1}

  #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
  m2 = ${m2}

  # energy ratio
  chi = ${chi}

  #diffusion coefficient #for structural stress coupling
  D = 0
[]

[AuxVariables]
  [./resid_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./resid_y]
      order = FIRST
      family = LAGRANGE
  []
  [./resid_z]
    order = FIRST
    family = LAGRANGE
  []
  [./resid_slipweakening_x]
      order = FIRST
      family = LAGRANGE
  [../]
  [./resid_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  [../]
  [./resid_slipweakening_z]
      order = FIRST
      family = LAGRANGE
  [../]
  [./disp_slipweakening_x]
      order = FIRST
      family = LAGRANGE
  []
  [./disp_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  []
  [./disp_slipweakening_z]
    order = FIRST
    family = LAGRANGE
  []
  [./vel_slipweakening_x]
    order = FIRST
    family = LAGRANGE
  []
  [./vel_slipweakening_y]
      order = FIRST
      family = LAGRANGE
  []
  [./vel_slipweakening_z]
    order = FIRST
    family = LAGRANGE
  []
  ###
  # [jump_x_aux]
  #   order = FIRST
  #   family = MONOMIAL
  # []
  # [jump_x_rate_aux]
  #   order = FIRST
  #   family = MONOMIAL
  # []
  # [traction_x_aux]
  #   order = FIRST
  #   family = MONOMIAL
  # []  
  ###
  #output CDB model properties
  [alpha_damagedvar_aux]
      order = FIRST
      family = MONOMIAL
  []
  [B_aux]
      order = FIRST
      family = MONOMIAL
  []
  [xi_aux]
      order = FIRST
      family = MONOMIAL
  [] 
[]

[Physics/SolidMechanics/CohesiveZone]
  [./czm_ik]
    boundary = 'Block100_Block200'
    strain = SMALL
    generate_output='traction_x traction_y traction_z jump_x jump_y jump_z'
  [../]
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = true
        generate_output = 'stress_xx stress_yy stress_xy'
        extra_vector_tags = 'restore_tag'
      []
    []
  []
[]

[Problem]
  extra_tag_vectors = 'restore_tag'
[]

[AuxKernels]
  [Displacment_x]
    type = ProjectionAux
    variable = disp_slipweakening_x
    v = disp_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Displacement_y]
    type = ProjectionAux
    variable = disp_slipweakening_y
    v = disp_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Displacement_z]
    type = ProjectionAux
    variable = disp_slipweakening_z
    v = disp_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Vel_x]
    type = CompVarRate
    variable = vel_slipweakening_x
    coupled = disp_x
    execute_on = 'TIMESTEP_END'
  []
  [Vel_y]
    type = CompVarRate
    variable = vel_slipweakening_y
    coupled = disp_y
    execute_on = 'TIMESTEP_END'
  []
  [Vel_z]
    type = CompVarRate
    variable = vel_slipweakening_z
    coupled = disp_z
    execute_on = 'TIMESTEP_END'
  []
  [Residual_x]
    type = ProjectionAux
    variable = resid_slipweakening_x
    v = resid_x
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_y]
    type = ProjectionAux
    variable = resid_slipweakening_y
    v = resid_y
    execute_on = 'TIMESTEP_BEGIN'
  []
  [Residual_z]
    type = ProjectionAux
    variable = resid_slipweakening_z
    v = resid_z
    execute_on = 'TIMESTEP_BEGIN'
  []
  [restore_x]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_x'
    variable = 'resid_x'
    execute_on = 'TIMESTEP_END'
  []
  [restore_y]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_y'
    variable = 'resid_y'
    execute_on = 'TIMESTEP_END'
  []
  [restore_z]
    type = TagVectorAux
    vector_tag = 'restore_tag'
    v = 'disp_z'
    variable = 'resid_z'
    execute_on = 'TIMESTEP_END'
  []
  ###
  # [get_jump_x_aux]
  #   type = MaterialRealAux
  #   property = jump_x
  #   variable = jump_x_aux
  #   boundary = 'Block100_Block200'
  # []
  # [get_jump_x_rate_aux]
  #   type = FDCompVarRate
  #   variable = jump_x_rate_aux
  #   coupled = jump_x
  #   execute_on = 'TIMESTEP_END'
  #   boundary = 'Block100_Block200'
  # []
  # [get_traction_x_aux]
  #   type = MaterialRealAux
  #   property = traction_x
  #   variable = traction_x_aux
  #   boundary = 'Block100_Block200'
  # []
  ###get CDB model properties
  [get_alpha_damagedvar]
      type = MaterialRealAux
      variable = alpha_damagedvar_aux
      property = alpha_damagedvar
      execute_on = 'TIMESTEP_END'
  []
  [get_B]
      type = MaterialRealAux
      variable = B_aux
      property = B
      execute_on = 'TIMESTEP_END'
  []
  [get_xi]
      type = MaterialRealAux
      variable = xi_aux
      property = xi
      execute_on = 'TIMESTEP_END'
  []
  ###
[]

[Kernels]
  [./inertia_x]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_x
  []
  [./inertia_y]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_y
  []
  [./inertia_z]
    type = InertialForce
    use_displaced_mesh = false
    variable = disp_z
  []
  [./Reactionx]
    type = StiffPropDamping
    variable = 'disp_x'
    component = '0'
  []
  [./Reactiony]
    type = StiffPropDamping
    variable = 'disp_y'
    component = '1'
  []
  [./Reactionz]
    type = StiffPropDamping
    variable = 'disp_z'
    component = '2'
  []
[]

[Materials]
  #damage breakage model
  [stress_medium]
      type = ComputeDamageBreakageStress3DSlipWeakening
      output_properties = 'B alpha_damagedvar xi Cd_mat deviatoric_strain_rate'
      use_strain_rate_dependent_Cd = ${use_strain_rate_dependent_Cd}
      m_exponent = ${m_exponent}
      strain_rate_hat = ${strain_rate_hat}
      cd_hat = ${cd_hat}
      outputs = exodus
  []
  [initial_damage_surround]
      type = InitialDamageCycleSim3DPlane
      sigma = ${sigma}
      peak_val = ${peak_val}
      len_of_fault_strike = ${len_of_fault_strike}
      len_of_fault_dip = ${len_of_fault_dip}
      nucl_center = ${fault_center}
      output_properties = 'initial_damage'      
      outputs = exodus
  []
  [dummy_material]
      type = GenericConstantMaterial
      prop_names = 'initial_breakage damage_perturbation'
      prop_values = '0 0'
  []
  [density]
      type = GenericConstantMaterial
      prop_names = density
      prop_values = ${density}
  []
  [./czm_mat]
      type = SlipWeakeningFrictionczm3dCDBM
      disp_slipweakening_x     = disp_slipweakening_x
      disp_slipweakening_y     = disp_slipweakening_y
      disp_slipweakening_z     = disp_slipweakening_z
      vel_slipweakening_x      = vel_slipweakening_x
      vel_slipweakening_y      = vel_slipweakening_y
      vel_slipweakening_z      = vel_slipweakening_z
      reaction_slipweakening_x = resid_slipweakening_x
      reaction_slipweakening_y = resid_slipweakening_y
      reaction_slipweakening_z = resid_slipweakening_z
      mu_s = ${mu_s}
      mu_d = ${mu_d}
      Dc = ${Dc}
      len = ${elem_size}
      use_cohesion = ${use_cohesion}
      cohesion_function = 'func_cohesion'
      boundary = 'Block100_Block200'
  [../]
  [./static_initial_strain_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_strain_tensor
      tensor_functions = 'func_initial_strain_xx   func_initial_strain_xy      func_initial_strain_xz 
                          func_initial_strain_xy   func_initial_strain_yy      func_initial_strain_yz
                          func_initial_strain_xz   func_initial_strain_yz      func_initial_strain_zz'
  [../]
  [./static_initial_stress_tensor]
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
        tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
  [../]
[]

[Functions]
  #cohesion 
  [./func_cohesion]
    type = ParsedFunction
    expression = ${cohesion_expression}
  []
  ###
  #the initial shear stress needs additional nucleation parameters
  [./func_initial_stress_xy]
      type = InitialShearStressCDBM
      peak_value = ${peak_shear_value}
      nucl_center_x = ${nucl_center_x}
      nucl_center_z = ${nucl_center_z}
      nucl_size = ${nucl_size}
      elem_size = ${elem_size}
      solution = init_sol_components
      from_variable = 'stress_01'
  []
  ###
  [./func_initial_strain_xx]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_00'
  []
  [./func_initial_strain_xy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_01'
  []
  [./func_initial_strain_xz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_02'
  []
  [./func_initial_strain_yy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_11'
  []
  [./func_initial_strain_yz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_12'
  []
  [./func_initial_strain_zz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'elastic_strain_22'
  []
  ###
  [./func_initial_stress_xx]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_00'
  []
  # [./func_initial_stress_xy]
  #   type = SolutionFunction
  #   solution = init_sol_components
  #   from_variable = 'stress_01'
  # []
  [./func_initial_stress_xz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_02'
  []
  [./func_initial_stress_yy]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_11'
  []
  [./func_initial_stress_yz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_12'
  []
  [./func_initial_stress_zz]
    type = SolutionFunction
    solution = init_sol_components
    from_variable = 'stress_22'
  []
[]

[UserObjects]
  [recompute_residual_tag]
      type = ResidualEvaluationUserObject
      vector_tag = 'restore_tag'
      force_preaux = true
      execute_on = 'TIMESTEP_END'
  []
  [./init_sol_components]
    type = SolutionUserObject
    mesh = '../static_solve/static_solve_out.e'
    system_variables = 'elastic_strain_00 elastic_strain_01 elastic_strain_02
                        elastic_strain_11 elastic_strain_12 elastic_strain_22
                        stress_00 stress_01 stress_02 stress_11 stress_12 stress_22'
    timestep = LATEST
    force_preaux = true
  [../]  
[]

[Executioner]
  type = Transient
  dt = ${dt}
  end_time = ${end_time}
  # num_steps = ${num_steps}
  [TimeIntegrator]
    type = CentralDifference
    solve_type = lumped
  []
[]

[Outputs]
  exodus = true
  show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z alpha_damagedvar_aux B_aux xi_aux stress_xx stress_yy stress_xy Cd_mat deviatoric_strain_rate'
  time_step_interval = ${exodus_time_step_interval}
  [csv]
    type = CSV
    execute_on = 'timestep_end'
    time_step_interval = ${csv_time_step_interval}
  []
  [out]
    type = Checkpoint
    time_step_interval = ${checkpoint_time_step_interval}
    num_files = ${checkpoint_num_files}
  []
[]    

[VectorPostprocessors]
  [main_fault]
    type = SideValueSampler
    variable = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z ' 
    boundary = 'Block100_Block200'
    sort_by = x
  []
[]