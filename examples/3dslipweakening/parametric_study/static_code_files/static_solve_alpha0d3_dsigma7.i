#parameters

##mesh parameters
###need to check gmsh file for changing the parameters here
bottom_nodes_coord =' -40000 -40000 -200000;
                      40000 -40000 -200000;
                      40000 40000  -200000;
                     -40000 40000  -200000'
##boundary loading parameters
confining_pressure = 120e6 #Pa, confining pressure
shear_traction = 67e6 #Pa, shear traction
##-------------------------##

##material properties##
density = 2670 #density
lambda_o = 3.204e10 #first lame constant
shear_modulus_o = 3.204e10 #second lame constant
# Cs = '${fparse shear_modulus_o / density }' #shear wave speed
# Cp = '${fparse (lambda_o + 2 * shear_modulus_o) / density }' #pressure wave speed
##-------------------------##

##CDB model parameters##
xi_0 = -0.8 #strain invariants ratio: onset of damage evolution
xi_d = -0.9 #strain invariants ratio: onset of breakage healing

###constant Cd
Cd_constant = 0 #coefficient gives positive damage evolution
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
peak_val = 0.3
len_of_fault_strike = 30000
len_of_fault_dip = 15000
fault_center = '0 0 -9500'
##-------------------------##

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../../mesh/v3_100m_freesurface_coarsened.msh'
    []
    [./sidesets]
        input = msh
        type = SideSetsFromNormalsGenerator
        normals = '-1 0 0
                    1 0 0
                    0 -1 0
                    0 1 0
                    0 0 -1
                    0 0 1'
        new_boundary = 'left right front back bottom top'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = ${bottom_nodes_coord}
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'

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
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
  []
  [disp_z]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
[]

[AuxKernels]
[]

[Kernels]
  [SolidMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = ${shear_modulus_o}
    lambda = ${lambda_o}
  [../]
  [strain]
    type = ComputeSmallStrain
  []
  [stress_medium]
    type = ComputeDamageBreakageStress3DStatic
    output_properties = 'stress elastic_strain'
    outputs = exodus
  []
  [dummy_material]
      type = GenericConstantMaterial
      prop_names = 'initial_breakage damage_perturbation density'
      prop_values = '0 0 ${density}'
  []
  [initial_damage_surround]
    type = InitialDamageCycleSim3DPlane
    sigma = ${sigma}
    peak_val = ${peak_val}
    len_of_fault_strike = ${len_of_fault_strike}
    len_of_fault_dip = ${len_of_fault_dip}
    nucl_center = ${fault_center}
  []
  [./static_initial_strain_tensor] #this is used in the ComputeDamageBreakageStress3DSlipWeakening
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                      func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                      func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
    eigenstrain_name = ini_stress
  [../]
  [./static_initial_stress_tensor] #this is used in the ComputeDamageBreakageStress3DSlipWeakening, SlipWeakeningFrictionczm3dCDBM
      type = GenericFunctionRankTwoTensor
      tensor_name = static_initial_stress_tensor
      tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                          func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                          func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
  [../]    
  [./comp_xi]
    type = ComputeXi
  []
[]

[BCs]
  [static_pressure_left]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = left
      function = func_pos_xx_stress
      displacements = 'disp_x disp_y disp_z'
  []  
  [static_pressure_right]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = right
      function = func_neg_xx_stress
      displacements = 'disp_x disp_y disp_z'
  [] 
  #
  [static_pressure_front]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = front
      function = func_pos_yy_stress
      displacements = 'disp_x disp_y disp_z'
  []  
  [static_pressure_back]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = back
      function = func_neg_yy_stress
      displacements = 'disp_x disp_y disp_z'
  [] 
  #
  [static_pressure_front_shear]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = front
      function = func_neg_xy_stress
      displacements = 'disp_x disp_y disp_z'
  []  
  [static_pressure_back_shear]
      type = FunctionNeumannBC
      variable = disp_x
      boundary = back
      function = func_pos_xy_stress
      displacements = 'disp_x disp_y disp_z'
  [] 
  [static_pressure_left_shear]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = left
      function = func_neg_xy_stress
      displacements = 'disp_x disp_y disp_z'
  []  
  [static_pressure_right_shear]
      type = FunctionNeumannBC
      variable = disp_y
      boundary = right
      function = func_pos_xy_stress
      displacements = 'disp_x disp_y disp_z'
  []   
  #
  [fix_node_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'corner_ptr'
    value = 0
  [../]
  [fix_node_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'corner_ptr'
    value = 0
  [../]
  [fix_node_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'corner_ptr'
    value = 0
  [../]
[]

[Functions]
  ###stress field###
  [./func_initial_stress_xx]
    type = ConstantFunction
    value = ${fparse confining_pressure}
  []
  [./func_pos_xx_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xx'
    scale_factor = '1'
  [../]
  [./func_neg_xx_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xx'
    scale_factor = '-1'
  [../]
  ##
  [./func_initial_stress_xy]
    type = ConstantFunction
    value = ${fparse shear_traction}
  []
  [./func_pos_xy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xy'
    scale_factor = '1'
  [../]
  [./func_neg_xy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_xy'
    scale_factor = '-1'
  [../]
  ##
  [./func_initial_stress_xz]
    type = ConstantFunction
    value = 0
  []
  ##
  [./func_initial_stress_yy]
    type = ConstantFunction
    value = ${fparse confining_pressure}
  []
  [./func_pos_yy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_yy'
    scale_factor = '1'
  [../]
  [./func_neg_yy_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_yy'
    scale_factor = '-1'
  [../]
  ##
  [./func_initial_stress_yz]
    type = ConstantFunction
    value = 0
  []
  [./func_initial_stress_zz]
    type = ConstantFunction
    value = 0
  []
  [./func_pos_zz_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_zz'
    scale_factor = '1'
  [../]
  [./func_neg_zz_stress]
    type = CompositeFunction
    functions = 'func_initial_stress_zz'
    scale_factor = '-1'
  [../]
[]

[Preconditioning]
  [./andy]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  solve_type = NEWTON
  type = Steady

  nl_abs_tol = 1E-12
  nl_rel_tol = 1E-10
  l_tol = 1E-7
  l_max_its = 200
  nl_max_its = 400
  line_search  = 'bt'
  # automatic_scaling = true
  verbose = true
  petsc_options_iname = '-ksp_type -pc_type -ksp_initial_guess_nonzero'
  petsc_options_value = 'gmres     hypre  True'
[]

[Outputs]
  exodus = true
  show = 'elastic_strain_00 elastic_strain_01 elastic_strain_02 elastic_strain_11 elastic_strain_12 elastic_strain_22 stress_00 stress_01 stress_02 stress_11 stress_12 stress_22'
[]   