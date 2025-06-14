#parameters

##material properties##
density = 2670 #density
lambda_o = 3.204e10 #first lame constant
shear_modulus_o = 3.204e10 #second lame constant
shear_modulus_d = 3.204e10 #second lame constant ## damaged shear modulus as input
Cs = '${fparse shear_modulus_d / density }' #shear wave speed
Cp = '${fparse (lambda_o + 2 * shear_modulus_d) / density }' #pressure wave speed
##-------------------------##

##Slip weakening parameters##
Dc = 0.4 #characteristic length (m)
q = 0.4 #damping ratio
mu_s = 0.677 #static friction coefficient
mu_d = 0.525 #dynamic friction coefficient
##-------------------------##

##CDB model parameters##
xi_0 = -0.6 #strain invariants ratio: onset of damage evolution
xi_d = -0.7 #strain invariants ratio: onset of breakage healing
Cd_constant = 0 #coefficient gives positive damage evolution
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

##initial stress parameters##
peak_shear_value = 70e6 #initial shear stress perturbation peak value
domain_shear_value = 70e6 #initial shear stress perturbation domain value
nucl_center_x = 0 #nucleation center x coordinate
nucl_center_y = -2500 #nucleation center y coordinate
nucl_size = 3000 #nucleation size
initial_normal_stress_value = -120e6 #initial normal stress value
initial_dip_stress_value  = 0 #initial dip stress value
initial_strike_stress_value = -120e6 #initial strike stress value
##-------------------------##

##model parameters##
elem_size = 200 #element size for initial shear stress perturbation
dt = 0.0025 #time step size
end_time = 4 #end time for simulation
exodus_time_step_interval = 1 #time step interval for output
checkpoint_time_step_interval = 80 #time step interval for checkpoint output
checkpoint_num_files = 2 #number of files for checkpoint output
##-------------------------##

# ##initial damage parameters##
# initial_damage_domain_size = '${fparse 2 * elem_size }' #initial damage domain size
# initial_damage_value = 0.7 #initial damage value within the fault zone
# ##-------------------------##

##coordinates ##
#x: strike direction
#y: dip direction (free surface is at y = 0, going down is negative) 
#z: normal direction
##-------------------------##
[Mesh]
    [msh]
      type = GeneratedMeshGenerator
      dim = 3
      xmin = -4000
      xmax = 4000
      ymin = -5000
      ymax = 0
      zmin = -2000
      zmax = 2000
      nx = 40
      ny = 25
      nz = 20
      subdomain_ids = 1
    []
    # [msh]
    #   type = GeneratedMeshGenerator
    #   dim = 3
    #   xmin = -15000
    #   xmax = 15000
    #   ymin = -20000
    #   ymax = 0
    #   zmin = -10000
    #   zmax = 10000
    #   nx = 150
    #   ny = 100
    #   nz = 100
    #   subdomain_ids = 1
    # []
    [./new_block_1]
      type = ParsedSubdomainMeshGenerator
      input = msh
      combinatorial_geometry = 'x >= -13000 & x <= 13000 & y > -15000 & z < 0'
      block_id = 2
    []
    [./new_block_2]
        type = ParsedSubdomainMeshGenerator
        input = new_block_1
        combinatorial_geometry = 'x > -13000 & x < 13000 & y > -15000 & z > 0'
        block_id = 3
    []      
    [./split_1]
        type = BreakMeshByBlockGenerator
        input = new_block_2
        split_interface = true
        block_pairs = '2 3'
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
[]
    
[GlobalParams]
    ##------------slip weakening------------##
    displacements = 'disp_x disp_y disp_z'
    
    #damping ratio
    q = ${q}
    
    #characteristic length (m) #this gives ~200m for resolve L_f
    Dc = ${Dc}

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
    #restoration force for damping (tag after solve)
    [./resid_damp_x]
        order = FIRST
        family = LAGRANGE
    [../]
    [./resid_damp_y]
        order = FIRST
        family = LAGRANGE
    [../] 
    [./resid_damp_z]
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
    [./accel_slipweakening_x]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_slipweakening_y]
        order = FIRST
        family = LAGRANGE
    []
    [./accel_slipweakening_z]
        order = FIRST
        family = LAGRANGE
    []
    [./mu_s]
        order = CONSTANT
        family = MONOMIAL
    []
    [./mu_d]
        order = CONSTANT
        family = MONOMIAL
    []
    [./ini_shear_stress]
        order = FIRST
        family = LAGRANGE
    []
    #
    [./tangent_jump]
        order = CONSTANT
        family = MONOMIAL
    []
    [./tangent_jump_rate]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [./elem_length]
        order = CONSTANT
        family = MONOMIAL
    [../]
    #
    [./jump_x]
        order = CONSTANT
        family = MONOMIAL        
    []
    [./jump_y]
        order = CONSTANT
        family = MONOMIAL        
    []
    [./jump_z]
        order = CONSTANT
        family = MONOMIAL         
    []
    #
    [jump_rate_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [jump_rate_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [jump_rate_z]
        order = CONSTANT
        family = MONOMIAL
    []
    #
    [traction_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [traction_z]
        order = CONSTANT
        family = MONOMIAL
    [] 
    #grad_alpha
    [./alpha_grad_x]
        order = CONSTANT
        family = MONOMIAL
    []
    [./alpha_grad_y]
        order = CONSTANT
        family = MONOMIAL
    []
    [./alpha_grad_z]
        order = CONSTANT
        family = MONOMIAL
    [] 
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

[Physics]
    [SolidMechanics]
        [QuasiStatic]
            [./all]
                strain = SMALL
                add_variables = true
                extra_vector_tags = 'restore_tag'
                eigenstrain_names = static_initial_strain_tensor
            [../]
        [../]
    [../]
[]

[Problem]
    extra_tag_vectors = 'restore_tag restore_dampx_tag restore_dampy_tag restore_dampz_tag'
[]

[AuxKernels]
    [Displacment_x]
        type = ProjectionAux
        variable = disp_slipweakening_x
        v = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Displacement_y]
        type = ProjectionAux
        variable = disp_slipweakening_y
        v = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [Displacement_z]
        type = ProjectionAux
        variable = disp_slipweakening_z
        v = disp_z
        execute_on = 'TIMESTEP_END'
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
    #
    [XJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [YJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [ZJump]
        type = MaterialRealVectorValueAux
        property = displacement_jump_global
        variable = jump_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    #
    [XJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [YJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [ZJumpRate]
        type = MaterialRealVectorValueAux
        property = displacement_jump_rate_global
        variable = jump_rate_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []    
    #
    [TractionX]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_x
        component = 0
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'        
    []
    [TractionY]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_y
        component = 1
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []
    [TractionZ]
        type = MaterialRealVectorValueAux
        property = traction_on_interface
        variable = traction_z
        component = 2
        execute_on = 'TIMESTEP_END'
        boundary = 'Block2_Block3'
    []        
    #
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
    #damping
    [restore_dampx]
        type = TagVectorAux
        vector_tag = 'restore_dampx_tag'
        v = 'disp_x'
        variable = 'resid_damp_x'
        execute_on = 'TIMESTEP_END'
    []
    [restore_dampy]
        type = TagVectorAux
        vector_tag = 'restore_dampy_tag'
        v = 'disp_y'
        variable = 'resid_damp_y'
        execute_on = 'TIMESTEP_END'
    []
    [restore_dampz]
        type = TagVectorAux
        vector_tag = 'restore_dampz_tag'
        v = 'disp_z'
        variable = 'resid_damp_z'
        execute_on = 'TIMESTEP_END'
    []
    [StaticFricCoeff]
        type = FunctionAux
        variable = mu_s
        function = func_static_friction_coeff_mus
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [DynamicFricCoeff]
        type = FunctionAux
        variable = mu_d
        function = func_dynamic_friction_coeff_mud
        execute_on = 'INITIAL TIMESTEP_BEGIN'
    []
    [elem_length]
        type = ConstantAux
        variable = elem_length
        value = 200
    []
    #get CDB model properties
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
        extra_vector_tags = restore_dampx_tag
    []
    [./Reactiony]
        type = StiffPropDamping
        variable = 'disp_y'
        component = '1'
        extra_vector_tags = restore_dampy_tag
    []
    [./Reactionz]
        type = StiffPropDamping
        variable = 'disp_z'
        component = '2'
        extra_vector_tags = restore_dampz_tag
    []
[]

[InterfaceKernels]
    [czm_interface_kernel_x]
        type = FarmsCZM
        variable = disp_x
        neighbor_var = disp_x
        boundary = 'Block2_Block3'
    []
    [czm_interface_kernel_y]
        type = FarmsCZM
        variable = disp_y
        neighbor_var = disp_y
        boundary = 'Block2_Block3'
    []
    [czm_interface_kernel_z]
        type = FarmsCZM
        variable = disp_z
        neighbor_var = disp_z
        boundary = 'Block2_Block3'
    []
[]

[Materials]
    #damage breakage model
    [stress_medium]
        type = ComputeDamageBreakageStress3DSlipWeakening
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi eps_p eps_e I1 I2 sts_total'
        outputs = exodus
    []
    #dummy material properties needed for stress_medium
    [initial_damage]
        type = GenericConstantMaterial
        prop_names = initial_damage
        prop_values = 0
    []
    # [initial_damage]
    #     type = NearFaultInitialDamage
    #     distance_from_fault = ${initial_damage_domain_size}
    #     initial_damage_value = ${initial_damage_value}
    #     output_properties = initial_damage
    #     outputs = exodus
    # []
    [initial_breakage]
        type = GenericConstantMaterial
        prop_names = initial_breakage
        prop_values = 0
    []
    [damage_perturbation]
        type = GenericConstantMaterial
        prop_names = damage_perturbation
        prop_values = 0
    []
    ##
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = ${density}
    []
    [elasticity]
        type = ComputeIsotropicElasticityTensor
        shear_modulus = ${shear_modulus_d}
        lambda = ${lambda_o}
    []
    [./czm_mat]
        type = FarmsSlipWeakeningCZM
        disp_slipweakening_x     = disp_slipweakening_x
        disp_slipweakening_y     = disp_slipweakening_y
        disp_slipweakening_z     = disp_slipweakening_z
        vel_slipweakening_x      = vel_slipweakening_x
        vel_slipweakening_y      = vel_slipweakening_y
        vel_slipweakening_z      = vel_slipweakening_z
        reaction_slipweakening_x = resid_x
        reaction_slipweakening_y = resid_y
        reaction_slipweakening_z = resid_z
        reaction_damp_x = resid_damp_x
        reaction_damp_y = resid_damp_y
        reaction_damp_z = resid_damp_z
        elem_length = elem_length
        mu_d = mu_d
        mu_s = mu_s
        boundary = 'Block2_Block3'
    [../]
    [./strain_from_initial_stress]
        type = ComputeEigenstrainFromInitialStress
        initial_stress = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
        eigenstrain_name = static_initial_strain_tensor
    [../] 
    # [./strain_from_initial_stress]
    #     type = ComputeDamageBreakageEigenstrainFromInitialStress
    #     initial_stress = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
    #                       func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
    #                       func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
    #     eigenstrain_name = static_initial_strain_tensor
    #     lambda_o = ${lambda_o}
    #     shear_modulus_o = ${shear_modulus_o}
    #     xi_o = ${xi_0}
    # [../] 
    [./static_initial_stress_tensor_slipweakening]
        type = GenericFunctionRankTwoTensor
        tensor_name = static_initial_stress_tensor_slipweakening
        tensor_functions = 'func_initial_stress_xx   func_initial_stress_xy      func_initial_stress_xz 
                            func_initial_stress_xy   func_initial_stress_yy      func_initial_stress_yz
                            func_initial_stress_xz   func_initial_stress_yz      func_initial_stress_zz'
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
    [func_static_friction_coeff_mus]
        type = ConstantFunction
        value = ${mu_s}
    []
    #mud constant value
    [func_dynamic_friction_coeff_mud]
        type = ConstantFunction
        value = ${mu_d}
    []
    #Note:restrict stress variation along the fault only
    #this function is used in czm only
    [./func_initial_stress_xx]
        type = ConstantFunction
        value = ${initial_strike_stress_value}
    []
    [./func_initial_stress_xy]
        type = ConstantFunction
        value = 0
    []
    [./func_initial_stress_xz]
        type = InitialShearStressCDBM
        peak_value = ${peak_shear_value}
        domain_value = ${domain_shear_value}
        nucl_center_x = ${nucl_center_x}
        nucl_center_y = ${nucl_center_y}
        nucl_size = ${nucl_size}
        elem_size = ${elem_size}
    []
    [./func_initial_stress_yy]
        type = ConstantFunction
        value = ${initial_dip_stress_value}
    []
    [./func_initial_stress_yz]
        type = ConstantFunction
        value = 0
    []
    [./func_initial_stress_zz]
        type = ConstantFunction
        value = ${initial_normal_stress_value}
    []
[]

[UserObjects]
    [recompute_residual_tag]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    #damping
    [recompute_residual_tag_dampx]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampx_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampy]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampy_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
    [recompute_residual_tag_dampz]
        type = ResidualEvaluationUserObject
        vector_tag = 'restore_dampz_tag'
        force_preaux = true
        execute_on = 'TIMESTEP_END'
    []
[]

[Executioner]
    type = Transient
    dt = ${dt}
    end_time = ${end_time}
    # num_steps = 1
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []
[]

[Outputs]
    exodus = true
    time_step_interval = ${exodus_time_step_interval}
    show = 'vel_slipweakening_x vel_slipweakening_y vel_slipweakening_z disp_slipweakening_x disp_slipweakening_y disp_slipweakening_z B_aux alpha_damagedvar_aux xi_aux'
    [checkpoint]
        type = Checkpoint
        time_step_interval = ${checkpoint_time_step_interval}
        num_files = ${checkpoint_num_files}
    []
[]

[BCs]
    ##non-reflecting bc
    #
    [./dashpot_bottom_x]
        type = NonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = bottom
    []
    [./dashpot_bottom_y]
        type = NonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = bottom
    []
    [./dashpot_bottom_z]
        type = NonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = bottom
    []
    #
    [./dashpot_left_x]
        type = NonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = left
    []
    [./dashpot_left_y]
        type = NonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = left
    []
    [./dashpot_left_z]
        type = NonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = left
    []
    #
    [./dashpot_right_x]
        type = NonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = right
    []
    [./dashpot_right_y]
        type = NonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = right
    []
    [./dashpot_right_z]
        type = NonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = right
    []
    #
    [./dashpot_front_x]
        type = NonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = front
    []
    [./dashpot_front_y]
        type = NonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = front
    []
    [./dashpot_front_z]
        type = NonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = front
    []
    #
    [./dashpot_back_x]
        type = NonReflectDashpotBC3d
        component = 0
        variable = disp_x
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = back
    []
    [./dashpot_back_y]
        type = NonReflectDashpotBC3d
        component = 1
        variable = disp_y
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = back
    []
    [./dashpot_back_z]
        type = NonReflectDashpotBC3d
        component = 2
        variable = disp_z
        disp_x = disp_x
        disp_y = disp_y
        disp_z = disp_z
        p_wave_speed = ${Cp}
        shear_wave_speed = ${Cs}
        boundary = back
    []
[]