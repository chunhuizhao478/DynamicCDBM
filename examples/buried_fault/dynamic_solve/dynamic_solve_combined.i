zmin = -10000.0
CBH_constant = 10000.0
end_time = 50.0
normal_traction_y = 127500000.0
shear_modulus_o = 30000000000.0
m2 = 1.0
normal_traction_x = 135000000.0
ymin = -10000.0
xi_d = -0.9
nucl_thickness = 200.0
shear_traction = 55000000.0
zmax = 10000.0
C_2 = 0.05
C_1 = 300.0
C_g = 1e-10
xi_0 = -0.8
dt = 0.0001
normal_traction_z = 120000000.0
Cd_constant = 10000.0
chi = 0.8
xmax = 10000.0
xmin = -10000.0
rho = 2700.0
e_damage = 0.3
CdCb_multiplier = 1000.0
m1 = 10.0
beta_width = 0.03
lambda_o = 30000000000.0
nucl_center_y = -5000.0
time_step_interval = 1000.0
nucl_center_z = 0.0
nucl_center_x = 0.0
duration = 0.1
nucl_distance = 400.0

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/buried_fault.msh'
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
        new_boundary = 'left right bottom top back front'
    []
    [./extranodeset1]
        type = ExtraNodesetGenerator
        coord = '${xmin}  ${ymin}  ${zmin};
                 ${xmax}  ${ymin}  ${zmin};
                 ${xmin}  ${ymin}  ${zmax};
                 ${xmax}  ${ymin}  ${zmax}'
        new_boundary = corner_ptr
        input = sidesets
    []  
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
    beta_width = ${beta_width} #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = ${C_g}
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = ${m1}
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = ${m2}
    
    # energy ratio
    chi = ${chi}

    #
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
    [alpha_grad_x]
    []
    [alpha_grad_y]
    []    
    [alpha_grad_z]
    []
    [vel_x]
    []  
    [vel_y]
    []
    [vel_z]
    []
    [initial_damage_aux]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[AuxKernels]
    [Vel_x]
        type = CompVarRate
        variable = vel_x
        coupled = disp_x
        execute_on = 'TIMESTEP_END'
    []
    [Vel_y]
        type = CompVarRate
        variable = vel_y
        coupled = disp_y
        execute_on = 'TIMESTEP_END'
    []
    [Vel_z]
      type = CompVarRate
      variable = vel_z
      coupled = disp_z
      execute_on = 'TIMESTEP_END'
    []
    [initial_damage]
        type = SolutionAux
        variable = initial_damage_aux
        solution = init_sol_components
        from_variable = initial_damage
    []
[]

[Kernels]
    [dispkernel_x]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = StressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
    []
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
[]

[AuxKernels]
[]

[Materials]
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        # outputs = exodus
    [] 
    [density]
        type = GenericConstantMaterial
        prop_names = density
        prop_values = ${rho}
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3D
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi'
        outputs = exodus
    [] 
    [initial_damage]
        type = ParsedMaterial
        property_name = initial_damage
        coupled_variables = initial_damage_aux
        expression = 'initial_damage_aux'
        outputs = exodus
    []
    [damage_perturb]
        type = DamagePerturbationSquare
        nucl_center = '${nucl_center_x} ${nucl_center_y} ${nucl_center_z}'
        e_damage = ${e_damage}
        thickness = ${nucl_thickness}
        length = ${nucl_distance}
        duration = ${duration}
        outputs = exodus
    []
[]  

[Functions]
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_combined_out.e'
      system_variables = 'disp_x disp_y disp_z initial_damage'
      timestep = LATEST
      force_preaux = true
    [../]
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]
  
[Executioner]
    type = Transient
    dt = ${dt}
    end_time = ${end_time}
    # num_steps = 10
    [TimeIntegrator]
        type = CentralDifference
        solve_type = lumped
        use_constant_mass = true
    []
[]

[Outputs] 
    exodus = true
    time_step_interval = ${time_step_interval}
    # [sample_snapshots]
    #     type = Exodus
    #     time_step_interval = 2000
    # []
    # [./checkpoint]
    #     type = Checkpoint
    #     wall_time_interval = 4000 # interval length in seconds
    # [../]    
[]

#We assume the simulation is loaded with compressive pressure and shear stress
[BCs]
    [pressure_right]
        type = ADPressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        factor = ${normal_traction_x}
    []
    [pressure_left]
        type = ADPressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        factor = ${normal_traction_x}
    []
    [pressure_front]
        type = ADPressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        factor = ${normal_traction_z}
    []
    [pressure_back]
        type = ADPressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        factor = ${normal_traction_z}  
    []
    [pressure_top]
        type = ADPressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = top
        factor = ${normal_traction_y}      
    []
    [pressure_bottom]
        type = ADPressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = bottom
        factor = ${normal_traction_y}            
    []
    #
    [pressure_shear_front]
        type = ADNeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        value = ${shear_traction}
    []
    [pressure_shear_back]
        type = ADNeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        value = -${shear_traction}  
    []
    [pressure_shear_left]
        type = ADNeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        value = -${shear_traction}
    []
    [pressure_shear_right]
        type = ADNeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        value = ${shear_traction}    
    []
[]

[ICs]
    [disp_x_ic]
      type = SolutionIC
      variable = disp_x
      solution_uo = init_sol_components
      from_variable = disp_x
    []
    [disp_y_ic]
      type = SolutionIC
      variable = disp_y
      solution_uo = init_sol_components
      from_variable = disp_y
    []
    [disp_z_ic]
      type = SolutionIC
      variable = disp_z
      solution_uo = init_sol_components
      from_variable = disp_z
    []
[]