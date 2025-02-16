#continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh.msh'
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
        coord = ' -15000 -15000 -15000;
                   15000 -15000 -15000;
                   15000 15000  -15000;
                  -15000 15000  -15000'
        new_boundary = corner_ptr
        input = sidesets
    []
    displacements = 'disp_x disp_y disp_z'
[]

[GlobalParams]

    displacements = 'disp_x disp_y disp_z'
    
    ##----continuum damage breakage model----##
    #initial lambda value (first lame constant) [Pa]
    lambda_o = 30e9
        
    #initial shear modulus value (second lame constant) [Pa]
    shear_modulus_o = 30e9
    
    #<strain invariants ratio: onset of damage evolution>: relate to internal friction angle, refer to "note_mar25"
    xi_0 = -0.8
    
    #<strain invariants ratio: onset of breakage healing>: tunable param, see ggw183.pdf
    xi_d = -0.9
    
    #<strain invariants ratio: maximum allowable value>: set boundary
    #Xu_etal_P15-2D
    #may need a bit space, use 1.5 as boundary
    xi_max = 1.8
    
    #<strain invariants ratio: minimum allowable value>: set boundary
    #Xu_etal_P15-2D
    xi_min = -1.8

    #if option 2, use Cd_constant
    Cd_constant = 0

    #<coefficient gives positive breakage evolution >: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    #The multiplier between Cd and Cb: Cb = CdCb_multiplier * Cd
    CdCb_multiplier = 1000

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 1e4

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 300

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_2 = 0.05

    #<coefficient gives width of transitional region>: see P(alpha), refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    beta_width = 0.03 #1e-3
    
    #<material parameter: compliance or fluidity of the fine grain granular material>: refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    C_g = 1e-10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Table 1
    m1 = 10
    
    #<coefficient of power law indexes>: see flow rule (power law rheology): refer to "Lyak_BZ_JMPS14_splitstrain" Equation 18
    m2 = 1
    
    # energy ratio
    chi = 0.7

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
        prop_names = 'density nonADdensity'
        prop_values = '2700 2700'
    []
    [stress_medium]
        type = ComputeDamageBreakageStress3D
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi eps_p'
        block = '1 3'
        # outputs = exodus
    [] 
    [stress_elastic]
        type = ComputeLinearElasticStress
        block = '2'
    []
    [elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 30e9
        shear_modulus = 30e9
    []
    [initial_damage_surround]
        type = InitialDamageCycleSim3DPlane
        sigma = 5e2
        peak_val = 0.7
        len_of_fault_strike = 8000
        len_of_fault_dip = 3000
        nucl_center = '0 0 -7500'
        output_properties = 'initial_damage'      
        # outputs = exodus
    []
    [initial_breakage_surround]
        type = InitialBreakageCycleSim3DPlane
        sigma = 5e2
        peak_val = 0.1
        len_of_fault_strike = 8000
        len_of_fault_dip = 3000
        nucl_center = '0 0 -7500'
        output_properties = 'initial_breakage'      
        # outputs = exodus
    []
    [dummy_damage_perturb]
        type = GenericConstantMaterial
        prop_names = 'damage_perturb'
        prop_values = '0'
        block = '1 2 3'
    []
    [dummy_material]
        type = GenericConstantMaterial
        prop_names = 'damage_perturbation initial_shear_stress shear_stress_perturbation'
        prop_values = '0 0 0'
    []
[]  

[Functions]
[]

[UserObjects]
    [./init_sol_components]
      type = SolutionUserObject
      mesh = '../static_solve/static_solve_out.e'
      system_variables = 'disp_x disp_y disp_z initial_damage'
      timestep = LATEST
      force_preaux = true
    [../]
[]

[Postprocessors]
    [./maxvelx]
        type = NodalExtremeValue
        variable = vel_x
    [../]
    [./maxvely]
        type = NodalExtremeValue
        variable = vel_y
    [../]
[../]
  
[Executioner]
    type = Transient
    dt = 1e-4
    end_time = 10.0
    # num_steps = 10
    [TimeIntegrator]
        type = CentralDifference
        solve_type = consistent
        # use_constant_mass = true
    []
[]

[Outputs] 
    #save the solution to a exodus file every 0.1 seconds
    exodus = false
    # time_step_interval = 1
    [./csv]
        type = CSV
        time_step_interval = 1
        show = 'maxvelx maxvely'
    [../]
[]

#We assume the simulation is loaded with compressive pressure and shear stress
[BCs]
    #Note: use neuamnnBC gives minimum waves than pressureBC
    [static_pressure_top]
        type = NeumannBC
        variable = disp_z
        boundary = top
        value = -50e6
        displacements = 'disp_x disp_y disp_z'
    []
    [static_pressure_bottom]
        type = NeumannBC
        variable = disp_z
        boundary = bottom
        value = 50e6
        displacements = 'disp_x disp_y disp_z'
    []     
    [static_pressure_left]
        type = NeumannBC
        variable = disp_x
        boundary = left
        value = 50e6
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right]
        type = NeumannBC
        variable = disp_x
        boundary = right
        value = -50e6
        displacements = 'disp_x disp_y disp_z'
    [] 
    [static_pressure_front]
        type = NeumannBC
        variable = disp_y
        boundary = front
        value = 50e6
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back]
        type = NeumannBC
        variable = disp_y
        boundary = back
        value = -50e6
        displacements = 'disp_x disp_y disp_z'
    []
    [static_pressure_front_shear]
        type = ADNeumannBC
        variable = disp_x
        boundary = front
        value = -20e6
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back_shear]
        type = ADNeumannBC
        variable = disp_x
        boundary = back
        value = 20e6
        displacements = 'disp_x disp_y disp_z'
    []      
    # fix ptr
    [./fix_cptr1_x]
        type = DirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_y]
        type = DirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_z]
        type = DirichletBC
        variable = disp_z
        boundary = corner_ptr
        value = 0
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