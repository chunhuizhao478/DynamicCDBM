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
    [initial_damage_aux]
        order = FIRST
        family = MONOMIAL
    []
    [correlated_randalpha_o]
        order = FIRST
        family = LAGRANGE
    []
    [initial_cd_aux]
        order = FIRST
        family = MONOMIAL
    []
[]

[AuxKernels]
    [get_initial_damage]
        type = ADMaterialRealAux
        variable = initial_damage_aux
        property = initial_damage
    []
[]

[Kernels]
    [dispkernel_x]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_x
        component = 0
    []
    [dispkernel_y]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_y
        component = 1
    []
    [dispkernel_z]
        type = ADStressDivergenceTensors
        displacements = 'disp_x disp_y disp_z'
        variable = disp_z
        component = 2
    []
[]

[Materials]
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        outputs = exodus
    [] 
    [stress]
        type = ADComputeDamageStressStaticDistribution
        lambda_o = 30e9
        shear_modulus_o = 30e9
        xi_o = -0.8
        chi = 0.7
        xi_d = -0.9
        outputs = exodus
        block = '1 3'
    [] 
    [stress_elastic]
        type = ADComputeLinearElasticStress
        outputs = exodus
        block = '2'
    []
    [elasticity_tensor]
        type = ADComputeIsotropicElasticityTensor
        lambda = 30e9
        shear_modulus = 30e9
    []
    [getxi]
        type = ADComputeXi
        outputs = exodus
    []
    [initial_damage_surround]
        type = ADInitialDamageCycleSim3DPlane
        sigma = 5e2
        peak_val = 0.7
        len_of_fault_strike = 8000
        len_of_fault_dip = 3000
        nucl_center = '0 0 -7500'
        output_properties = 'initial_damage'      
        outputs = exodus
    []
    [initial_breakage_surround]
        type = ADInitialBreakageCycleSim3DPlane
        sigma = 5e2
        peak_val = 0.1
        len_of_fault_strike = 8000
        len_of_fault_dip = 3000
        nucl_center = '0 0 -7500'
        output_properties = 'initial_breakage'      
        outputs = exodus
    []
[]  

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]
  
[Executioner]
    type = Steady
    solve_type = 'NEWTON'
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-10
    nl_max_its = 20
    nl_abs_tol = 1e-12
    # this is very robust, use as default
    petsc_options_iname = '-ksp_type -pc_type -ksp_initial_guess_nonzero'
    petsc_options_value = 'gmres     hypre  True'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
[]  

[Outputs]
    exodus = true       
[]

[BCs]
    #Note: use neuamnnBC gives minimum waves than pressureBC
    [static_pressure_top]
        type = ADNeumannBC
        variable = disp_z
        boundary = top
        value = -50e6
        displacements = 'disp_x disp_y disp_z'
    []
    [static_pressure_bottom]
        type = ADNeumannBC
        variable = disp_z
        boundary = bottom
        value = 50e6
        displacements = 'disp_x disp_y disp_z'
    []     
    [static_pressure_left]
        type = ADNeumannBC
        variable = disp_x
        boundary = left
        value = 50e6
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_right]
        type = ADNeumannBC
        variable = disp_x
        boundary = right
        value = -50e6
        displacements = 'disp_x disp_y disp_z'
    [] 
    [static_pressure_front]
        type = ADNeumannBC
        variable = disp_y
        boundary = front
        value = 50e6
        displacements = 'disp_x disp_y disp_z'
    []  
    [static_pressure_back]
        type = ADNeumannBC
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
        type = ADDirichletBC
        variable = disp_x
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_y]
        type = ADDirichletBC
        variable = disp_y
        boundary = corner_ptr
        value = 0
    []
    [./fix_cptr1_z]
        type = ADDirichletBC
        variable = disp_z
        boundary = corner_ptr
        value = 0
    []     
[]