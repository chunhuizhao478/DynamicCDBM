#continuum damage-breakage model dynamics

##########################################################################################################################################
#Mesh section
#FileMeshGenerator: read mesh file
#SideSetsFromNormalsGenerator: generate side sets from normals
#ExtraNodesetGenerator: generate extra nodeset - here we use it to define corner points associated with the bottom boundary
##########################################################################################################################################
[Mesh]
    [./msh]
        type = FileMeshGenerator
        file = '../meshfile/mesh_singleblock_large.msh'
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
        coord = ' -30000 -30000 -30000'
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
[]

[Materials]
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        outputs = exodus
    [] 
    [stress_elastic]
        type = ComputeLinearElasticStress
        outputs = exodus
    []
    [elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        lambda = 32.04e9
        shear_modulus = 32.04e9
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
    nl_rel_tol = 1e-6
    nl_max_its = 100
    nl_abs_tol = 1e-8
    # this is very robust, use as default
    # petsc_options_iname = '-ksp_type -pc_type -ksp_initial_guess_nonzero'
    # petsc_options_value = 'gmres     hypre  True'
    petsc_options_iname = '-pc_type -pc_factor_shift_type'
    petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    # automatic_scaling = true
[]  

[Outputs]
    exodus = true       
[]

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
    #
    [static_pressure_front_shear]
        type = NeumannBC
        variable = disp_x
        boundary = front
        value = -30e6
        displacements = 'disp_x disp_y disp_z'
    []
    [static_pressure_back_shear]
        type = NeumannBC
        variable = disp_x
        boundary = back
        value = 30e6
        displacements = 'disp_x disp_y disp_z'
    []
    [static_pressure_left_shear]
        type = NeumannBC
        variable = disp_y
        boundary = left
        value = -30e6
        displacements = 'disp_x disp_y disp_z'
    []
    [static_pressure_right_shear]
        type = NeumannBC
        variable = disp_y
        boundary = right
        value = 30e6
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