#explicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        # file =  '../meshfile/cdbm_largeplane_coarse.msh'
        file = '../meshfile/cdbm_largeplane_buried_coarse_newmethod.msh'
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
    [./bottomline]
        type = ExtraNodesetGenerator
        input = sidesets
        new_boundary = 'bottomline'
        coord = '-11000 -17000 0; -10000 -17000 0; -9000 -17000 0; 
                 -8000 -17000 0; -7000 -17000 0; -6000 -17000 0;
                 -5000 -17000 0; -4000 -17000 0; -3000 -17000 0;
                 -2000 -17000 0; -1000 -17000 0; 0 -17000 0;
                 1000 -17000 0; 2000 -17000 0; 3000 -17000 0;
                 4000 -17000 0; 5000 -17000 0; 6000 -17000 0;
                 7000 -17000 0; 8000 -17000 0; 9000 -17000 0;
                 10000 -17000 0; 11000 -17000 0'
    []
    # [./extranodeset1]
    #     type = ExtraNodesetGenerator
    #     coord = '-11000 -17000 2250'
    #     new_boundary = corner_ptr
    #     input = sidesets
    # [] 
    # [./extranodeset2]
    #     type = ExtraNodesetGenerator
    #     coord = '11000 -17000 2250'
    #     new_boundary = corner_ptr2
    #     input = extranodeset1
    # [] 
    # [./extranodeset3]
    #     type = ExtraNodesetGenerator
    #     coord = '11000 0 2250'
    #     new_boundary = corner_ptr3
    #     input = extranodeset2
    # []    
    # allow_renumbering = false
[]

[GlobalParams]
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

[AuxKernels]
[]

[Materials]
    [strain]
        type = ADComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        outputs = exodus
    [] 
    [stress_nucleation]
        type = ADComputeDamageStressStaticDistribution
        lambda_o = 30e9
        shear_modulus_o = 30e9
        xi_o = -0.8
        chi = 0.7
        xi_d = -0.9
        outputs = exodus
    []
    # [stress]
    #     type = ADComputeLinearElasticStress
    # []
    # [elasticity]
    #     type = ADComputeIsotropicElasticityTensor
    #     lambda = 30e9
    #     shear_modulus = 30e9
    # []
    [getxi]
        type = ADComputeXi
        outputs = exodus
    []
    #parameters definition:
    #nucl_center: center of the nucleation region
    #fault_plane: the whole fault plane (xmin xmax ymin ymax zmin zmax)
    #nucl_distance: the distance along strike and dip direction of the nucleation region
    #nucl_thickness: the thickness of the nucleation region along normal direction
    #nucl_damage: the damage value within the nucleation region (say 0.8)
    #e_damage: the peak damage value (say 0.7) and the peak value for exponential decay along the normal direction
    #e_sigma: the standard deviation for the exponential decay along the normal direction
    [initialdamage]
        type = ADInitialDamageBenchmark
        nucl_center = '-7000 -7500 0'
        fault_plane = '-11000 11000 -15000 0 -500 500'
        nucl_distance = 2500
        nucl_thickness = 200
        nucl_damage = 0.7
        e_damage = 0.7
        e_sigma = 1e3
        outputs = exodus
    [] 
    [initialbreakage]
        type = ADInitialBreakageBenchmark
        nucl_center = '-7000 -7500 0'
        fault_plane = '-11000 11000 -15000 0 -500 500'
        nucl_distance = 2500
        nucl_thickness = 200
        nucl_breakage = 0.1
        e_breakage = 0.1
        e_sigma = 1e3
        outputs = exodus
    []
    #block 1: nucleation patch
    #block 3: fault plane
    #block 2: the rest of the domain
    # [initialdamage_block1]
    #     type = ADGenericConstantMaterial
    #     prop_names = 'initial_damage'
    #     prop_values = 0.7
    #     block = '1 3'
    #     outputs = exodus
    # []
    # [initialdamage_block2]
    #     type = ADGenericConstantMaterial
    #     prop_names = 'initial_damage'
    #     prop_values = 0.0
    #     block = '2'
    #     outputs = exodus
    # []
    # [initialbreakage_block1]
    #     type = ADGenericConstantMaterial
    #     prop_names = 'initial_breakage'
    #     prop_values = 0.1
    #     block = '1 3'
    #     outputs = exodus
    # []
    # [initialbreakage_block2]
    #     type = ADGenericConstantMaterial
    #     prop_names = 'initial_breakage'
    #     prop_values = 0.0
    #     block = '2'
    #     outputs = exodus
    # []
[]  

[Functions]
[]

[Preconditioning]
    [smp]
      type = SMP
      full = true
    []
[]
  
[Executioner]
    type = Steady
    solve_type = NEWTON
    l_max_its = 30
    l_tol = 1e-6
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-8
    nl_max_its = 30
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    line_search = bt
[]  

[Outputs]
    exodus = true   
    #show = 'initial_damage xi_initial'
    csv = true
[]

#We assume the simulation is loaded with compressive pressure and shear stress
#Let's go through the boundary conditions
#basically we have 6 DOFs that needs to be fixed: translations along x y z,
#and rotation about x y z axis. The first corner point is fixed in all 3 directions (ux = 0, uy = 0, uz = 0),
#this means the following: (1) three translations are fixed (2) rotations about y axis is fixed (ux = 0, uz = 0) (3) rotation about z axis is fixed (ux = 0, uy = 0) [note (2)(3) are combined effect with second node]
#now we fix a middle point only for y and z directions (uy = 0, uz = 0), this means we fix the rotation about x axis (ux = 0)
#thus all the 6 DOFs are fixed
[BCs]
    [pressure_right]
        type = ADPressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        factor = 135e6
    []
    [pressure_left]
        type = ADPressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        factor = 135e6
    []
    [pressure_front]
        type = ADPressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        factor = 120e6
    []
    [pressure_back]
        type = ADPressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        factor = 120e6        
    []
    [pressure_top]
        type = ADPressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = top
        factor = 127.5e6         
    []
    [pressure_bottom]
        type = ADPressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = bottom
        factor = 127.5e6              
    []
    #
    [pressure_shear_front]
        type = ADNeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        value = 50e6
    []
    [pressure_shear_back]
        type = ADNeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        value = -50e6   
    []
    [pressure_shear_left]
        type = ADNeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        value = -50e6
    []
    [pressure_shear_right]
        type = ADNeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        value = 50e6     
    []
    #
    [fix_bottom_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = bottom
    []
    #
    [fix_bottomline_x]
        type = ADDirichletBC
        variable = disp_x
        value = 0
        boundary = bottomline
    []
    [fix_bottomline_z]
        type = ADDirichletBC
        variable = disp_z
        value = 0
        boundary = bottomline
    []
    # #
    # [fix_ptr_x]
    #     type = ADDirichletBC
    #     variable = disp_x
    #     value = 0
    #     boundary = corner_ptr
    # []
    # [fix_ptr_y]
    #     type = ADDirichletBC
    #     variable = disp_y
    #     value = 0
    #     boundary = corner_ptr
    # []
    # [fix_ptr_z]
    #     type = ADDirichletBC
    #     variable = disp_z
    #     value = 0
    #     boundary = corner_ptr
    # []
    # #
    # [fix_ptr2_y]
    #     type = ADDirichletBC
    #     variable = disp_y
    #     value = 0
    #     boundary = corner_ptr2
    # []
    # [fix_ptr2_z]
    #     type = ADDirichletBC
    #     variable = disp_z
    #     value = 0
    #     boundary = corner_ptr2
    # []
    # #
    # [fix_ptr3_z]
    #     type = ADDirichletBC
    #     variable = disp_z
    #     value = 0
    #     boundary = corner_ptr3
    # []
[]