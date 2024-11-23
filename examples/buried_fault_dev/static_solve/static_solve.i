#explicit continuum damage-breakage model dynamics

[Mesh]
    [./msh]
        type = FileMeshGenerator
        file =  '../meshfile/cdbm_tpv2053d_buried_small.msh'
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
        coord = '-10000  -10000  -10000;
                  10000  -10000   10000'
        new_boundary = corner_ptr
        input = sidesets
    [] 
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
        outputs = exodus
    []
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
        nucl_center = '0 -5000 0'
        fault_plane = '-2500 2500 -7500 -2500 -500 500'
        nucl_distance = 400
        nucl_thickness = 200
        nucl_damage = 0.8
        e_damage = 0.7
        e_sigma = 2.5e2
        outputs = exodus
    [] 
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
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-10
    nl_max_its = 30
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    line_search = bt
[]  

[Outputs]
    exodus = true   
    #show = 'initial_damage xi_initial'
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
    [fix_ptr_x]
        type = ADDirichletBC
        variable = disp_x
        value = 0
        boundary = corner_ptr
    []
    [fix_ptr_y]
        type = ADDirichletBC
        variable = disp_y
        value = 0
        boundary = corner_ptr
    []
    [fix_ptr_z]
        type = ADDirichletBC
        variable = disp_z
        value = 0
        boundary = corner_ptr
    []
[]