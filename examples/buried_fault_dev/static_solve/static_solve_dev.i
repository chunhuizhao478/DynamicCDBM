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
                  10000  -10000   10000;'
        new_boundary = corner_ptr
        input = sidesets
    [] 
    [./extranodeset2]
        type = ExtraNodesetGenerator
        coord = '0  -10000  -10000;'
        new_boundary = corner_ptr2
        input = extranodeset1
        use_closest_node=true
    [] 
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
    CdCb_multiplier = 0

    #<coefficient of healing for breakage evolution>: refer to "Lyakhovsky_Ben-Zion_P14" (10 * C_B)
    # CBCBH_multiplier = 0.0
    CBH_constant = 0

    #<coefficient of healing for damage evolution>: refer to "ggw183.pdf"
    C_1 = 0

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
    chi = 0.8

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

[AuxKernels]
[]

[Materials]
    [strain]
        type = ComputeSmallStrain
        displacements = 'disp_x disp_y disp_z'
        outputs = exodus
    [] 
    [stress_nucleation]
        type = ComputeDamageBreakageStress3D
        alpha_grad_x = alpha_grad_x
        alpha_grad_y = alpha_grad_y
        alpha_grad_z = alpha_grad_z
        output_properties = 'B alpha_damagedvar xi'
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
        type = InitialDamageBenchmark
        nucl_center = '0 -5000 0'
        fault_plane = '-2500 2500 -7500 -2500 -500 500'
        nucl_distance = 400
        nucl_thickness = 200
        nucl_damage = 0.7
        e_damage = 0.7
        e_sigma = 2.5e2
        outputs = exodus
    []
    #we close the damage perturbation for the dynamic solve
    #we nucleate the rupture by increase the boundary shear loading
    #e_damage = 0
    [damage_perturb]
        type = DamagePerturbationSquare
        nucl_center = '0 -5000 0'
        e_damage = 0
        thickness = 200
        length = 400
        duration = 1e-1
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
    # solve_type = 'NEWTON'
    solve_type = 'PJFNK'
    # start_time = -1e-12
    # end_time = 1e100
    # num_steps = 1
    l_max_its = 100
    l_tol = 1e-7
    nl_rel_tol = 1e-8
    nl_max_its = 50
    nl_abs_tol = 1e-10
    petsc_options_iname = '-ksp_type -pc_type'
    petsc_options_value = 'gmres     hypre'
    # petsc_options_iname = '-pc_type -pc_factor_shift_type'
    # petsc_options_value = 'lu       NONZERO'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
    # petsc_options_value = '101                asm      lu'
    # petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    # petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    # nl_forced_its = 3
    # line_search = 'bt'
    # dt = 1e-8
    # steady_state_detection = true
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
        type = Pressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        factor = 135e6
    []
    [pressure_left]
        type = Pressure
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        factor = 135e6
    []
    [pressure_front]
        type = Pressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        factor = 120e6
    []
    [pressure_back]
        type = Pressure
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        factor = 120e6        
    []
    [pressure_top]
        type = Pressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = top
        factor = 127.5e6         
    []
    [pressure_bottom]
        type = Pressure
        variable = disp_y
        displacements = 'disp_x disp_y disp_z'
        boundary = bottom
        factor = 127.5e6              
    []
    #
    [pressure_shear_front]
        type = NeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = front
        value = 55e6
    []
    [pressure_shear_back]
        type = NeumannBC
        variable = disp_x
        displacements = 'disp_x disp_y disp_z'
        boundary = back
        value = -55e6   
    []
    [pressure_shear_left]
        type = NeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = left
        value = -55e6
    []
    [pressure_shear_right]
        type = NeumannBC
        variable = disp_z
        displacements = 'disp_x disp_y disp_z'
        boundary = right
        value = 55e6     
    []
    #
    [fix_ptr_x]
        type = DirichletBC
        variable = disp_x
        value = 0
        boundary = corner_ptr
    []
    [fix_ptr_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = corner_ptr
    []
    [fix_ptr_z]
        type = DirichletBC
        variable = disp_z
        value = 0
        boundary = corner_ptr
    []
    #
    [fix_ptr2_y]
        type = DirichletBC
        variable = disp_y
        value = 0
        boundary = corner_ptr2
    []
    [fix_ptr2_z]
        type = DirichletBC
        variable = disp_z
        value = 0
        boundary = corner_ptr2
    []
[]