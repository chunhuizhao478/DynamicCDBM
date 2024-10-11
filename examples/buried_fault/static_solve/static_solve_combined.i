zmin = -10000.0
normal_traction_y = 127500000.0
shear_modulus_o = 30000000000.0
fault_plane_xmin = -2500.0
normal_traction_x = 135000000.0
ymin = -10000.0
fault_plane_ymin = -7500.0
fault_plane_xmax = 2500.0
nucl_thickness = 200.0
shear_traction = 55000000.0
zmax = 10000.0
nucl_damage = 0.7
e_sigma = 250.0
xi_0 = -0.8
fault_plane_ymax = -2500.0
normal_traction_z = 120000000.0
fault_plane_zmin = -500.0
xmax = 10000.0
xmin = -10000.0
e_damage = 0.3
lambda_o = 30000000000.0
nucl_center_y = -5000.0
nucl_center_z = 0.0
nucl_center_x = 0.0
fault_plane_zmax = 500.0
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
        lambda_o = ${lambda_o}
        shear_modulus_o = ${shear_modulus_o}
        xi_o = ${xi_0}
        outputs = exodus
    []
    [getxi]
        type = ADComputeXi
        outputs = exodus
    []
    [initialdamage]
        type = ADInitialDamageBenchmark
        nucl_center = '${nucl_center_x} ${nucl_center_y} ${nucl_center_z}'
        fault_plane = '${fault_plane_xmin} ${fault_plane_xmax} ${fault_plane_ymin} ${fault_plane_ymax} ${fault_plane_zmin} ${fault_plane_zmax}'
        nucl_distance = ${nucl_distance}
        nucl_thickness = ${nucl_thickness}
        nucl_damage = ${nucl_damage}
        e_damage = ${e_damage}
        e_sigma = ${e_sigma}
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
    l_max_its = 10
    l_tol = 1e-6
    nl_rel_tol = 1e-8
    nl_abs_tol = 1e-10
    nl_max_its = 10
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type  -ksp_initial_guess_nonzero -ksp_pc_side -ksp_max_it -ksp_rtol -ksp_atol'
    petsc_options_value = 'gmres        hypre      boomeramg                   True        right       1500        1e-7      1e-9    '
    automatic_scaling = true
    line_search = basic
[]  

[Outputs]
    exodus = true   
    #show = 'initial_damage xi_initial'
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