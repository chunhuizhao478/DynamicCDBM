import gmsh
import math
import sys

def create_fault_geometry(lc,
                          lc_fault,
                          Fault_length,
                          Fault_width,
                          Fault_thickness,
                          R_nucl,
                          thickness_nucl,
                          Xmax,
                          Xmin,
                          Ymax,
                          Ymin,
                          Zmax,
                          Zmin,
                          output_file: str = "fault.msh", 
                          show_gui: bool = False):
    """
    Create fault geometry following the exact structure of the original .geo file.
    
    Args:
        output_file: Name of the output mesh file
        show_gui: Whether to display the Gmsh GUI
    """
    # Initialize Gmsh
    gmsh.initialize()
    
    # Create a new model
    gmsh.model.add("fault_model")
    
    # Define parameters
    lc = lc
    lc_fault = lc_fault
    Fault_length = Fault_length
    Fault_width = Fault_width
    Fault_thickness = Fault_thickness
    Fault_dip = math.pi/2

    # Nucleation parameters
    R_nucl = R_nucl
    thickness_nucl = thickness_nucl
    
    # Domain bounds
    Xmax = Xmax
    Xmin = Xmin
    Ymax = Ymax
    Ymin = Ymin
    
    # Derived parameters
    Width_nucl = Ymin / 2
    Zmin = Zmin
    Zmax = Zmax
    
    # Create boxes
    box1 = gmsh.model.occ.addBox(Xmin, 0, Zmin, 2 * Xmax, Ymin, 2 * Zmax)
    box2 = gmsh.model.occ.addBox(-Fault_length/2, Ymin/2 - Fault_width/2, -Fault_thickness/2, 
                                 Fault_length, Fault_width, Fault_thickness)
    box3 = gmsh.model.occ.addBox(-R_nucl/2, Width_nucl - R_nucl/2, -thickness_nucl/2, 
                                 R_nucl, R_nucl, thickness_nucl)
    box4 = gmsh.model.occ.addBox(-Fault_length/2, Ymin/2 - Fault_width/2, -thickness_nucl/2, 
                                 Fault_length, Fault_width, thickness_nucl)
    
    # Boolean operations
    volumes, _ = gmsh.model.occ.fragment([(3, box1), (3, box2), (3, box3), (3, box4)], [])
    
    # Synchronize
    gmsh.model.occ.synchronize()
    
    # Define mesh fields
    field_tag = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(field_tag, "VIn", lc_fault)
    gmsh.model.mesh.field.setNumber(field_tag, "VOut", lc / 4)
    gmsh.model.mesh.field.setNumber(field_tag, "XMin", -Fault_length/2)
    gmsh.model.mesh.field.setNumber(field_tag, "XMax", Fault_length/2)
    gmsh.model.mesh.field.setNumber(field_tag, "YMin", Ymin/2 - Fault_width/2)
    gmsh.model.mesh.field.setNumber(field_tag, "YMax", Ymin/2 + Fault_width/2)
    gmsh.model.mesh.field.setNumber(field_tag, "ZMin", -Fault_thickness/2)
    gmsh.model.mesh.field.setNumber(field_tag, "ZMax", Fault_thickness/2)
    
    # Set background mesh
    gmsh.model.mesh.field.setAsBackgroundMesh(field_tag)
    
    # Physical volumes
    all_volumes = gmsh.model.getEntities(3)
    for i, (dim, tag) in enumerate(all_volumes):
        physical_tag = gmsh.model.addPhysicalGroup(dim, [tag])
        gmsh.model.setPhysicalName(dim, physical_tag, f"Volume_{i+1}")
    
    print(f"Number of volumes created: {len(all_volumes)}")

    # Mesh generation
    gmsh.model.mesh.generate(3)
    
    # Write mesh to file
    gmsh.write(output_file)
    
    # Show GUI if requested
    if show_gui:
        if '-nopopup' not in sys.argv:
            gmsh.fltk.run()
    
    # Finalize Gmsh
    gmsh.finalize()
