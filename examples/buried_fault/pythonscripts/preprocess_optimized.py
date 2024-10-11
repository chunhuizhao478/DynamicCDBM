import os
import argparse
from typing import Dict, Set, List, Union
import numpy as np
from pathlib import Path
import gmsh
import math
import sys

class Preprocess:
    """Class to handle preprocessing of MOOSE input files."""
    
    def __init__(self, 
                 mesh_file: str,
                 params_file: str,
                 params_static_file: str,
                 params_dynamic_file: str, 
                 static_input_file: str, 
                 dynamic_input_file: str,
                 static_input_file_combined: str,
                 dynamic_input_file_combined: str):
        """Initialize the Preprocess class with file paths."""
        # Convert all file paths to Path objects for better path handling
        self.mesh_file = Path(mesh_file)
        self.params_file = Path(params_file)
        self.params_static_file = Path(params_static_file)
        self.params_dynamic_file = Path(params_dynamic_file)
        self.static_input_file = Path(static_input_file)
        self.dynamic_input_file = Path(dynamic_input_file)
        self.static_input_file_combined = Path(static_input_file_combined)
        self.dynamic_input_file_combined = Path(dynamic_input_file_combined)

        # Define parameter lists as frozen sets for immutability and faster lookups
        self.static_solve_list: Set[str] = frozenset([
            'xmin', 'xmax', 'ymin', 'zmin', 'zmax', 'lambda_o', 'shear_modulus_o', 'xi_0',
            'nucl_center_x', 'nucl_center_y', 'nucl_center_z', 'fault_plane_xmin', 'fault_plane_xmax',
            'fault_plane_ymin', 'fault_plane_ymax', 'fault_plane_zmin', 'fault_plane_zmax', 'nucl_distance',
            'nucl_thickness', 'nucl_damage', 'e_damage', 'e_sigma', 'normal_traction_x', 'normal_traction_z',
            'normal_traction_y', 'shear_traction'
        ])

        self.dynamic_solve_list: Set[str] = frozenset([
            'xmin', 'xmax', 'ymin', 'zmin', 'zmax', 'lambda_o', 'shear_modulus_o', 'xi_0', 'xi_d', 
            'Cd_constant', 'CdCb_multiplier', 'CBH_constant', 'C_1', 'C_2', 'beta_width', 'C_g', 
            'm1', 'm2', 'chi', 'rho', 'nucl_thickness', 'nucl_distance', 'duration', 'dt', 'end_time', 
            'time_step_interval', 'normal_traction_x', 'normal_traction_z', 'normal_traction_y', 
            'shear_traction', 'nucl_center_x', 'nucl_center_y', 'nucl_center_z', 'e_damage'
        ])

        # Parameter validation rules
        self.validation_rules = {
            'lambda_o': (0, float('inf')),
            'shear_modulus_o': (0, float('inf')),
            'xi_0': (-np.sqrt(3), np.sqrt(3)),
            'xi_d': (-np.sqrt(3), np.sqrt(3)),
            'Cd_constant': (0, float('inf')),
            'CdCb_multiplier': (0, float('inf')),
            'CBH_constant': (0, float('inf')),
            'C_1': (0, float('inf')),
            'C_2': (0, float('inf')),
            'beta_width': (0, float('inf')),
            'C_g': (0, float('inf')),
            'm1': (0, float('inf')),
            'm2': (0, float('inf')),
            'chi': (0, 1),
            'dt': (0, float('inf')),
            'end_time': (0, float('inf')),
            'time_step_interval': (0, float('inf')),
            'nucl_distance': (0, float('inf')),
            'nucl_thickness': (0, float('inf')),
            'e_damage': (0, float('inf')),
            'e_sigma': (0, float('inf'))
        }

    @staticmethod
    def read_file(file_path: Path) -> str:
        """Read and return the content of a file."""
        try:
            return file_path.read_text()
        except IOError as e:
            raise IOError(f"Error reading file {file_path}: {e}")

    def write_file(self, file_path: Path, content: str) -> None:
        """Write content to a file."""
        try:
            file_path.write_text(content)
            print(f"File created: {file_path}")
        except IOError as e:
            raise IOError(f"Error writing to file {file_path}: {e}")

    def combine_files(self) -> None:
        """Combine parameter files with their respective input files."""
        for params_file, input_file, output_file in [
            (self.params_static_file, self.static_input_file, self.static_input_file_combined),
            (self.params_dynamic_file, self.dynamic_input_file, self.dynamic_input_file_combined)
        ]:
            params_content = self.read_file(params_file)
            input_content = self.read_file(input_file)
            combined_content = f"{params_content}\n\n{input_content}"
            self.write_file(output_file, combined_content)

    def delete_files(self) -> None:
        """Delete all generated files if they exist."""
        for file_path in [self.params_file, self.params_static_file, self.params_dynamic_file,
                         self.static_input_file_combined, self.dynamic_input_file_combined]:
            if file_path.exists():
                file_path.unlink()
                print(f"File deleted: {file_path}")

    def validate_param_value(self, param_name: str, param_value: str) -> float:
        """Validate parameter values against defined rules."""
        try:
            value = float(param_value)
            if param_name in self.validation_rules:
                min_val, max_val = self.validation_rules[param_name]
                if not min_val < value < max_val:
                    raise ValueError(
                        f"The value of '{param_name}' must be between {min_val} and {max_val}")
            return value
        except ValueError as e:
            raise ValueError(f"Invalid value for {param_name}: {e}")

    def validate_geometry_values(self, params: Dict[str, float]) -> None:
        """Validate geometry parameters."""
        required_params = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']
        missing_params = set(required_params) - set(params.keys())
        if missing_params:
            raise ValueError(f"Missing required parameters: {', '.join(missing_params)}")

        for dim in ['x', 'y', 'z']:
            min_val = params[f'{dim}min']
            max_val = params[f'{dim}max']
            if min_val >= max_val:
                raise ValueError(f"{dim}min must be less than {dim}max")

    def calculate_fault_parameters(self, params: Dict[str, Union[str, float]]) -> Dict[str, float]:
        """Calculate fault plane parameters."""
        fault_params = {}
        dimensions = ['x', 'y', 'z']
        lengths = ['strike', 'dip', 'normal']
        
        for dim, length in zip(dimensions, lengths):
            center = float(params[f'nucl_center_{dim}'])
            length_val = float(params[f'len_fault_{length}'])
            half_length = 0.5 * length_val
            
            fault_params[f'fault_plane_{dim}min'] = center - half_length
            fault_params[f'fault_plane_{dim}max'] = center + half_length
            
        return fault_params

    def build_input_file(self, user_params: Dict[str, str]) -> Path:
        """Build the input file from user parameters."""
        # Validate all parameters
        validated_params = {
            name: self.validate_param_value(name, value)
            for name, value in user_params.items()
        }
        
        self.validate_geometry_values(validated_params)
        
        # Calculate additional parameters
        fault_params = self.calculate_fault_parameters(user_params)
        validated_params.update(fault_params)

        # Build mesh
        self.build_mesh(user_params)
        
        # Generate content
        content = '\n'.join(f"{name} = {value}" for name, value in validated_params.items())
        
        self.write_file(self.params_file, content)
        return self.params_file

    def generate_static_dynamic_files(self) -> None:
        """Generate static and dynamic solve input files."""
        params_content = self.read_file(self.params_file)
        
        # Parse parameters
        params_dict = dict(
            line.split('=', 1) for line in params_content.splitlines() 
            if '=' in line
        )
        params_dict = {k.strip(): v.strip() for k, v in params_dict.items()}
        
        # Generate content for each file type
        for solve_type, params_list, output_file in [
            ('static', self.static_solve_list, self.params_static_file),
            ('dynamic', self.dynamic_solve_list, self.params_dynamic_file)
        ]:
            content = []
            missing_params = []
            
            for param in params_list:
                if param in params_dict:
                    content.append(f"{param} = {params_dict[param]}")
                else:
                    missing_params.append(param)
            
            if missing_params:
                raise ValueError(f"Missing {solve_type} parameters: {', '.join(missing_params)}")
            
            self.write_file(output_file, '\n'.join(content))

    def args_parsing(self) -> None:
        """Parse command line arguments and build input file."""
        parser = argparse.ArgumentParser(description="Build MOOSE input file with parameters.")
        parser.add_argument("params", nargs='+', help="Parameter names and values in pairs")
        args = parser.parse_args()

        if len(args.params) % 2 != 0:
            raise ValueError("Parameters must be provided in name-value pairs")

        params_dict = dict(zip(args.params[::2], args.params[1::2]))
        self.build_input_file(params_dict)

    def create_fault_geometry(self,lc,
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

    def build_mesh(self, user_params):
        """Build mesh using the command line parameters"""
        self.create_fault_geometry(lc = float(user_params["lc"]),
                            lc_fault = float(user_params["lc_fault"]),
                            Fault_length = float(user_params["len_fault_strike"]),
                            Fault_width = float(user_params["len_fault_dip"]),
                            Fault_thickness = float(user_params["len_fault_normal"]),
                            R_nucl = float(user_params["nucl_distance"]),
                            thickness_nucl = float(user_params["nucl_thickness"]),
                            Xmax = float(user_params["xmax"]),
                            Xmin = float(user_params["xmin"]),
                            Ymax = float(user_params["ymax"]),
                            Ymin = float(user_params["ymin"]),
                            Zmax = float(user_params["zmax"]),
                            Zmin = float(user_params["zmin"]),
                            output_file = str(self.mesh_file), 
                            show_gui = False)

def main():
    """Main function to run the preprocessing pipeline."""
    base_dir = Path('./examples/buried_fault')
    mesh_dir = base_dir / 'meshfile'
    params_dir = base_dir / 'params'
    static_dir = base_dir / 'static_solve'
    dynamic_dir = base_dir / 'dynamic_solve'

    preprocess_builder = Preprocess(
        mesh_file=mesh_dir / 'buried_fault.msh',
        params_file=params_dir / 'params.i',
        params_static_file=params_dir / 'params_static_solve.i',
        params_dynamic_file=params_dir / 'params_dynamic_solve.i',
        static_input_file=static_dir / 'static_solve.i',
        dynamic_input_file=dynamic_dir / 'dynamic_solve.i',
        static_input_file_combined=static_dir / 'static_solve_combined.i',
        dynamic_input_file_combined=dynamic_dir / 'dynamic_solve_combined.i'
    )

    preprocess_builder.delete_files()
    preprocess_builder.args_parsing()
    preprocess_builder.generate_static_dynamic_files()
    preprocess_builder.combine_files()


if __name__ == "__main__":
    main()