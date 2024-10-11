import os
import argparse
import numpy as np

class preprocess:

    def __init__(self, 
                 params_file,
                 params_static_file,
                 params_dynamic_file, 
                 static_input_file, 
                 dynamic_input_file,
                 static_input_file_combined,
                 dynamic_input_file_combined):
        
        self.params_file = params_file
        
        self.params_static_file = params_static_file
        self.params_dynamic_file = params_dynamic_file
        
        #don't delete these two
        self.static_input_file = static_input_file
        self.dynamic_input_file = dynamic_input_file
        
        self.static_input_file_combined = static_input_file_combined
        self.dynamic_input_file_combined = dynamic_input_file_combined

        # Lists of static and dynamic parameters based on existing param file
        self.static_solve_list = [
            'xmin', 'xmax', 'ymin', 'zmin', 'zmax', 'lambda_o', 'shear_modulus_o', 'xi_o',
            'nucl_center_x', 'nucl_center_y', 'nucl_center_z', 'fault_plane_xmin', 'fault_plane_xmax',
            'fault_plane_ymin', 'fault_plane_ymax', 'fault_plane_zmin', 'fault_plane_zmax', 'nucl_distance',
            'nucl_thickness', 'nucl_damage', 'e_damage', 'e_sigma', 'normal_traction_x', 'normal_traction_z',
            'normal_traction_y', 'shear_traction'
        ]

        self.dynamic_solve_list = [
            'xmin', 'xmax', 'ymin', 'zmin', 'zmax', 'lambda_o', 'shear_modulus_o', 'xi_0', 'xi_d', 'Cd_constant', 'CdCb_multiplier', 'CBH_constant', 'C_1', 'C_2',
            'beta_width', 'C_g', 'm1', 'm2', 'chi', 'rho', 'thickness', 'length', 'duration', 'dt',
            'end_time', 'time_step_interval', 'normal_traction_x', 'normal_traction_z',
            'normal_traction_y', 'shear_traction','nucl_center_x', 'nucl_center_y', 'nucl_center_z','e_damage'
        ]

    ##combine files##
    #-----------------------------------------------------------------------------#

    def read_file(self, file_path):
        """Reads and returns the content of a file."""
        with open(file_path, 'r') as file:
            return file.read()

    def combine_files(self):
        """Combines the content of params_file and original_file."""

        ## static
        #---------------------------------------------------------------------------------------#
        # Read contents of both files
        params_content_static = self.read_file(self.params_static_file)
        original_file_content_static = self.read_file(self.static_input_file)

        # Combine contents (params at the top, followed by original_file)
        combined_content_static = params_content_static + '\n\n' + original_file_content_static

        # Write the combined content to the output file
        with open(self.static_input_file_combined, 'w') as file:
            file.write(combined_content_static)
        
        print(f"Combined file created: {self.static_input_file_combined}")

        ## dynamic
        #---------------------------------------------------------------------------------------#
        # Read contents of both files
        params_content_dynamic = self.read_file(self.params_dynamic_file)
        original_file_content_dynamic = self.read_file(self.dynamic_input_file)

        # Combine contents (params at the top, followed by original_file)
        combined_content_dynamic = params_content_dynamic + '\n\n' + original_file_content_dynamic

        # Write the combined content to the output file
        with open(self.dynamic_input_file_combined, 'w') as file:
            file.write(combined_content_dynamic)
        
        print(f"Combined file created: {self.dynamic_input_file_combined}")        

    def delete_files(self):
        """Deletes the files if it exists."""
        if os.path.exists(self.params_file):
            os.remove(self.params_file)
            print(f"Combined file deleted: {self.params_file}")
        if os.path.exists(self.params_static_file):
            os.remove(self.params_static_file)
            print(f"Combined file deleted: {self.params_static_file}")
        if os.path.exists(self.params_dynamic_file):
            os.remove(self.params_dynamic_file)
            print(f"Combined file deleted: {self.params_dynamic_file}")
        if os.path.exists(self.static_input_file_combined):
            os.remove(self.static_input_file_combined)
            print(f"Combined file deleted: {self.static_input_file_combined}")
        if os.path.exists(self.dynamic_input_file_combined):
            os.remove(self.dynamic_input_file_combined)
            print(f"Combined file deleted: {self.dynamic_input_file_combined}")

    ## Argument parsing for generating param.i ##
    #-----------------------------------------------------------------------------#

    def validate_param_value(self, param_name, param_value):
        """Validates that the parameter value is a float and greater than zero."""
        try:
            value = float(param_value)
        except ValueError:
            raise ValueError(f"Error: The value of '{param_name}' must be a valid float")
        """check material parameter values"""
        if param_name == 'lambda_o' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}' must be greater than zero.")
        if param_name == 'shear_modulus_o' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}' must be greater than zero.")
        if param_name == 'xi_0' and ( value < -np.sqrt(3) or value > np.sqrt(3) ):
            raise ValueError(f"Error: The value of '{param_name}' must be in the range (-sqrt(3), sqrt(3)).")
        if param_name == 'xi_d' and ( value < -np.sqrt(3) or value > np.sqrt(3) ):
            raise ValueError(f"Error: The value of '{param_name}' must be in the range (-sqrt(3), sqrt(3)).")                
        if param_name == 'Cd_constant' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")        
        if param_name == 'CdCb_multiplier' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")          
        if param_name == 'CBH_constant' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'C_1' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'C_2' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'beta_width' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")        
        if param_name == 'C_g' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")          
        if param_name == 'm1' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'm2' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'chi' and ( value < 0 or value > 1 ):
            raise ValueError(f"Error: The value of '{param_name}'  must be in the range (0, 1).")
        if param_name == 'dt' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'end_time' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'time_step_interval' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        """check initial damamge geometry parameter values"""
        if param_name == 'nucl_distance' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'nucl_thickness' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")            
        if param_name == 'e_damage' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")
        if param_name == 'e_sigma' and value <= 0:
            raise ValueError(f"Error: The value of '{param_name}'  must be greater than zero.")  
        return value
    
    def validate_check_geometry_values(self, user_params):
        """Check min and max values are reasonable and ensure required parameters exist."""
        
        # List of required parameters
        required_params = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']
        
        # Check if all required parameters are in user_params
        for param in required_params:
            if param not in user_params:
                raise ValueError(f"Error: Missing required parameter '{param}'")
        
        # Perform the min-max validation checks
        if user_params['xmin'] >= user_params['xmax']:
            raise ValueError(f"Error: The value of xmin must be less than xmax.")
        
        if user_params['ymin'] >= user_params['ymax']:
            raise ValueError(f"Error: The value of ymin must be less than ymax.")
        
        if user_params['zmin'] >= user_params['zmax']:
            raise ValueError(f"Error: The value of zmin must be less than zmax.")
        
    def add_additional_params(self, user_params, content):

        #fault parameters
        ##calculate values
        user_params["fault_plane_xmin"] = float(user_params["nucl_center_x"]) - 0.5 * float(user_params["len_fault_strike"])
        user_params["fault_plane_xmax"] = float(user_params["nucl_center_x"]) + 0.5 * float(user_params["len_fault_strike"])
        user_params["fault_plane_ymin"] = float(user_params["nucl_center_y"]) - 0.5 * float(user_params["len_fault_dip"])
        user_params["fault_plane_ymax"] = float(user_params["nucl_center_y"]) + 0.5 * float(user_params["len_fault_dip"])
        user_params["fault_plane_zmin"] = float(user_params["nucl_center_z"]) - 0.5 * float(user_params["len_fault_normal"])
        user_params["fault_plane_zmax"] = float(user_params["nucl_center_z"]) + 0.5 * float(user_params["len_fault_normal"])

        ##write in content
        content += f"fault_plane_xmin = {user_params['fault_plane_xmin']}\n"
        content += f"fault_plane_xmax = {user_params['fault_plane_xmax']}\n"
        content += f"fault_plane_ymin = {user_params['fault_plane_ymin']}\n"
        content += f"fault_plane_ymax = {user_params['fault_plane_ymax']}\n"
        content += f"fault_plane_zmin = {user_params['fault_plane_zmin']}\n"
        content += f"fault_plane_zmax = {user_params['fault_plane_zmax']}\n"

        return content, user_params

    def build_input_file(self, user_params):
        """
        Creates a new file with parameters formatted as 'param_name = value'.
        
        Note:the parameters here can contain additional ones that is not needed
        to generate static and dynamic files      

        """
        content = ""

        # Create lines like 'param_name = value'
        for param_name, param_value in user_params.items():
            # Validate each parameter value
            validated_value = self.validate_param_value(param_name, param_value)
            content += f"{param_name} = {validated_value}\n"

        # Check reasonable geometry ranges
        self.validate_check_geometry_values(user_params)

        # Add additional parameters (from provided length)
        content, user_params = self.add_additional_params(user_params=user_params,
                                                 content=content)

        # Write the new content to the output file
        with open(self.params_file, 'w') as file:
            file.write(content)

        return self.params_file

    def args_pasing(self):
        "Process MOOSE input file with command-line parameters."
        parser = argparse.ArgumentParser(description="Build MOOSE input file with command-line parameters.")
        
        # Add parameters dynamically, for now assuming we have several parameters
        parser.add_argument("params", nargs='+', help="Parameter names and values")

        # Parse the arguments
        args = parser.parse_args()

        # Expecting the format of params as pairs of (param_name, param_value)
        if len(args.params) % 2 != 0:
            raise ValueError("Parameters and values should be provided in pairs.")

        # Create a dictionary from the provided params
        user_params = {args.params[i]: args.params[i + 1] for i in range(0, len(args.params), 2)}

        # Instantiate the builder and generate the new input file
        new_file_path = self.build_input_file(user_params)

        print(f"New input file created at: {new_file_path}")

    ## generate param_static_solve.i and param_dynamic_solve.i using param.i ##
    #-----------------------------------------------------------------------------#
    def generate_static_dynamic_files(self):
        """
        Generates param_static_solve.i and param_dynamic_solve.i using params.i.
        
        Note: the parameters here are needed to match self.static_solve_list and 
        self.dynamic_solve_list exactly.
        """
        all_params = self.read_file(self.params_file)

        static_content = ""
        dynamic_content = ""

        assigned_static_params = set()
        assigned_dynamic_params = set()

        # Generate static and dynamic solve input files
        for line in all_params.splitlines():
            if "=" in line:
                param_name, param_value = line.split('=')
                param_name = param_name.strip()

                if param_name in self.static_solve_list:
                    static_content += f"{param_name} = {param_value.strip()}\n"
                    assigned_static_params.add(param_name)

                if param_name in self.dynamic_solve_list:
                    dynamic_content += f"{param_name} = {param_value.strip()}\n"
                    assigned_dynamic_params.add(param_name)

        # Check for missing parameters
        missing_static_params = set(self.static_solve_list) - assigned_static_params
        missing_dynamic_params = set(self.dynamic_solve_list) - assigned_dynamic_params

        if missing_static_params:
            raise ValueError(f"Missing static parameters: {', '.join(missing_static_params)}")

        if missing_dynamic_params:
            raise ValueError(f"Missing dynamic parameters: {', '.join(missing_dynamic_params)}")

        # Write to static solve file
        with open(self.params_static_file, 'w') as file:
            file.write(static_content)
        print(f"Static solve file generated: {self.params_static_file}")

        # Write to dynamic solve file
        with open(self.params_dynamic_file, 'w') as file:
            file.write(dynamic_content)
        print(f"Dynamic solve file generated: {self.params_dynamic_file}")
    #-----------------------------------------------------------------------------#

if __name__ == "__main__":

    # Define the file paths
    params_file = '../params/params.i'

    params_static_file = '../params/params_static_solve.i'
    static_input_file = '../static_solve/static_solve.i'
    static_input_file_combined = '../static_solve/static_solve_combined.i'

    params_dynamic_file = '../params/params_dynamic_solve.i'
    dynamic_input_file = '../dynamic_solve/dynamic_solve.i'
    dynamic_input_file_combined = '../dynamic_solve/dynamic_solve_combined.i'

    # execute
    preprocess_builder = preprocess(params_file=params_file,
                                    params_static_file=params_static_file,
                                    params_dynamic_file=params_dynamic_file,
                                    static_input_file=static_input_file,
                                    dynamic_input_file=dynamic_input_file,
                                    static_input_file_combined=static_input_file_combined,
                                    dynamic_input_file_combined=dynamic_input_file_combined)
    preprocess_builder.delete_files()
    preprocess_builder.args_pasing()
    preprocess_builder.generate_static_dynamic_files()
    preprocess_builder.combine_files()
