# This file is the resource file for the python program 'branch.py'
# All entries must be valid python code

def set_parameters():

	upper_level = "a3F4p_y4F3*"
	
	file_name = "branches"
	reference_level = "a3F4s_b4F4"
	reference_file = "Fe011720b.003.I"
	normal_level = "d5____a4P1"
	
	identified_lines_file = "FeII.GN"
	level_file = "fe_ii.lev"
	lifetime_unc = 0.10
	spectrum_files = '*.II'

	return (file_name, upper_level, reference_level, reference_file,
			normal_level, identified_lines_file, level_file, lifetime_unc, spectrum_files)
