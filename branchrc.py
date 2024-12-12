# This file is the resource file for the python program 'branch.py'
# All entries must be valid python code

def set_parameters():
	
    file_name = "branches"
    upper_level = "a3F4p_y4F3*"
    reference_level = "a3F4s_b4F4"
    reference_file = "Fe011720b.003.I"
    normal_level = "d5____a4P1"
    identified_lines_file = "FeII.lines"
    level_file = "fe_ii.lev"
    lifetime_unc = 0.10
    spectrum_files = '*.II'
    calc_file = "FeII_waveno.E1"
    discrim = 0.2
    plotwin_length = 32

    return (file_name, upper_level, reference_level, reference_file,
			normal_level, identified_lines_file, level_file, lifetime_unc, spectrum_files,
            calc_file,discrim,plotwin_length)
