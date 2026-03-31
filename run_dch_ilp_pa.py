# -*- coding: utf-8 -*-

"""     
    Before executing this code, make sure that the extant 
    files are already created by sim.py, and match with
    folders and files names in folder_file_pair in below.

"""

""" 
    For adjusting the number of clusters for each extant file change 
    the resolution parameter (gamma).
"""

from dch_ilp_pa import run_dch_ilp_pa
        
folder_file_pair = {
    "clusters_for_table_2": [f"nx_p=0.5_n={n}" for n in [50, 100, 250, 500, 1000]],
    "deneme": [f"nx_p=0.5_n={n}" for n in [12, 20, 50, 100 ]],
    "clusters_for_table_4": [f"nx_p=0.5_n={n}" for n in [200]],
}

time_limits = {
    "clusters_for_table_2": [4],
    "clusters_for_table-3": [1, 2, 4, 8, 16, 32, 64],
    "clusters_for_table_4": [1.7, 17, 51, 68, 85, 136],
}


for folder in folder_file_pair.keys():
    for file in folder_file_pair[folder]:
        time_limit_list = time_limits[folder]
        run_dch_ilp_pa(folder, file, time_limit_list, N=1, gamma=0.5, p=0.5)

