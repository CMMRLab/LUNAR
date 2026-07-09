# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
July 9, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import traceback


# Function to write all logged array processing data to a csv file
def output(filename, array_data):
    # Go through and find all unique wildcards and keys
    all_wildcards = set(); all_other_keys = set()
    for file, output_dict in array_data.items():
        wildcards = output_dict.get('wildcards', [])
        all_wildcards.add( len(wildcards) )
        
        # Get all other keys
        for key in output_dict:
            if key == 'wildcards': continue
            all_other_keys.add( key )

    # Start making unique data structure to hold all possible combinations
    # of wildards and outupts - initialize all with None's in case value
    # is missing.
    max_wildards = max(all_wildcards)
    columns = ['wildcards[{}]'.format(i) for i in range(max_wildards)] # [column names for csv file]
    try:
        from natsort import natsorted
        all_other_keys = natsorted(list(all_other_keys))
    except: all_other_keys = sorted(list(all_other_keys))
    columns.extend( all_other_keys )
    csv_data = {} # {'filename':{'wildcards[0]':'wildcard', 'FWHM':value, ... Ncomputes},  ... Nfiles}
    for file in array_data:
        csv_data[file] = {i:None for i in columns}
        
    # Go through and update None's to actually values
    for file, output_dict in array_data.items():
        
        # Upate wildcards
        wildcards = output_dict.get('wildcards', [])
        for n, wildcard in enumerate(wildcards):
            name = 'wildcards[{}]'.format(n)
            csv_data[file][name] = wildcard
            
        # Get all other keys
        for key, value in output_dict.items():
            if key == 'wildcards': continue
            name = key
            csv_data[file][name] = value
            
    # Write results to csv file
    try:
        with open(filename+'.csv', 'w') as f:
            # Join with comma's and write titles
            titles = ', '.join(['filename'] + columns)
            f.write('{}\n'.format(titles))
            
            # Write rows
            for file in csv_data:
                data = csv_data[file]
                outputs = [str(file)]
                for column in columns:
                    outputs.append(str(data[column]))
                f.write('{}\n'.format(', '.join(outputs)))
    except:
        print(f'ERROR - could not write {filename}.csv. Likely open in another program.')
        stack_trace_string = traceback.format_exc()
        print(f'\nAn error occurred: {stack_trace_string}')
    return


# Test writing the results out
if __name__ == "__main__": 
    array_file = 'free_volume_array_logging'
    array_data = {'file1': {'wildcards': [1],      'lx': 8.8250066741959, 'ly': 10.81838029137422, 'lz': 4.303315176034795, 'density': 0.7205598974470918, 'cell_volume': 410.847303993000, 'atom_volume': 332.8623990684031, 'free_volume': 77.9849049245973, 'percent_free_volume': 18.98148148148148}, 
                  'file2': {'wildcards': [2, 'I'], 'lx': 8.6985220886048, 'ly': 9.302360895942131, 'lz': 4.886883752006046, 'density': 0.7486517877835156, 'cell_volume': 395.430954767472, 'atom_volume': 337.8544695952852, 'free_volume': 57.5764851721868, 'percent_free_volume': 14.56043956043956},
                  'file3': {'wildcards': [3],      'lx': 8.1516182650192, 'ly': 8.919329635986871, 'lz': 5.783466685719964, 'density': 0.7040220195629541, 'cell_volume': 420.498340968635, 'atom_volume': 367.9360483475561, 'free_volume': 52.5622926210794, 'percent_free_volume': 12.50000000000000},
                  'file4': {'wildcards': [4, 1.1], 'lx': 6.5808990175037, 'ly': 9.129101893126181, 'lz': 5.600504702304101, 'density': 0.8798529247962354, 'cell_volume': 336.465428355733, 'atom_volume': 324.2433935650339, 'free_volume': 12.2220347906997, 'percent_free_volume': 3.632478632478632}}
    
    output(array_file, array_data)