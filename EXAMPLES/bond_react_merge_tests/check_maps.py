# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
June 2, 2026
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931
"""
##############################
# Import Necessary Libraries #
##############################
import glob
import os
import re




# Function to read edgeIDs and Equivalences from map file
def read_map_file(filename):
    edge_ids = []
    edge_types = {}

    equivalences = []
    equiv_types = {}

    section = None
    with open(filename, "r") as f:
        for line in f:
            stripped = line.strip()

            if not stripped:
                continue

            # Detect sections
            if stripped == "EdgeIDs":
                section = "edge"
                continue

            elif stripped.startswith("Equivalences"):
                section = "equiv"
                continue

            elif stripped.split()[0] in {"InitiatorIDs", "CreateIDs", "DeleteIDs"}:
                section = None
                continue

            # -----------------
            # EdgeIDs section
            # -----------------
            if section == "edge":
                try:
                    atom_id = int(stripped.split()[0])
                    edge_ids.append(atom_id)

                    if "#" in line:
                        comment = line.split("#", 1)[1].strip()
                        atom_type = comment.split()[0]
                        edge_types[atom_id] = atom_type

                except (ValueError, IndexError):
                    pass

            # ----------------------
            # Equivalences section
            # ----------------------
            elif section == "equiv":
                try:
                    parts = stripped.split()

                    pre_id = int(parts[0])
                    post_id = int(parts[1])

                    pair = (pre_id, post_id)

                    equivalences.append(pair)

                    if "#" in line:
                        comment = line.split("#", 1)[1]
                        m = re.search(r"types:\s*(\S+)\s*->\s*(\S+)", comment)

                        if m:
                            equiv_types[pair] = (m.group(1), m.group(2))

                except (ValueError, IndexError):
                    pass

    return edge_ids, edge_types, equivalences, equiv_types

# Function to check against validated map
def run_check(maps_to_check, skip_hs=True):
    # Tallies
    passed = 0
    failed = 0
    
    # Run all tests looking for multiple validated map files
    all_tests = {} # {(validated_files_string, map_to_check_string): {validated_file_string:(failed_edge_ids, failed_equivs, failed_types)}}
    for validated_files, map_to_check in maps_to_check.items():
        tmp = {}
        for validated in glob.glob(validated_files):
            if not os.path.isfile(validated):
                print('  ERROR validated map file does not exist: {}'.format(validated))
                continue
    
            if not os.path.isfile(map_to_check):
                print('  ERROR map file to check does not exist: {}'.format(map_to_check))
                continue
            
            # Read files
            valid_edge_ids, valid_edge_types, valid_equivs, valid_equiv_types = read_map_file(validated)
            check_edge_ids, check_edge_types, check_equivs, check_equiv_types = read_map_file(map_to_check)
            
            # Check for edge ids
            failed_edge_ids = [] # [edgeID1, edgeID2, ...]
            for edge_id in check_edge_ids:
                if edge_id not in valid_edge_ids:
                    failed_edge_ids.append( edge_id )
                    
            # Check for equivalences
            failed_equivs = [] # [(preID, postID), (preID, postID) ...]
            failed_types  = [] # [(preType, postType), (preType, postType), ...]
            for pair in check_equivs:
                t1, t2 = check_equiv_types.get(pair, ('', ''))
                if t1.startswith('h') and skip_hs: continue
                if t1.startswith('H') and skip_hs: continue
                if pair not in valid_equivs:
                    failed_equivs.append( pair )
                    failed_types.append( (t1, t2) )
                    
            tmp[validated] = (failed_edge_ids, failed_equivs, failed_types)
        all_tests[(validated_files, map_to_check)] = tmp
        
    # Use the closests test to log
    for (validated_files, map_to_check) in all_tests:
        tests = all_tests[(validated_files, map_to_check)]
        if not tests: continue
        best_validated = min(tests, key=lambda k: len(tests[k][1]))        
        failed_edge_ids, failed_equivs, failed_types = tests[best_validated]
        print('checking "{}" -> "{}" mapping.'.format(best_validated, map_to_check))
        
        # Check for edge ids
        if not failed_edge_ids:
            print('  EdgeID mapping is correct!')
        else:
            failed_ids = [(str(i)) for i in failed_edge_ids]
            print('  EdgeID mapping failed with these: [{}]'.format( ', '.join(failed_ids) ))
            
        # Check for equivalences
        if not failed_equivs:
            print('  Equivalence mapping is correct!')
        else:
            print('  Equivalences mapping failed with these:')
            for (i, j), (ti, tj) in zip(failed_equivs, failed_types):
                print('  * ({}, {}) # {} {}'.format(i, j, ti, tj))
        
        # Full test tally
        if not failed_equivs and not failed_edge_ids:
            passed += 1
        else:
            failed += 1
            
    return passed, failed


##################################################
### Perform some automated checks on map files ###
##################################################
if __name__ == "__main__":  
    # Tally to keep track of results
    passed, failed = 0, 0
    
    # Check test1 maps
    print('\n\n')
    print('**************')
    print('* test1 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test1/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test1/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test1/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test1/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    
    # Check test2 maps
    print('\n\n')
    print('**************')
    print('* test2 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test2/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test2/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test2/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test2/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test3 maps
    print('\n\n')
    print('**************')
    print('* test3 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test3/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test3/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test3/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test3/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test4 maps
    print('\n\n')
    print('**************')
    print('* test4 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test4/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test4/run/pre1-post1_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test5 maps
    print('\n\n')
    print('**************')
    print('* test5 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test5/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test5/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test5/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test5/run/pre2-post2_rxn-map_commented.txt'
    maps_to_check['test5/validated_map_files/pre3-post3_rxn-map*.txt'] = 'test5/run/pre3-post3_rxn-map_commented.txt'
    maps_to_check['test5/validated_map_files/pre4-post4_rxn-map*.txt'] = 'test5/run/pre4-post4_rxn-map_commented.txt'
    maps_to_check['test5/validated_map_files/pre5-post5_rxn-map*.txt'] = 'test5/run/pre5-post5_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test6 maps
    print('\n\n')
    print('**************')
    print('* test6 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test6/validated_map_files/DPET_rxn*.txt'] = 'test6/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test6/validated_map_files/DPEI_rxn*.txt'] = 'test6/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test7 maps
    print('\n\n')
    print('**************')
    print('* test7 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test7/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test7/run/pre1-post1_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test8 maps
    print('\n\n')
    print('**************')
    print('* test8 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test8/validated_and_large_files/pre1-post1_rxn-map*.txt'] = 'test8/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test8/validated_and_large_files/pre2-post2_rxn-map*.txt'] = 'test8/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test9 maps
    print('\n\n')
    print('**************')
    print('* test9 maps *')
    print('**************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test9/validated_and_large_files/pre1-post1_rxn-map*.txt'] = 'test9/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test9/validated_and_large_files/pre2-post2_rxn-map*.txt'] = 'test9/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test10 maps
    print('\n\n')
    print('***************')
    print('* test10 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test10/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test10/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test10/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test10/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    
    # Check test11 maps
    print('\n\n')
    print('***************')
    print('* test11 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test11/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test11/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test11/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test11/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test12 maps
    print('\n\n')
    print('***************')
    print('* test12 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test12/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test12/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test12/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test12/run/pre2-post2_rxn-map_commented.txt'
    maps_to_check['test12/validated_map_files/pre3-post3_rxn-map*.txt'] = 'test12/run/pre3-post3_rxn-map_commented.txt'
    maps_to_check['test12/validated_map_files/pre4-post4_rxn-map*.txt'] = 'test12/run/pre4-post4_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test13 maps
    print('\n\n')
    print('***************')
    print('* test13 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test13/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test13/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test13/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test13/run/pre2-post2_rxn-map_commented.txt'
    maps_to_check['test13/validated_map_files/pre3-post3_rxn-map*.txt'] = 'test13/run/pre3-post3_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test14 maps
    print('\n\n')
    print('***************')
    print('* test14 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test14/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test14/run/pre1-post1_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test15 maps
    print('\n\n')
    print('***************')
    print('* test15 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test15/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test15/run/pre1-post1_rxn-map_commented.txt'
    maps_to_check['test15/validated_map_files/pre2-post2_rxn-map*.txt'] = 'test15/run/pre2-post2_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test16 maps
    print('\n\n')
    print('***************')
    print('* test16 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test16/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test16/run/pre1-post1_rxn-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
    
    # Check test17 maps
    print('\n\n')
    print('***************')
    print('* test17 maps *')
    print('***************')
    maps_to_check = {} # {'validated_map_file':'map_file_to_check'}
    maps_to_check['test17/validated_map_files/pre1-post1_rxn-map*.txt'] = 'test17/run/rxn1-map_commented.txt'
    passed_check, failed_check = run_check(maps_to_check)
    passed += passed_check
    failed += failed_check
            
    
    # Final breakdown
    print('\n\n')
    print('*****************')
    print('* Final results *')
    print('*****************')
    total = passed + failed
    print('Passed : {} of {}'.format(passed, total))
    print('Failed : {} of {}'.format(failed, total))

