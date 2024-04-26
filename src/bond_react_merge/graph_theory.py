# -*- coding: utf-8 -*-
"""
@author: Josh Kemppainen
Revision 1.0
January 9th, 2023
Michigan Technological University
1400 Townsend Dr.
Houghton, MI 49931


All graph-theory functions are taken from here:
https://www.python.org/doc/essays/graphs/

Each funtion takes in a graph like this:
    
    graph = {'A': ['B', 'C'],
             'B': ['C', 'D'],
             'C': ['D'],
             'D': ['C'],
             'E': ['F'],
             'F': ['C']}
    
    # Which presents this directed graph
    A -> B -> D
    |    /   ^
    |   /   /
    | /   /
    V   /
    C <
    
Examples of Functions:
    
    print(find_path(graph, 'A', 'D'))
    ['A', 'B', 'C', 'D']
    
    
    print(find_all_paths(graph, 'A', 'D'))
    [['A', 'B', 'C', 'D'], ['A', 'B', 'D'], ['A', 'C', 'D']]
    
    print(find_shortest_path(graph, 'A', 'D'))
    ['A', 'B', 'D']
    
    
Only minor changes were made to each of the functions
in this file to allow them to work on python3. Otherwise
all other logic is given credit from the python website.
"""


# Function to find an arbitrary path
def find_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None


# Function to find all path's
def find_all_paths(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return [path]
    if start not in graph:
        return []
    paths = []
    for node in graph[start]:
        if node not in path:
            newpaths = find_all_paths(graph, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths


# Function to find shortest path
def find_shortest_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = find_shortest_path(graph, node, end, path)
            if newpath:
                if not shortest or len(newpath) < len(shortest):
                    shortest = newpath
    return shortest
