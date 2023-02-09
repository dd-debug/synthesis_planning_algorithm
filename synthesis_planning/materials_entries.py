'''
Created on Feb 1, 2023

@author: jiadongc@umich.edu
'''
import os 
import json

from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.ext.matproj import MPRester
# For people who use the new MaterialsProject API:
'''from mp_api import MPRester'''

from synthesis_planning import settings


def getOrigStableEntriesList(els,filename = None):
    """
    Guery stable entries from MaterialsProject databse based on chemical elements. 
    Save the entries in the first time query.
    Args:
        els (list): list of strings of chemical elements
        filename (string): path and filename of the saved entries file
    Return:
        list of queried entries
    """
    directory = os.getcwd()
    directory += "/data/entries_files"
    if not os.path.exists(directory):
        os.makedirs(directory)
    s_els = list(els).copy()
    
    s_els.sort()
    if filename == None:
        filename = '-'.join(s_els)

    cache = os.path.join(directory, filename)
    if os.path.exists(cache):
        print('loading from cache. stable entries','-'.join(s_els))
        with open(cache, 'r') as f:
            dict_entries = json.load(f)
        list_entries = []
        for e in dict_entries:
            list_entries.append(ComputedEntry.from_dict(e))
        return list_entries
    else:
        print('Reading from database.stable entries','-'.join(s_els))
        with MPRester(settings.MPI_KEY) as MPR:
            entries = MPR.get_entries_in_chemsys(s_els)
        pd = PhaseDiagram(entries)
        newentries=[]
        for e in pd.stable_entries:
            newentries.append(e)
        dict_entries = []
        for e in newentries:
            dict_entries.append(e.as_dict())
        with open(cache,'w') as f:
            json.dump(dict_entries,f)
        return newentries

def getEntriesList(els,filename = None):
    """
    Guery all entries from MaterialsProject databse based on chemical elements. 
    Save the entries in the first time query.
    Args:
        els (list): list of strings of chemical elements
        filename (string): path and filename of the saved entries file
    Return:
        list of queried entries
    """
    directory = os.getcwd()
    directory += "/data/entries_files"
    s_els = list(els).copy()
    s_els.sort()
    if filename == None:
        filename = '-'.join(s_els) + "_all_entries"
    cache = os.path.join(directory, filename)
    if os.path.exists(cache):
        print('loading from cache. all entries','-'.join(s_els))
        with open(cache, 'r') as f:
            dict_entries = json.load(f)
        list_entries = []
        for e in dict_entries:
            list_entries.append(ComputedEntry.from_dict(e))
        return list_entries
    else:
        print('Reading from database. all entries','-'.join(s_els))
        with MPRester(settings.MPI_KEY) as MPR:
            entries = MPR.get_entries_in_chemsys(s_els)

        dict_entries = []
        for e in entries:
            dict_entries.append(e.as_dict())
        with open(cache,'w') as f:
            json.dump(dict_entries,f)
        return entries