'''
Created on Feb 1, 2023

@author: jiadongc@umich.edu
'''
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, CompoundPhaseDiagram

import os
import math
import pandas as pds

from synthesis_planning.materials_entries import getOrigStableEntriesList
# from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from synthesis_planning.reactions import get_possible_reactions, Reaction

class SynthesisPathways():
    def __init__(self, formula : str,
                 exclude_reactants = [], 
                 selected_reactions_to_csv = True,
                 normalize_to_pretty_formula = True):
        """
        Predict optimal synthesis pathways for a target by minimizing reaction
        compound convex hull complexity and save most driving force for the process 
        from competing phases to the target.

        Args:
            formula (str): string of the target formula
            exclude_reactants (list): list of strings of formulas to avoid in 
                reactants. Reactions that have reactants in this list will not 
                be considered.
            selected_reactions_to_csv (bool / string): whether save predicted reactions
                to a csv file or not. bool input type corresponds to default filenames,
                while string input type corresponds to designated filename.
            normalize_to_pretty_formula (bool): whether change the target formula
                to a pretty integer formula, i.e., Ba0.3Li0.3B0.3O0.9 changes to BaLiBO3
        """
        if normalize_to_pretty_formula:
            formula = norm_formula(formula)
        self.formula = formula
        self.exclude_reactants = exclude_reactants
        els = [str(el) for el in Composition(formula).elements]
        entries = getOrigStableEntriesList(els)
        self.pd = PhaseDiagram(entries)
        
        target = self.remove_original_target_entry_if_exists_in_hull(entries, formula)
        if not target:
            # If target is not a stable material but promising to synthesize, 
            # its should be near the convex hull, representing a near-hull 
            # metastable stability. We make a fake entry for the target with 
            # a slightly below the convex hull free energy for simplicity. 
            # This allows us to do estimate to target materials without structures.
            target = self.make_stable_entry_from_comp(Composition(formula))
        self.target = target
        
        self.reactions = get_possible_reactions(entries, [target])

        print(f"all possible pairwise reactions: {len(self.reactions)}")
        selected_reactions = self.get_target_deepest_reactions()
        self.all_pairwise_reactions = self.get_all_pairwise_reactions_info()
        self.selected_reactions = self.save_most_energy_for_last_step(selected_reactions)
        
        if selected_reactions_to_csv:
            self.turn_to_csv(selected_reactions_to_csv)
    
    def turn_to_csv(self, filename):
        """
        Save selected predicted reactions to default / designated csv file.
        Args:
            filename (bool / string)
        """
        self.df = pds.DataFrame([react.__repr__() for react in self.selected_reactions], 
                               columns = ['target', 'reactants',
                                          'reaction_energy(eV/atom)',
                                          'inverse_hull_energy(eV/atom)',
                                          'reaction',
                                          'competing phases'])
        if not isinstance(filename, str):
            directory = os.getcwd()
            directory += "/data/results_files"
            if not os.path.exists(directory):
                os.mkdir(directory)
            filename = self.formula + "_result.csv"
            self.df.to_csv(directory + "/" + filename)
        else:
            self.df.to_csv(filename)
    
    
    def save_most_energy_for_last_step(self, selected_reactions):
        """
        Sort the selected predicted reactions based on their inverse hull energies.
        A more negative inverse hull energy represents the driving force from competing
        phases to the target is greater.
        Args:
            selected_reactions (list): list of Reaction objects
        Return:
            list of sorted Reaction objects
        """
        if selected_reactions:
            selected_reactions = sorted(selected_reactions, 
                                        key = lambda x:x.invE)
        return selected_reactions

    def get_target_deepest_reactions(self):
        """
        Select reactions where the target is the deepest entry in the reaction compound
        convex hull. This means the reaction driving force from reactants to the target 
        is the largest.
        Return:
            list of Reaction objects0
        """
        pd = self.pd
        selected_reactions = []
        for reaction in self.reactions:

            reactants = reaction._reactant_entries
            # we only consider pairwise reactions
            comp1 = reactants[0].composition
            comp2 = reactants[1].composition
            # exclude reactions with specific reactants
            if any(reactant.name in self.exclude_reactants for reactant in reactants):
                continue
            new_entries, products = self.construct_kinks_entries(pd, comp1, comp2, self.target)

            cpd = CompoundPhaseDiagram(new_entries,[comp1,comp2])
            lowest_entry, depth = self.get_lowest_entry_and_energy(cpd)

            if norm_formula(lowest_entry.name) == norm_formula(self.target.name):
                invE = get_inverse_hull_energy(lowest_entry, cpd)
                selected_reactions.append(
                    Reaction(
                         self.target,
                         reactants, 
                         depth, 
                         invE,
                         reaction,
                         new_entries,
                         products,                         
                    ))
                
        return selected_reactions
    
    def get_all_pairwise_reactions_info(self):
        pd = self.pd
        all_pairwise_reactions = []
        for reaction in self.reactions:
            # print(reaction)
            reactants = reaction._reactant_entries
            # we only consider pairwise reactions
            comp1 = reactants[0].composition
            comp2 = reactants[1].composition
            # exclude reactions with specific reactants
            if any(reactant.name in self.exclude_reactants for reactant in reactants):
                continue
            new_entries, products = self.construct_kinks_entries(pd, comp1, comp2, self.target)
            
            cpd = CompoundPhaseDiagram(new_entries,[comp1,comp2])
            depth = reaction.calculated_reaction_energy/Composition(self.target.name).num_atoms
            # if norm_formula(lowest_entry.name) == norm_formula(self.target.name):
            
            for e in cpd.stable_entries:
                if norm_formula(e.name) == norm_formula(self.target.name):
                    target_entry = e
            invE = get_inverse_hull_energy(target_entry, cpd)
            all_pairwise_reactions.append(
                Reaction(
                     self.target,
                     reactants, 
                     depth, 
                     invE,
                     reaction,
                     new_entries,
                     products,                         
                ))
                
        return all_pairwise_reactions
    
    def remove_original_target_entry_if_exists_in_hull(self, entries, formula):
        """
        Remove the entry of the formula if it is a stable material in the 
        MaterialsProject database
        Args: 
            entries (list): list of ComputedEntry objects
            formula (string): string of a material's formula
        Return:
            whether the material is stable. If stable, return the target entry, 
            else, return False
        """
        target = False
        for entry in entries:
            if norm_formula(entry.name) == norm_formula(formula):
                target = entry
        if target:
            entries.remove(target)
        return target
    
    def make_stable_entry_from_comp(self, comp):
        """
        Make a fake stable entry that is slightly below the current convex hull
        Args:
            comp (Composition): a pymatgen Composition object
        Return:
            a fake convex hull stable ComputedEntry
        """
        
        energy = self.pd.get_hull_energy(comp)
        new_entry=ComputedEntry(comp, energy-0.001)
        return new_entry
    
    def construct_kinks_entries(self, pd, comp1, comp2, target):
        """
        Transfer kinks to entries.The kinks are intersections of the 'comp1-comp2' 
        reaction compound convex hull slice plane with tie lines or equilibrium 
        phases along the higher dimensional convex hull 'pd'. These kinks are 
        potential competing phases or decomposition reactions during the synthesis.
        Args:
            pd (PhaseDiagram): a pymatgen PhaseDiagram object that formed by all
                elements in comp1 and comp2
            comp1 (Composition): a pymatgen Composition object of a precursor
            comp2 (Composition): a pymatgen Composition object of a precursor
            target (ComputedEntry): a pymatgen ComputedEntry the target
        Return:
            new_entries (list): a list of kinks entries
            products (dict): {kink entry: decomposition entries at the kink} dictionary
        """

        cricomps = pd.get_critical_compositions(comp1, comp2)

        new_entries = []
        products = {}
        for comp in cricomps:
            energy = pd.get_hull_energy(comp)
            entry = ComputedEntry(comp, energy)
            products[entry] = list(self.pd.get_decomposition(comp).keys())
            new_entries.append(entry)
        new_entries.append(target) # in case the target is not in pd
        products[target] = [target]
        return new_entries, products
    
    def get_lowest_entry_and_energy(self, cpd):
        """
        Select the entry that has the most negative reaction energy versus the 
        terminal entries in the compound convex hull.
        Args:
            cpd (CompoundPhaseDiagram): a pymatgen CompoundPhaseDiagram object
                elements in comp1 and comp2
        Return:
            lowest_entry (ComputedEntry): a pymatgen ComputedEntry
            depth (float): reaction energy of the lowest_entry
        """
        depth = math.inf
        for e in cpd.stable_entries: #not new_entries
            ener = cpd.get_form_energy_per_atom(e)
            if ener < depth:
                depth = ener
                lowest_entry = e
        return lowest_entry,depth
    
def get_inverse_hull_energy(cpd_entry, cpd):
    """
    get the inverse hull energy of the entry. Inverse hull energy refers to the
    reaction energy of the entry versus to its nearby composition entries in the
    compound convex hull 'cpd'. 
    Args:
        cpd_entry (TransformedPDEntry / ComputedEntry)
        cpd (CompoundPhaseDiagram / PhaseDiagram)
    Return:
        invE (float): inverse hull energy of cpd_entry
    """

    mod_entries = []
    for e in cpd.stable_entries:
        if e.name != cpd_entry.name:
            mod_entries.append(e)
    mod_cpd = PhaseDiagram(mod_entries,cpd.species_mapping.values())

    invE = mod_cpd.get_decomp_and_e_above_hull(cpd_entry, 
                                               allow_negative=True)[1]
    return invE
    
def norm_formula(formula):
    """
    Transform the formula to a pretty integer formula, such as 'LiO0.5' to 'Li2O'
    Args:
        formula (string): 
    Return:
        string of a pretty integer formula
    """
    return Composition(formula).get_integer_formula_and_factor()[0]
    
    
    
    
    
    
    
    
    
    
    