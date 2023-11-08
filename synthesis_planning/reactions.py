'''
Created on Feb 1, 2023

@author: jiadongc@umich.edu
'''

from itertools import combinations

from pymatgen.analysis.reaction_calculator import ComputedReaction

class SkipReaction(Exception):
    pass
   
def get_possible_reactions(precursors, product):
    """
    get all possible pairwise combinatorial reactions based on precursors and products
    Args:
        precursors (list): list of reactant entries
        product (list): list of product entries 
    Return:
        a list of reactions of pymatgen ComputedReaction objects
    """    
    combs = list(combinations(precursors, 2))
    reactions = []
    for reactants in combs:
        try:
            reactants = list(reactants)
            reaction = ComputedReaction(reactants, product)
            if len(reaction.reactants) == 2 and len(reaction.products)==1:
                reactions.append(reaction)
        except:
            SkipReaction('Reaction can not be compositionally balanced')
    return reactions

class Reaction():
    def __init__(self, target, reactants,
                 reaction_energy, inverse_hull_energy,
                 reaction, all_entries, all_phases):
        """
        Gather critical information of a chemical reaction
        Args:
            target (ComputedEntry): product of the reaction
            reactants (list): list of reactants entries of the reaction
            reaction_energy (float): reaction energy of the reaction
            inverse_hull_energy (float): inverse hull energy of the reaction
            reaction (ComputedReaction): a pymatgen ComputedReaction object of the reaction
            all_entries (list): list of ComputedEntry objects in the reaction compound 
                convex hull, including kinks (potential competing phases and decomposition 
                reactions). 
            all_phases (dict): {kink entry: decomposition entries at the kink} dictionary
        """
        self.target = target
        self.reactants = reactants
        self.reactE = reaction_energy
        self.invE = inverse_hull_energy
        self.reaction = reaction
        self.all_entries = all_entries
        self.all_phases = all_phases
        self.competing_phases_names = self.get_competing_phases_names()
    
    def get_competing_phases_names(self):
        """
        Get potential competing phases. Competing phases that will appear together 
        due to the decomposition reactions are in the same sub list. 
        e.g.,  for reaction LiBO2 + BaO -> BaLiBO3, there is one competing decomposition 
        reaction: LiBO2 + BaO -> Ba2Li(BO2)5 + Li6B4O9. So the return result will be
        [...,['Ba2Li(BO2)5', 'Li6B4O9'],...] 
        Return:
            list of list of strings.
        """
        competing_phases_names = []
        for entry in self.all_phases:
            cphases = []
            for j in self.all_phases[entry]:
                cphases.append(j.name)
            if len(cphases) == 1:
                if cphases[0] != self.target.name and\
                   cphases[0] not in [react.name for react in self.reactants]:
                    competing_phases_names += cphases
            else:
                competing_phases_names.append(cphases)
        return competing_phases_names
    
    def display(self):
        """display the reaction information"""
        print(self.__str__())
    
    def __repr__(self):
        """Return: a list of reaction information"""
        outputs = [
            self.target.name,
            [i.name for i in self.reactants],
            self.reactE,
            self.invE,
            self.reaction.__str__(),
            self.competing_phases_names
            ]
        return outputs
    
    def __str__(self):
        """Return: a string of reaction information"""
        outputs = [f"target: {self.target.name}",
        f"reactants: {[i.name for i in self.reactants]}",
        f"reaction energy: {self.reactE}",
        f"inverse hull energy: {self.invE}",
        f"reaction: {self.reaction.__str__()}",
        f"competing phases: {self.competing_phases_names}"
        ]
        return "\n"+"\n".join(outputs)
        