'''
Created on Feb 1, 2023

@author: jiadongc
'''
from synthesis_planning.synthesis_pathways import SynthesisPathways
from synthesis_planning.interfacial_pdplotter import InterReactions, Inter_PDPlotter

from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.analysis.phase_diagram import CompoundPhaseDiagram

target = "BaAl2Si2O8"

# get the optimal synthesis recipe for a target material
sp = SynthesisPathways(target,
                       exclude_reactants = ["O2"],
                       selected_reactions_to_csv = True)

# display the selected reactions
for reaction in sp.selected_reactions:
    reaction.display()

# Visualize interfacial reaction compound phase diagram for the optimal reaction
reaction = sp.selected_reactions[0]
interfacial_reactions = InterReactions(reaction)
Inter_PDPlotter(
    interfacial_reactions,
    emphasize_entries = [reaction.target]
    ).show()


