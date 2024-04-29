# Synthesis planning algorithm

## Overview

Synthesis is a vital component of computational materials discovery. While high-throughput computation accelerates the identification of new ‘stable’ materials with functional properties, the actual realization of these materials is limited by their synthesis. This synthesis planning algorithm offer a physics motivated way to determine optimal synthesis recipes.

Our algorithm introduces a conceptual description of the convex hull to navigate optimal reaction pathways. The overarching principle is to identify precursors that save substantial reaction energy for the process from competing phases to target products, while avoiding low-energy geometrical subjects in the convex hull that may represent impurities or decomposition byproducts. 

Based on our algorithms, over 2000 recipes are high-throughput generated for potential high-component battery cathodes and electrolytes, such as Li/Na/K-based 4-component phosphates, borates, and redox-active/non-active transition metal oxides. We validate our theoretical framework with an automated robotic laboratory.

## Prerequisites

### Pymatgen

This algorithm has a dependency on `pymatgen` package of the Materials Project database using Python 3. You can install `pymatgen` by either

1. install the required packages in requirements.txt

   ```bash
   pip install -r requirements.txt
   ```

2. Go [here](https://pymatgen.org/installation.html) and follow the instructions to install your `pymatgen`.

### Pymatgen API Key

To use this algorithm, you need to generate an API key. This algorithm is using the legacy Materials Project API by default, but you can switch to new Materials Project API if needed.

- Go [here](https://legacy.materialsproject.org/open) to get a legacy Materials Project API
- Go [here](https://next-gen.materialsproject.org/api) to get a new Materials Project API

After you get a API Key, copy it and go to `synthesis_planning.settings` python file,  paste its string to `MPI_KEY` global variable:

```python
MPI_KEY = 'Your Materials Project API key'
```

## Tutorials and examples

Please run `BaLiBO3_example.py` to see how to use this algorithm. You use `synthesis_pathways` module to predict selected optimal reactions. Then use `interfacial_pdplotter` module to visualize and analyze reaction compound phase diagram.

## Citation

This synthesis planning algorithm is created by Jiadong Chen, Wenhao Sun in University of Michigan. If you use this synthesis planning algorithm, we kindly ask you to cite the following publication:

**Chen, J., Cross, S. R., Miara, L. J., Cho, J. J., Wang, Y., & Sun, W. (2024). Navigating phase diagram complexity to guide robotic inorganic materials synthesis. *Nature Synthesis*, 1-9.**