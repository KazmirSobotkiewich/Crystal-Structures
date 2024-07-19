# Crystal-Structures

[![python](https://img.shields.io/badge/Python-3.9-3776AB.svg?style=flat&logo=python&logoColor=white)](https://www.python.org)

The following are the highlights of my senior research thesis completed from 2019-2020 at McGill University in collaboration with the Siwick Research Group.

## Vanadium Dioxide (VO<sub>2</sub>)

At a temperature of ~ 340 K, the crystal VO<sub>2</sub> undergoes an insulator-metal transition, which is of signficiant interest in the application of dielectrics, semiconductors, as well as in magnetic and optical uses. The metal phase above 340 K is called the Rutile $R$ phase, which transitions into a semi-conductor when cooled below 340 K, called the Monoclinic $M$<sub>1</sub> phase.

The goal of this project is to create predictive models of VO<sub>2</sub> from theory in order to improve the satistical precision of experimental measurements in ultrafast laser spectroscopy.

## Creating 2D cross sections of VO<sub>2</sub>'s molecular structure

![alt text](https://github.com/KazmirSobotkiewich/Crystal-Structures/blob/main/Cross-Section-VO2-R.png) ![alt text](https://github.com/KazmirSobotkiewich/Crystal-Structures/blob/main/Cross-Section-VO2-M1.png)

These images show the molecular structure of VO<sub>2</sub> in the Rutile $R$ phase (left) and the Monoclinic $M$<sub>1</sub> phase (right). Vanadium atoms can be seen in red nearly perfectly octohedrally-coordinates in both phases, while Oxygen atoms can be seen in white. In the $M$<sub>1</sub> phase, V-V dimerization occurs which can be seen by noting the tilting of the red V atoms in alternating directions.

The code used to create these plots can be found in [ElectrostaticPotential_R.py](https://github.com/KazmirSobotkiewich/Crystal-Structures/blob/main/ElectrostaticPotential_R.py) and [ElectrostaticPotential_M1.py](https://github.com/KazmirSobotkiewich/Crystal-Structures/blob/main/ElectrostaticPotential_M1.py).

## Theory:

### Crystal lattice:

The crystal's structure is known by using the lattice vectors, the recriprocal lattice vectors, the location of the atoms using the lattice vectors as the basis and the elastic atomic scattering factors of electrons for the corresponding neutral atoms.

### Structure factor:

$`S(\textbf{G}) = \sum\limits_{j}f_j(\textbf{G})e^{-i(\textbf{G} \cdot \textbf{r}_j)}`$ 

$f_j$ are the atomic form factors and $`\textbf{G} = h\:\textbf{b}_1 + k\:\textbf{b}_2 + l\:\textbf{b}_3`$.

The Miller indices are integers labeled $h$, $k$ and $l$ which represent points in the reciprocal lattice basis, where $\textbf{b}_1$, $\textbf{b}_2$ and $\textbf{b}_3$ are the reciprocal lattice vectors.

### Electrostatic potential:

$`\phi (\textbf{r},t) = \sum\limits_{\{\textbf{G}\}} |S(\textbf{G},t)|e^{i \chi (\textbf{G},t)} \cos (\textbf{G} \cdot \textbf{r})`$

$\chi (\textbf{G},t)$ is the phase of the generally complex structure factor.
