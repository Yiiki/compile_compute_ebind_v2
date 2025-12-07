compute binding energy as

εim = E6 – Ee + Vbulk – Vcorner

where E6 is the dopant eigen energy level

Ee is the reference band edge (Ec for donor, Ev for acceptor) in the pure bulk

Vcorner is the potential of the most distant mesh point to the dopant in the doped supercell

Vbulk is the corresponding potential in the pure bulk

Note, use single remote point is dangerous and unstable. 

On the one hand, Ecut2 requires high to get converge results from DFT calculations;

On the other hand, using approximated method, such as patching method would be destroyed due to the local error.

Therefore, a shell average method is adopted here.
