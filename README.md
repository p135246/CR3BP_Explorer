# CR3BP_Explorer

The goal of this project is to implement various algorithms related to space mission design in Mathematica 10.4 and create an easy to use interface primarily for demonstration purposes. 


Version from 08.12.18 --- Trajectories with prescribed itineraries in the CR3BP from [1, Chapter 4]
=====================

* Algorithms from [1, Chapter 4] to find *trajectories with prescribed itineraries in the restricted 3-body problem* are implemented in the library *CR3BP.wl*. The only missing part is the algorithm to compute Lyapunov orbits of energy E+dE from the Lyapunov orbit of energy E. The following is the reason why it should be implemented in a future release: 

  * The shooting algorithm with initial values coming from the linearized problem converges only for energies E slightly above the energy E(L) of the Lagrange point L. In order to get Lyapunov orbits at both L_1 and L_2 simultaneously for given E, we have to choose E>E(L_2)>E(L_1) slightly above E(L_2), so that both necks are open; however, depending on E(L_2)-E(L_1), it might happen that E is not slightly above E(L_1) anymore. Therefore, we have to find the Lyapunov orbit for an energy E(L_2)>E'>E(L_1) by shooting and then extend it to E = E'+dE by the algorithm in question. 

* The library *CR3BP.wl* is used by the Mathematica 10.4 notebook *explorer.nb*, which is a functioning prototype of the promised user interface. One can look for Lyapunov orbits, draw stable and unstable manifolds and search for their intersections at various Poincare sections. The explorer keeps a list of computed data whose visualization can be turned on and off. The following are some nice examples:
  * For \mu="Sun-Jupiter",  the Lagrange point L_2 and dE=0.0015, one can find 1:1 homoclinic transverse interesection of stable and unstable manifolds at U_4 in the exterior realm.

* KNOWN ERRORS
  * The definitions of U_2 and U_3 are messed up.
  * When no Poincare section is found for given conditions, the program gives an error (it cannot handle empty lists). 

* TODO
  * Write manual to *explorer.nb*
  * Implement the algorithm for Lyapunov orbits with higher energy.
  * Find some more examples, e.g. parameters E and \mu for and a complete homoclinic-heteroclinic chain!
  * Long term - Implement the methods of weak stability boundary from [2]

References
==========

[1] W.S. Koon, M.W. Lo, J.E. Marsden, S.D. Ross: Dynamical Systems, the Three-Body Problem
and Space Mission Design (link http://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBk.pdf)

[2] E. Belbruno, Capture Dynamics and Chaotic Motions in Celestial Mechanics