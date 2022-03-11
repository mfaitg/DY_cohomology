# DY_cohomology
This repository contains the GAP programs used in section 5.4 of the paper "Davydov-Yetter cohomology and relative homological algebra" by M. Faitg, A.M. Gainutdinov and C. Schweigert, [arXiv:2202.12287](https://arxiv.org/abs/2202.12287). We wrote these programs in order to compute the 3rd and 4th Davydov-Yetter (DY) cohomology groups of the category of finite-dimensional modules over the restricted quantum group of sl(2) at a 4th root of unity, which is computationally too difficult to be done by hand.

The idea is just to encode all the objects into linear-algebraic data so that the problems are reduced to solve (quite big) systems of linear equations. Note that one could have used any software implementing linear algebra to do this; we chose GAP4 because of our familiarity with it.

The programs are extensively commented with detailed explanations.

### *dimension_formula.g*
This file contains the GAP program used in section 5.4.1 of the paper to compute the dimensions of the third and fourth DY cohomology groups for restricted quantum sl(2) at a 4th root of unity, thanks to the dimension formula in Corollary 4.9.

### *basis_H3.g*
This file contains the GAP program used in section 5.4.2 of the paper to compute a basis of the third DY cohomology group for restricted quantum sl(2) at a 4th root of unity.
