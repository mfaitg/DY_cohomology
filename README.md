# DY_cohomology
This repository contains the GAP programs used in section 5.4 of the paper [arXiv:2202.12287](https://arxiv.org/abs/2202.12287) to compute the Davydov-Yetter (DY) cohomology of the category of finite-dimensional modules over the restricted quantum group of sl(2) at a 4th root of unity.

We encode all the objects into linear-algebraic data and then the problems are reduced to solve (quite big) systems of linear equations. One could have used any software implementing linear algebra to do this; we chose GAP4 because of our familiarity with it.

The programs are extensively commented with detailed explanations.

### *dimension_formula.g*
This file contains the GAP program used in section 5.4.1 of the paper to compute the dimensions of the third and fourth DY cohomology groups for restricted quantum sl(2) at a 4th root of unity, thanks to the dimension formula in Corollary 4.9.

### *basis_H3.g*
This file contains the GAP program used in section 5.4.2 of the paper to compute a basis of the third DY cohomology group for restricted quantum sl(2) at a 4th root of unity.
