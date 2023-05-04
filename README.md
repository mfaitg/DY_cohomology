# DY_cohomology

## Introduction
This repository contains the [GAP](https://www.gap-system.org/) programs used in section 5.4 of the paper "Davydov-Yetter cohomology and relative homological algebra" by M. Faitg, A.M. Gainutdinov and C. Schweigert, [arXiv:2202.12287](https://arxiv.org/abs/2202.12287). It also contains the GAP program used in section 6.4 of the new version of the paper, which is not yet on the arXiv.

We wrote these programs in order to compute the 3rd and 4th Davydov-Yetter (DY) cohomology groups of the category of finite-dimensional modules over the restricted quantum group of sl(2) at a 3rd and 4th root of unity, which is computationally too difficult to be done by hand.

The idea is just to encode all the objects into linear-algebraic data so that the problems are reduced to solve (quite big) systems of linear equations. Note that one could have used any software implementing linear algebra to do this; we chose GAP4 because of our familiarity with it.

The programs are extensively commented with detailed explanations.

## Description of the files
If you want to read or modify the files you can open them with any text editor, like Notepad on Windows.

### *dimension_formula_DY_barUi.g*
This file contains the GAP program used in section 5.4.1 of the paper to compute the dimensions of the 3rd and 4th DY cohomology groups for restricted quantum sl(2) at a 4th root of unity, thanks to the dimension formula in Corollary 4.10.

### *basis_H3_DY_barUi.g*
This file contains the GAP program used in section 5.4.2 of the paper to compute a basis of the third DY cohomology group for restricted quantum sl(2) at a 4th root of unity. We know that this cohomology group is a 3-dimensional vector space thanks to the previous program.

### *dimension_formula_DY_uj.g*
This file contains the GAP program used in section 6.4 of the new version of the paper to compute the dimensions of the 3rd DY cohomology group for small quantum sl(2) at a 3rd root of unity, thanks to the dimension formula in Corollary 4.10.

## How to run these files
[GAP4](https://www.gap-system.org/Download/) must be installed on your computer. To run the files, use the Read command in the GAP4 terminal:

> Read("path/fileName.g");

For instance if you are on a Windows computer and if the file is in Desktop, this will be

> Read("C:/Users/UserName/Desktop/fileName.g");

where you have to replace UserName by your actual user name and fileName by the actual file name.

Alternatively, you can simply copy/paste the source code of the program in the GAP4 terminal.
