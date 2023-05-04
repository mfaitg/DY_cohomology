###########################################
## INTRODUCTION, EXPLANATIONS, NOTATIONS ##
###########################################

## WHAT IS THIS FILE? ##
# This file is attached to the paper "Davydov-Yetter cohomology and relative homological algebra", by M. Faitg, A.M. Gainutdinov and C. Schweigert.
# It contains the GAP program to compute the dimension of H3_DY for the category of finite-diemnsional modules over the small quantum group of sl(2) at the third root of unity thanks to 
# the dimension formula of Corollary 4.10 of the paper. More precisely we use the reformulation of this formula explained in section 6.4 of the paper.


## OUTPUT OF THIS PROGRAM ##
# When you run it in the GAP terminal, the program prints the dimension of the 3rd DY cohomology group of the category of finite-dimensional modules
# over the small quantum group of sl(2) at the third root of unity.


## ENCODING MODULES ## 
# The second tensor power of the small quantum group of sl(2) at a third root of unity j is generated as an algebra by E1, F1, K1, E2, F2, K2 
# such that each group of the variable with the same index satisfy the commutation relations of u_j(sl2) and each groups of variables commutes with the other. 
# Let M be a (finite-dimensional) module over this algebra. We encode M by a collection of matrices denoted by E1_M, F1_M, K1_M, E2_M, F2_M, K2_M, 
# which are the representations of these generators in a specified basis of M.

## NOTATIONS ##
# Recall that the relatively projective cover of the trivial module C over the second tensor power of u_j(sl2) is the block Q1 (see Prop. 5.4 of the paper). It has dimension 18. 
# The quotient of Q1 by the copy of C in its socle will be denoted by T (see the paper, where we note that the dual of the kernel of the projection Q1 --> C is isomorphic to T).
# The tensor product of T with T will be denoted by TT, the tensor product of T with Q1 will be denoted by TQ1 and so on.


## WHAT WE HAVE TO COMPUTE ##
# According to the dimension formula of the paper:
# dim(H^3_DY) = dim Inv(TTT) - dim Inv(TTQ1) + dim Inv(TT)
# So we have to:
# - define the representation matrices for Q1 and T,
# - compute the representations matrices for the tensor products TT, TTT, TTQ1
# - find the subspaces of invariant elements in these tensor products.


## TENSOR PRODUCT OF MODULES ##
# Let M and N be two modules over the second tensor power of u_j(sl2). They are encoded by collections of matrices E1_M, F1_M, K1_M, E2_M, F2_M, K2_M and E1_N, F1_N, K1_N, E2_N, F2_N, K2_N 
# respectively. Then the actions on the tensor product MN of M and N is defined by the coproduct :
# E1_MN = KroneckerProduct(I_M, E1_N) + KroneckerProduct(E1_M, K1_N), where I_M is the identity on M.
# F1_MN = KroneckerProduct(F1_M, I_N) + KroneckerProduct(K1_M^(-1), F1_N)
# K1_MN = KroneckerProduct(K1_M, K1_N)
# E2_MN = KroneckerProduct(I_M, E2_N) + KroneckerProduct(E2_M, K2_N), where I_N is the identity on N.
# F2_MN = KroneckerProduct(F2_M, I_N) + KroneckerProduct(K2_M^(-1), F2_N)
# K2_MN = KroneckerProduct(K2_M, K2_N)


## INVARIANT SUBSPACE ##
# Let M be a module over a Hopf algebra A. The invariant subspace Inv(M) is the subspace of elements on which the A-action is trivial: a.m = eps(a)m, where eps is the counit of A.
# Now consider a module M over the second tensor power of u_j(sl2), encoded by a collection of matrices E1_M, F1_M, K1_M, E2_M, F2_M, K2_M. 
# An element m in M is invariant if and only if E1_M.m = 0, F1_M.m = 0, K1_M.m = m, E2_M.m = 0, F2_M.m = 0, K2_M.m = m. Hence in order to obtain Inv(M), 
# we ask the software to compute the subspaces ker(E1_M), ker(F1_M), ker(K1_M - I_M), ker(E2_M), ker(F2_M), ker(K2_M - I_M) (where I_M is the identity on M)
# and then we ask it to compute the intersection of these subspaces.


#############
## PROGRAM ##
#############

j := E(3);; #3rd root of unity

## REPRESENTATIONS MATRICES FOR Q1 ##
# Dimension 18 #
# Recall that Q1 denotes the relatively projective cover of the trivial module C.

# Basis for Q1
vv_1 :=   [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0x0_1 := [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1x0_1 := [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0x1_1 := [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1x1_1 := [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0v_1 :=  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1v_1 :=  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0v_2 :=  [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1v_2 :=  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
vx0_1 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0];;
vx1_1 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];;
vx0_2 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];;
vx1_2 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];;
vv_2 :=   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];;
x0x0_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];;
x1x0_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0];;
x0x1_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];;
x1x1_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];;

zero18 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;

# Representations matrices :
I_Q1 := IdentityMat(18, CF(3)); # CF(3) is the cyclotomic field of the third root of unity over the Rationals.

E1_Q1 := TransposedMat([ x1v_1, vx0_1, x0x0_1 + x0x0_2, vx1_1, x0x1_1 + x0x1_2, zero18, x0v_1, vv_2, x0v_2, zero18, zero18, x1x0_2, x1x1_2, zero18, zero18, x0x0_2, zero18, x0x1_2 ]);;
F1_Q1 := TransposedMat([ x0v_2, x1x0_1, vx0_2, x1x1_1, vx1_2, x1v_1, vv_2, x1v_2, zero18, x0x0_2, x0x1_2, zero18, zero18, zero18, x1x0_2, zero18, x1x1_2, zero18 ]);;
K1_Q1 := TransposedMat([ vv_1, j*x0x0_1, j^(-1)*x1x0_1, j*x0x1_1, j^(-1)*x1x1_1, j*x0v_1, j^(-1)*x1v_1, j*x0v_2, j^(-1)*x1v_2, vx0_1, vx1_1, vx0_2, vx1_2, vv_2, j*x0x0_2, j^(-1)*x1x0_2, j*x0x1_2, j^(-1)*x1x1_2 ]);;

E2_Q1 := TransposedMat([ -j*vx1_1, -x0v_1, -x1v_1, -j^(-1)*x0x0_1 - j^(-1)*x0x0_2, -j^(-1)*x1x0_1 - j^(-1)*x1x0_2, zero18, zero18, -j*x0x1_2, -j*x1x1_2, zero18, -j^(-1)*vx0_1, -vv_2, -j^(-1)*vx0_2, zero18, zero18, zero18, -j^(-1)*x0x0_2, -j^(-1)*x1x0_2 ]);;
F2_Q1 := TransposedMat([ -vx0_2, -j*x0x1_1, -j*x1x1_1, -j^(-1)*x0v_2, -j^(-1)*x1v_2, -x0x0_2, -x1x0_2, zero18, zero18, -j*vx1_1, -j^(-1)*vv_2, -j*vx1_2, zero18, zero18, -j*x0x1_2, -j*x1x1_2, zero18, zero18 ]);;
K2_Q1 := TransposedMat([ vv_1, j*x0x0_1, j*x1x0_1, j^(-1)*x0x1_1, j^(-1)*x1x1_1, x0v_1, x1v_1, x0v_2, x1v_2, j*vx0_1, j^(-1)*vx1_1, j*vx0_2, j^(-1)*vx1_2, vv_2, j*x0x0_2, j*x1x0_2, j^(-1)*x0x1_2, j^(-1)*x1x1_2 ]);;

#N.B. : to find these matrices we used the diagrams given in section 4.4.2 of arXiv:hep-th/0504093 (see the remarks in the proof of Prop. 6.8 of our paper).

#We delete the following variables to re-use their names below with other values:
Unbind(vv_1); Unbind(x0x0_1); Unbind(x1x0_1); Unbind(x0x1_1); Unbind(x1x1_1); Unbind(x0v_1); Unbind(x1v_1); Unbind(x0v_2);
Unbind(x1v_2); Unbind(vx0_1); Unbind(vx1_1); Unbind(vx0_2); Unbind(vx1_2); Unbind(vv_2); Unbind(x0x0_2); Unbind(x1x0_2);
Unbind(x0x1_2); Unbind(x1x1_2);


## REPRESENTATIONS MATRICES FOR T = Q1/C ##
# Dimension 17 #
# T has the same basis than Q1, except that we kill vv_2, which generates the submodule isomorphic to C.

vv_1 :=   [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0x0_1 := [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1x0_1 := [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0x1_1 := [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1x1_1 := [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0v_1 :=  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1v_1 :=  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x0v_2 :=  [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];;
x1v_2 :=  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0];;
vx0_1 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0];;
vx1_1 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];;
vx0_2 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0];;
vx1_2 :=  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0];;
x0x0_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];;
x1x0_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0];;
x0x1_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0];;
x1x1_2 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];;

zero17 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;

# Representations matrices :
I_T := IdentityMat(17, CF(3));;

E1_T := TransposedMat([ x1v_1, vx0_1, x0x0_1 + x0x0_2, vx1_1, x0x1_1 + x0x1_2, zero17, x0v_1, zero17, x0v_2, zero17, zero17, x1x0_2, x1x1_2, zero17, x0x0_2, zero17, x0x1_2 ]);;
F1_T := TransposedMat([ x0v_2, x1x0_1, vx0_2, x1x1_1, vx1_2, x1v_1, zero17, x1v_2, zero17, x0x0_2, x0x1_2, zero17, zero17, x1x0_2, zero17, x1x1_2, zero17 ]);;
K1_T := TransposedMat([ vv_1, j*x0x0_1, j^(-1)*x1x0_1, j*x0x1_1, j^(-1)*x1x1_1, j*x0v_1, j^(-1)*x1v_1, j*x0v_2, j^(-1)*x1v_2, vx0_1, vx1_1, vx0_2, vx1_2, j*x0x0_2, j^(-1)*x1x0_2, j*x0x1_2, j^(-1)*x1x1_2 ]);;

E2_T := TransposedMat([ -j*vx1_1, -x0v_1, -x1v_1, -j^(-1)*x0x0_1 - j^(-1)*x0x0_2, -j^(-1)*x1x0_1 - j^(-1)*x1x0_2, zero17, zero17, -j*x0x1_2, -j*x1x1_2, zero17, -j^(-1)*vx0_1, zero17, -j^(-1)*vx0_2, zero17, zero17, -j^(-1)*x0x0_2, -j^(-1)*x1x0_2 ]);;
F2_T := TransposedMat([ -vx0_2, -j*x0x1_1, -j*x1x1_1, -j^(-1)*x0v_2, -j^(-1)*x1v_2, -x0x0_2, -x1x0_2, zero17, zero17, -j*vx1_1, zero17, -j*vx1_2, zero17, -j*x0x1_2, -j*x1x1_2, zero17, zero17 ]);;
K2_T := TransposedMat([ vv_1, j*x0x0_1, j*x1x0_1, j^(-1)*x0x1_1, j^(-1)*x1x1_1, x0v_1, x1v_1, x0v_2, x1v_2, j*vx0_1, j^(-1)*vx1_1, j*vx0_2, j^(-1)*vx1_2, j*x0x0_2, j*x1x0_2, j^(-1)*x0x1_2, j^(-1)*x1x1_2 ]);;

# These matrices are derived from the matrices for Q1.


## REPRESENTATION MATRICES FOR TT = T \otimes T ##
# Dimension 17^2 = 289 #
Print("Computing the representations matrices for TT...\n\n");

I_TT := IdentityMat(17^2);;
E1_TT := KroneckerProduct(I_T,E1_T) + KroneckerProduct(E1_T,K1_T);;
F1_TT := KroneckerProduct(F1_T,I_T) + KroneckerProduct(K1_T^(-1),F1_T);;
K1_TT := KroneckerProduct(K1_T,K1_T);;
E2_TT := KroneckerProduct(I_T,E2_T) + KroneckerProduct(E2_T,K2_T);;
F2_TT := KroneckerProduct(F2_T,I_T) + KroneckerProduct(K2_T^(-1),F2_T);;
K2_TT := KroneckerProduct(K2_T,K2_T);;


## INVARIANT SUBSPACE IN TT ##
Print("Computing the invariant subspace in TT...\n\n");

# Remark: the command NullspaceMat returns a basis of the LEFT kernel of a matrix. Since we want to compute RIGHT kernels,
# we have to apply NullspaceMat to the transpose of a matrix.

InvE1_TT := Subspace(CF(3)^(17^2), NullspaceMat(TransposedMat(E1_TT)),"basis");; #Invariant subspace under the action of E1
InvF1_TT := Subspace(CF(3)^(17^2), NullspaceMat(TransposedMat(F1_TT)),"basis");; #Invariant subspace under the action of F1
InvK1_TT := Subspace(CF(3)^(17^2), NullspaceMat(TransposedMat(K1_TT) - IdentityMat(17^2)),"basis");; #Invariant subspace under the action of K1
InvE2_TT := Subspace(CF(3)^(17^2), NullspaceMat(TransposedMat(E2_TT)),"basis");; #Invariant subspace under the action of E1
InvF2_TT := Subspace(CF(3)^(17^2), NullspaceMat(TransposedMat(F2_TT)),"basis");; #Invariant subspace under the action of F1
InvK2_TT := Subspace(CF(3)^(17^2), NullspaceMat(TransposedMat(K2_TT) - IdentityMat(17^2)),"basis");; #Invariant subspace under the action of K1

# The subspace of invariant elements is the intersection of these subspaces:
Inv_TT := Intersection(InvE1_TT, InvF1_TT, InvK1_TT, InvE2_TT, InvF2_TT, InvK2_TT);;

# We no longer need the subspaces InvE1_TT, ..., InvK2_TT thus we can delete them to reduce the memory space:
Unbind(InvE1_TT); Unbind(InvF1_TT); Unbind(InvK1_TT); Unbind(InvE2_TT); Unbind(InvF2_TT); Unbind(InvK2_TT);


## INVARIANT SUBSPACE IN TTT ##
# Dimension 17^3 = 4913 #
# The dimension is big, so we organize the computation in such a way that it requires less memory space.

Print("Computing the invariant subspace in TTT...\n\n");

E1_TTT := KroneckerProduct(I_TT,E1_T) + KroneckerProduct(E1_TT,K1_T);;
Inv_TTT := Subspace(CF(3)^(17^3), NullspaceMat(TransposedMat(E1_TTT)),"basis");;
Unbind(E1_TTT);

F1_TTT := KroneckerProduct(F1_TT,I_T) + KroneckerProduct(K1_TT^(-1),F1_T);;
Inv_TTT := Intersection(Inv_TTT, Subspace(CF(3)^(17^3), NullspaceMat(TransposedMat(F1_TTT)),"basis"));;
Unbind(F1_TTT);

K1_TTT := KroneckerProduct(K1_TT,K1_T);;
Inv_TTT := Intersection(Inv_TTT, Subspace(CF(3)^(17^3), NullspaceMat(TransposedMat(K1_TTT) - IdentityMat(17^3)),"basis"));;
Unbind(K1_TTT);

E2_TTT := KroneckerProduct(I_TT,E2_T) + KroneckerProduct(E2_TT,K2_T);;
Inv_TTT := Intersection(Inv_TTT, Subspace(CF(3)^(17^3), NullspaceMat(TransposedMat(E2_TTT)),"basis"));;
Unbind(E2_TTT);

F2_TTT := KroneckerProduct(F2_TT,I_T) + KroneckerProduct(K2_TT^(-1),F2_T);;
Inv_TTT := Intersection(Inv_TTT, Subspace(CF(3)^(17^3), NullspaceMat(TransposedMat(F2_TTT)),"basis"));;
Unbind(F2_TTT);

K2_TTT := KroneckerProduct(K2_TT,K2_T);;
Inv_TTT := Intersection(Inv_TTT, Subspace(CF(3)^(17^3), NullspaceMat(TransposedMat(K2_TTT) - IdentityMat(17^3)),"basis"));;
Unbind(K2_TTT);


## INVARIANT SUBSPACE IN TTQ1 ##
# Dimension 17^2*18 = 5202 #
# The dimension is big, so we organize the computation in such a way that it requires less memory space.

Print("Computing the invariant subspace in TTQ1...\n\n");

E1_TTQ1 := KroneckerProduct(I_TT,E1_Q1) + KroneckerProduct(E1_TT,K1_Q1);;
Inv_TTQ1 := Subspace(CF(3)^(17^2*18), NullspaceMat(TransposedMat(E1_TTQ1)),"basis");;
Unbind(E1_TTQ1);

F1_TTQ1 := KroneckerProduct(F1_TT,I_Q1) + KroneckerProduct(K1_TT^(-1),F1_Q1);;
Inv_TTQ1 := Intersection(Inv_TTQ1, Subspace(CF(3)^(17^2*18), NullspaceMat(TransposedMat(F1_TTQ1)),"basis"));;
Unbind(F1_TTQ1);

K1_TTQ1 := KroneckerProduct(K1_TT,K1_Q1);;
Inv_TTQ1 := Intersection(Inv_TTQ1, Subspace(CF(3)^(17^2*18), NullspaceMat(TransposedMat(K1_TTQ1) - IdentityMat(17^2*18)),"basis"));;
Unbind(K1_TTQ1);

E2_TTQ1 := KroneckerProduct(I_TT,E2_Q1) + KroneckerProduct(E2_TT,K2_Q1);;
Inv_TTQ1 := Intersection(Inv_TTQ1, Subspace(CF(3)^(17^2*18), NullspaceMat(TransposedMat(E2_TTQ1)),"basis"));;
Unbind(E2_TTT);

F2_TTQ1 := KroneckerProduct(F2_TT,I_Q1) + KroneckerProduct(K2_TT^(-1),F2_Q1);;
Inv_TTQ1 := Intersection(Inv_TTQ1, Subspace(CF(3)^(17^2*18), NullspaceMat(TransposedMat(F2_TTQ1)),"basis"));;
Unbind(F2_TTQ1);

K2_TTQ1 := KroneckerProduct(K2_TT,K2_Q1);;
Inv_TTQ1 := Intersection(Inv_TTQ1, Subspace(CF(3)^(17^2*18), NullspaceMat(TransposedMat(K2_TTQ1) - IdentityMat(17^2*18)),"basis"));;
Unbind(K2_TTQ1);


#########################
## PRINTING THE RESULT ##
#########################
Print("RESULTS:\n");
Print("dim(H3_DY) = dim(Inv(TTT)) - dim(Inv(TTQ1)) + dim(Inv(TT)) = ", Dimension(Inv_TTT), " - ", Dimension(Inv_TTQ1), 
      " + ", Dimension(Inv_TT), " = ", Dimension(Inv_TTT) - Dimension(Inv_TTQ1) + Dimension(Inv_TT), "\n");
