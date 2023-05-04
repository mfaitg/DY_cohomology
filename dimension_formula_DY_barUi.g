###########################################
## INTRODUCTION, EXPLANATIONS, NOTATIONS ##
###########################################

## WHAT IS THIS FILE? ##
# This is the GAP program used in section 5.4.1 of arXiv2202.12287 to compute the dimension of H3_DY and H4_DY for the 
# category of finite-diemnsional modules over the restricted quantum group of sl(2) at the fourth root of unity thanks to 
# the dimension formula of Corollary 4.10. More precisely we use the reformulation of this formula in terms of invariant
# subspaces as given in equation (83) of the paper.


## OUTPUT OF THIS PROGRAM ##
# When you run it in the GAP terminal, the program prints the dimension of the 3rd and of the 4th DY cohomology groups 
# of the category of finite-dimensional modules over the restricted quantum group of sl(2) at the fourth root of unity.


## ENCODING MODULES ## 
# Recall that the Drinfeld double of the restricted quantum group of sl(2) at a fourth root of unity is generated as an
# algebra by E, F, K, a, b, c, d. Let M be a (finite-dimensional) module over this Drinfeld double. We encode M by a 
# collection of matrices denoted by E_M, F_M, K_M, a_M, b_M, c_M, d_M, which are the representations of these generators
# in a specified basis of M.
# Note that due to the q-determinant relation we have a = d^3 + i*b*c*d^3, where i^2 = -1. We use this remark below to
# define the representations matrices for a.


## NOTATIONS ##
# In the sequel, the relatively projective cover of the trivial module C over the Drinfeld double of restricted quantum sl(2)
# will be denoted by R. It is of dimension 8. The quotient of R by the copy of C in its socle will be denoted by Q (see the
# paper, where we note that the dual of the kernel of the projection R --> C is isomorphic to Q).
# The tensor product of Q with Q will be denoted by QQ, the tensor product of Q with R will be denoted by QR and so on.


## WHAT WE HAVE TO COMPUTE ##
# According to the dimension formula in equation (83) of the paper:
# dim(H^3_DY) = dim Inv(QQQ) - dim Inv(QQR) + dim Inv(QQ)
# dim(H^4_DY) = dim Inv(QQQQ) - dim Inv(QQQR) + dim Inv(QQQ)
# So we have to:
# - define the representation matrices for R and Q,
# - compute the representations matrices for the tensor products QQ, QQQ, QQR, QQQQ, QQQR,
# - find the subspaces of invariant elements in these tensor products.


## TENSOR PRODUCT OF MODULES ##
# Let M and N be two modules over the Drinfeld center. They are encoded by collections of matrices E_M, F_M, K_M, a_M, 
# b_M, c_M, d_M and E_N, F_N, K_N, a_N, b_N, c_N, d_N respectively. Then the actions on the tensor product MN of M and N
# is defined by the coproduct :
# E_MN = KroneckerProduct(I_M, E_N) + KroneckerProduct(E_M, K_N), where I_M is the identity on M.
# F_MN = KroneckerProduct(F_M, I_N) + KroneckerProduct(K_M^(-1), F_N)
# K_MN = KroneckerProduct(K_M, K_N)
# a_MN = KroneckerProduct(a_M, a_N) + KroneckerProduct(b_M, b_N)
# b_MN = KroneckerProduct(a_M, b_N) + KroneckerProduct(b_M, d_N)
# c_MN = KroneckerProduct(c_M, a_N) + KroneckerProduct(d_M, c_N)
# d_MN = KroneckerProduct(c_M, b_N) + KroneckerProduct(d_M, d_N)


## INVARIANT SUBSPACE ##
# Let M be a module over a Hopf algebra H. The invariant subspace Inv(M) is the subspace of elements on which the H-action
# is trivial: h.m = eps(h)m, where eps is the counit of H.
# Now consider a module M over the Drinfeld double of restricted quantum sl(2), encoded by a collection of matricesE_M, F_M,
# K_M, a_M, b_M, c_M, d_M. An element m in M is invariant if and only if
# E_M.m = 0, F_M.m = 0, K_M.m = m, b_M.m = 0, c_M.m = 0, d_M.m = m
# (the equation a_M.m = m is not required since a_M can be expressed as a_M = d_M^3 + i*b_M*c_M*d_M^3)
# Hence in order to obtain Inv(M), we ask the software to compute the subspaces 
# ker(E_M), ker(F_M), ker(K_M - I_M), ker(b_M), ker(c_M), ker(d_M - I_M)   (where I_M is the identity on M)
# and then we ask it to compute the intersection of these subspaces.



#############
## PROGRAM ##
#############

i := E(4);; #4th root of unity

## REPRESENTATIONS MATRICES FOR R ##
# Dimension 8 #
# Recall that R denotes the relatively projective cover of the trivial module C.

# Basis for R (notations are consistent with those in the paper)
phi0   := [1, 0, 0, 0, 0, 0, 0, 0];;
phi2   := [0, 1, 0, 0, 0, 0, 0, 0];;
bphi0  := [0, 0, 1, 0, 0, 0, 0, 0];;
bphi2  := [0, 0, 0, 1, 0, 0, 0, 0];;
cphi0  := [0, 0, 0, 0, 1, 0, 0, 0];;
cphi2  := [0, 0, 0, 0, 0, 1, 0, 0];;
bcphi0 := [0, 0, 0, 0, 0, 0, 1, 0];;
bcphi2 := [0, 0, 0, 0, 0, 0, 0, 1];;

zero8 := [0, 0, 0, 0, 0, 0, 0, 0];;

# Representations matrices :
I_R := IdentityMat(8, CF(4)); # CF(4) is the cyclotomic field of the fourth root of unity over the Rationals.
E_R := TransposedMat([ (1/2)*(cphi0 + cphi2), -(1/2)*(cphi0 + cphi2), (i/2)*(bcphi0 - bcphi2), (i/2)*(bcphi0 - bcphi2), zero8, zero8, zero8, zero8 ]);;
F_R := TransposedMat([ (i/2)*(bphi0 + bphi2), -(i/2)*(bphi0 + bphi2), zero8, zero8, (1/2)*(bcphi2 - bcphi0), (1/2)*(bcphi2 - bcphi0), zero8, zero8 ]);;
b_R := TransposedMat([ bphi0, bphi2, zero8, zero8, bcphi0, bcphi2, zero8, zero8 ]);;
c_R := TransposedMat([ cphi0, cphi2, bcphi0, bcphi2, zero8, zero8, zero8, zero8 ]);;
d_R := TransposedMat([ phi0, -phi2, -i*bphi0, i*bphi2, -i*cphi0, i*cphi2, -bcphi0, bcphi2 ]);;
K_R := d_R^2;; #because K acts as d^2 on this particular module, see the paper.
a_R := d_R^3 + i*b_R*c_R*d_R^3;;


## REPRESENTATIONS MATRICES FOR Q = R/C ##
# Dimension 7 #
# Q has the same basis than RC, except that we kill bcphi2, which generates the submodule isomorphic to C.

phi0   := [1, 0, 0, 0, 0, 0, 0];;
phi2   := [0, 1, 0, 0, 0, 0, 0];;
bphi0  := [0, 0, 1, 0, 0, 0, 0];;
bphi2  := [0, 0, 0, 1, 0, 0, 0];;
cphi0  := [0, 0, 0, 0, 1, 0, 0];;
cphi2  := [0, 0, 0, 0, 0, 1, 0];;
bcphi0 := [0, 0, 0, 0, 0, 0, 1];;

zero7 := [0, 0, 0, 0, 0, 0, 0];;

# Representations matrices :
I_Q := IdentityMat(7, CF(4)); # CF(4) is the cyclotomic field of the fourth root of unity over the Rationals.
E_Q := TransposedMat([ (1/2)*(cphi0 + cphi2), -(1/2)*(cphi0 + cphi2), (i/2)*bcphi0, (i/2)*bcphi0, zero7, zero7, zero7]);;
F_Q := TransposedMat([ (i/2)*(bphi0 + bphi2), -(i/2)*(bphi0 + bphi2), zero7, zero7, -(1/2)*bcphi0, -(1/2)*bcphi0, zero7]);;
b_Q := TransposedMat([ bphi0, bphi2, zero7, zero7, bcphi0, zero7, zero7]);;
c_Q := TransposedMat([ cphi0, cphi2, bcphi0, zero7, zero7, zero7, zero7]);;
d_Q := TransposedMat([ phi0, -phi2, -i*bphi0, i*bphi2, -i*cphi0, i*cphi2, -bcphi0]);;
K_Q := d_Q^2;; #because K acts as d^2 on this particular module, see the paper.
a_Q := d_Q^3 + i*b_Q*c_Q*d_Q^3;;


## REPRESENTATION MATRICES FOR QQ = Q \otimes Q ##
# Dimension 49 #
Print("Computing the representations matrices for QQ...\n\n");

I_QQ := IdentityMat(49);;
E_QQ := KroneckerProduct(I_Q,E_Q) + KroneckerProduct(E_Q,K_Q);;
F_QQ := KroneckerProduct(F_Q,I_Q) + KroneckerProduct(K_Q^(-1),F_Q);;
K_QQ := KroneckerProduct(K_Q,K_Q);;
a_QQ := KroneckerProduct(a_Q,a_Q) + KroneckerProduct(b_Q,c_Q);;
b_QQ := KroneckerProduct(a_Q, b_Q) + KroneckerProduct(b_Q,d_Q);;
c_QQ := KroneckerProduct(c_Q,a_Q) + KroneckerProduct(d_Q, c_Q);;
d_QQ := KroneckerProduct(c_Q,b_Q) + KroneckerProduct(d_Q,d_Q);;


## INVARIANT SUBSPACE IN QQ ##
Print("Computing the invariant subspace in QQ...\n\n");

# Remark: the command NullspaceMat returns a basis of the LEFT kernel of a matrix. Since we want to compute RIGHT kernels,
# we have to apply NullspaceMat to the transpose of a matrix.

InvE_QQ := Subspace(CF(4)^49, NullspaceMat(TransposedMat(E_QQ)),"basis");; #Invariant subspace under the action of E
InvF_QQ := Subspace(CF(4)^49, NullspaceMat(TransposedMat(F_QQ)),"basis");; #Invariant subspace under the action of F
InvK_QQ := Subspace(CF(4)^49, NullspaceMat(TransposedMat(K_QQ) - IdentityMat(49)),"basis");; #Invariant subspace under the action of K
Invb_QQ := Subspace(CF(4)^49, NullspaceMat(TransposedMat(b_QQ)),"basis");; #Invariant subspace under the action of b
Invc_QQ := Subspace(CF(4)^49, NullspaceMat(TransposedMat(c_QQ)),"basis");; #Invariant subspace under the action of c
Invd_QQ := Subspace(CF(4)^49, NullspaceMat(TransposedMat(d_QQ) - IdentityMat(49)),"basis");; #Invariant subspace under the action of d

# The subspace of invariant elements is the intersection of these subspaces, since E,F,K,b,c,d generate the Drinfeld double:
Inv_QQ := Intersection(InvE_QQ, InvF_QQ, InvK_QQ, Invb_QQ, Invc_QQ, Invd_QQ);;

# We no longer need the subspaces InvE_QQ, ..., Invd_QQ thus we can delete them to reduce the memory space:
Unbind(InvE_QQ); Unbind(InvF_QQ); Unbind(InvK_QQ); Unbind(Invb_QQ); Unbind(Invc_QQ); Unbind(Invd_QQ);


## The computation is completely similar for the other modules. ##


## REPRESENTATIONS MATRICES FOR QQQ = Q \otimes Q \otimes Q ##
# Dimension 343 #
Print("Computing the representations matrices for QQQ...\n\n");

I_QQQ := IdentityMat(343);;
E_QQQ := KroneckerProduct(I_QQ,E_Q) + KroneckerProduct(E_QQ,K_Q);;
F_QQQ := KroneckerProduct(F_QQ,I_Q) + KroneckerProduct(K_QQ^(-1),F_Q);;
K_QQQ := KroneckerProduct(K_QQ,K_Q);;
a_QQQ := KroneckerProduct(a_QQ,a_Q) + KroneckerProduct(b_QQ,c_Q);;
b_QQQ := KroneckerProduct(a_QQ, b_Q) + KroneckerProduct(b_QQ,d_Q);;
c_QQQ := KroneckerProduct(c_QQ,a_Q) + KroneckerProduct(d_QQ,c_Q);;
d_QQQ := KroneckerProduct(c_QQ,b_Q) + KroneckerProduct(d_QQ,d_Q);;


## INVARIANT SUBSPACE IN QQQ ##
Print("Computing the invariant subspace in QQQ...\n\n");

InvE_QQQ := Subspace(CF(4)^343, NullspaceMat(TransposedMat(E_QQQ)),"basis");;
InvF_QQQ := Subspace(CF(4)^343, NullspaceMat(TransposedMat(F_QQQ)),"basis");;
InvK_QQQ := Subspace(CF(4)^343, NullspaceMat(TransposedMat(K_QQQ) - IdentityMat(343)),"basis");;
Invb_QQQ := Subspace(CF(4)^343, NullspaceMat(TransposedMat(b_QQQ)),"basis");;
Invc_QQQ := Subspace(CF(4)^343, NullspaceMat(TransposedMat(c_QQQ)),"basis");;
Invd_QQQ := Subspace(CF(4)^343, NullspaceMat(TransposedMat(d_QQQ) - IdentityMat(343)),"basis");;

Inv_QQQ := Intersection(InvE_QQQ, InvF_QQQ, InvK_QQQ, Invb_QQQ, Invc_QQQ, Invd_QQQ);;
Unbind(InvE_QQQ); Unbind(InvF_QQQ); Unbind(InvK_QQQ); Unbind(Invb_QQQ); Unbind(Invc_QQQ); Unbind(Invd_QQQ);


## REPRESENTATIONS MATRICES FOR QQR = Q \otimes Q \otimes R ##
# Dimension 392 #
Print("Computing the representations matrices for QQR...\n\n");

E_QQR := KroneckerProduct(I_QQ,E_R) + KroneckerProduct(E_QQ,K_R);;
F_QQR := KroneckerProduct(F_QQ,I_R) + KroneckerProduct(K_QQ^(-1),F_R);;
K_QQR := KroneckerProduct(K_QQ,K_R);;
b_QQR := KroneckerProduct(a_QQ, b_R) + KroneckerProduct(b_QQ,d_R);;
c_QQR := KroneckerProduct(c_QQ,a_R) + KroneckerProduct(d_QQ,c_R);;
d_QQR := KroneckerProduct(c_QQ,b_R) + KroneckerProduct(d_QQ,d_R);;


## INVARIANT SUBSPACE IN QQR ##
Print("Computing the invariant subspace in QQR...\n\n");

InvE_QQR := Subspace(CF(4)^392, NullspaceMat(TransposedMat(E_QQR)),"basis");;
InvF_QQR := Subspace(CF(4)^392, NullspaceMat(TransposedMat(F_QQR)),"basis");;
InvK_QQR := Subspace(CF(4)^392, NullspaceMat(TransposedMat(K_QQR) - IdentityMat(392)),"basis");;
Invb_QQR := Subspace(CF(4)^392, NullspaceMat(TransposedMat(b_QQR)),"basis");;
Invc_QQR := Subspace(CF(4)^392, NullspaceMat(TransposedMat(c_QQR)),"basis");;
Invd_QQR := Subspace(CF(4)^392, NullspaceMat(TransposedMat(d_QQR) - IdentityMat(392)),"basis");;

Inv_QQR := Intersection(InvE_QQR, InvF_QQR, InvK_QQR, Invb_QQR, Invc_QQR, Invd_QQR);;
Unbind(InvE_QQR); Unbind(InvF_QQR); Unbind(InvK_QQR); Unbind(Invb_QQR); Unbind(Invc_QQR); Unbind(Invd_QQR);


## REPRESENTATIONS MATRICES FOR QQQQ = Q \otimes Q \otimes Q \otimes Q ##
# Dimension 2401 #
Print("Computing the representations matrices for QQQQ...\n\n");

I_QQQQ := IdentityMat(2401);;
E_QQQQ := KroneckerProduct(I_QQQ,E_Q) + KroneckerProduct(E_QQQ,K_Q);;
F_QQQQ := KroneckerProduct(F_QQQ,I_Q) + KroneckerProduct(K_QQQ^(-1),F_Q);;
K_QQQQ := KroneckerProduct(K_QQQ,K_Q);;
a_QQQQ := KroneckerProduct(a_QQQ,a_Q) + KroneckerProduct(b_QQQ,c_Q);;
b_QQQQ := KroneckerProduct(a_QQQ, b_Q) + KroneckerProduct(b_QQQ,d_Q);;
c_QQQQ := KroneckerProduct(c_QQQ,a_Q) + KroneckerProduct(d_QQQ,c_Q);;
d_QQQQ := KroneckerProduct(c_QQQ,b_Q) + KroneckerProduct(d_QQQ,d_Q);;


## INVARIANT SUBSPACE IN QQQQ ##
Print("Computing the invariant subspace in QQQQ...\n\n");

InvE_QQQQ := Subspace(CF(4)^2401, NullspaceMat(TransposedMat(E_QQQQ)),"basis");;
InvF_QQQQ := Subspace(CF(4)^2401, NullspaceMat(TransposedMat(F_QQQQ)),"basis");;
InvK_QQQQ := Subspace(CF(4)^2401, NullspaceMat(TransposedMat(K_QQQQ) - IdentityMat(2401)),"basis");;
Invb_QQQQ := Subspace(CF(4)^2401, NullspaceMat(TransposedMat(b_QQQQ)),"basis");;
Invc_QQQQ := Subspace(CF(4)^2401, NullspaceMat(TransposedMat(c_QQQQ)),"basis");;
Invd_QQQQ := Subspace(CF(4)^2401, NullspaceMat(TransposedMat(d_QQQQ) - IdentityMat(2401)),"basis");;

Inv_QQQQ := Intersection(InvE_QQQQ, InvF_QQQQ, InvK_QQQQ, Invb_QQQQ, Invc_QQQQ, Invd_QQQQ);;
Unbind(InvE_QQQQ); Unbind(InvF_QQQQ); Unbind(InvK_QQQQ); Unbind(Invb_QQQQ); Unbind(Invc_QQQQ); Unbind(Invd_QQQQ);


## REPRESENTATIONS MATRICES FOR QQQR = Q \otimes Q \otimes Q \otimes R ##
# Dimension 2744 #
Print("Computing the representations matrices for QQQR...\n\n");

I_QQQR := IdentityMat(2744);;
E_QQQR := KroneckerProduct(I_QQQ,E_R) + KroneckerProduct(E_QQQ,K_R);;
F_QQQR := KroneckerProduct(F_QQQ,I_R) + KroneckerProduct(K_QQQ^(-1),F_R);;
K_QQQR := KroneckerProduct(K_QQQ,K_R);;
b_QQQR := KroneckerProduct(a_QQQ,b_R) + KroneckerProduct(b_QQQ,d_R);;
c_QQQR := KroneckerProduct(c_QQQ,a_R) + KroneckerProduct(d_QQQ,c_R);;
d_QQQR := KroneckerProduct(c_QQQ,b_R) + KroneckerProduct(d_QQQ,d_R);;


## INVARIANT SUBSPACE IN QQQR ##
Print("Computing the invariant subspace in QQQR...\n\n");

InvE_QQQR := Subspace(CF(4)^2744, NullspaceMat(TransposedMat(E_QQQR)),"basis");;
InvF_QQQR := Subspace(CF(4)^2744, NullspaceMat(TransposedMat(F_QQQR)),"basis");;
InvK_QQQR := Subspace(CF(4)^2744, NullspaceMat(TransposedMat(K_QQQR) - IdentityMat(2744)),"basis");;
Invb_QQQR := Subspace(CF(4)^2744, NullspaceMat(TransposedMat(b_QQQR)),"basis");;
Invc_QQQR := Subspace(CF(4)^2744, NullspaceMat(TransposedMat(c_QQQR)),"basis");;
Invd_QQQR := Subspace(CF(4)^2744, NullspaceMat(TransposedMat(d_QQQR) - IdentityMat(2744)),"basis");;

Inv_QQQR := Intersection(InvE_QQQR, InvF_QQQR, InvK_QQQR, Invb_QQQR, Invc_QQQR, Invd_QQQR);;
Unbind(InvE_QQQR); Unbind(InvF_QQQR); Unbind(InvK_QQQR); Unbind(Invb_QQQR); Unbind(Invc_QQQR); Unbind(Invd_QQQR);



#############
## RESULTS ##
#############
Print("RESULTS:\n");
Print("dim(H3_DY) = dim(Inv(QQQ)) - dim(Inv(QQR)) + dim(Inv(QQ)) = ", Dimension(Inv_QQQ), " - ", Dimension(Inv_QQR), 
      " + ", Dimension(Inv_QQ), " = ", Dimension(Inv_QQQ) - Dimension(Inv_QQR) + Dimension(Inv_QQ), "\n");
Print("dim(H4_DY) = dim(Inv(QQQQ)) - dim(Inv(QQQR)) + dim(Inv(QQQ)) = ", Dimension(Inv_QQQQ), " - ", Dimension(Inv_QQQR), 
      " + ", Dimension(Inv_QQQ), " = ", Dimension(Inv_QQQQ) - Dimension(Inv_QQQR) + Dimension(Inv_QQQ), "\n");
