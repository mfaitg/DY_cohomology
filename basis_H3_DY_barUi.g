###########################################
## INTRODUCTION, EXPLANATIONS, NOTATIONS ##
###########################################

## WHAT IS THIS FILE? ##
# This is the GAP program used in section 5.4.2 of arXiv2202.12287 to compute a basis of explicit cocycles in the third
# Davydov-Yetter (DY) cohomology group of the category of finite-dimensional modules over the restricted quantum group 
# of sl(2) at a fourth root of unity. We use the reformulation of the DY complex in terms of the centralizers of the iterated
# coproducts with the Cartier differential, see equations (52) and (53) in the paper.


## OUTPUT OF THIS PROGRAM ##
# When you run it in the GAP terminal, the program prints three explicit DY cocycles which are a basis of the third 
# cohomology group H3_DY. They are displayed in a readable form but since we are in a text environment there is no symbol
# for the tensor product. Thus we display the tensor product of two vectors u and v as [u,v].
# For instance the first cocycle is displayed as:
# (1)*[1, EK2, EK3] + (-1)*[K2, EK2, EK3] + (-1)*[E, K, EK3] + (1)*[E, K3, EK3] + (-1)*[E, EK, 1] + (1)*[E, EK, K2]
# The numbers between brackets are the coefficients of the tensors and for instance [1, EK2, EK3] means 1 tensor EK2 tensor EK3,
# EK2 and EK3 are elements of the monomial basis of restricted quantum sl(2) (meaning respectively EK^2 and EK^3).


## NOTATIONS ##
# In the comments We use LaTeX-type notations, like \bar U_i(sl2) for the restricted quantum group of sl(2) at the 4th root
# of unity i (i^2 = -1) or \Delta for the coproduct of \bar U_i(sl2). The definition of \bar U_i(sl2) is recalled in §5.4 
# of the paper.


## ORGANIZATION OF THE PROGRAM ##
# The program is divided into the following parts:
#   - Description of \bar U_i(sl2) with linear-algebraic data
#   - Computation of the cochain space in degree 2
#   - Computation of the cochain space in degree 3
#   - Construction of the degree 2 differential d2_DY and computation of its image
#   - Construction of the degree 3 differential d3_DY and computation of its kernel
#   - Computation of explicit cocycles in H3_DY
#   - Display of the results.
# Each part contains lots of detailed comments.



#############################################################
## DESCRIPTION OF \bar U_i(sl2) WITH LINEAR-ALGEBRAIC DATA ##
#############################################################
Print("Constructing the multiplication and comultiplication tensors...\n\n");

i := E(4); #Primitive 4th root of unity.
Qi := GaussianRationals; #Cyclotomic field containing i.

# We represent \bar U_i(sl2) as the 16-dimensional vector space over Qi: Qi^16.
# Let B be the canonical basis of Qi^16:
B := CanonicalBasis(Qi^16);

# We decide that the elements of B correspond to the monomial basis of \bar U_i(sl2) as follows:
# B[1] = 1,    B[2] = K,     B[3] = K^2,     B[4] = K^3,
# B[5] = E,    B[6] = EK,    B[7] = EK^2,    B[8] = EK^3,
# B[9] = F,    B[10] = FK,   B[11] = FK^2,   B[12] = FK^3,
# B[13] = EF,  B[14] = EFK,  B[15] = EFK^2,  B[16] = EFK^3

zero16 := [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];; #Zero element of \bar U_i(sl2)
I_16 := IdentityMat(16);; #Identity matrix of size 16 over Qi.

# The bialgebra structure of \bar U_i(sl2) will be represented by:
#   - A multiplication tensor m : \bar U_i(sl2) x \bar U_i(sl2) --> \bar U_i(sl2), i.e. a matrix of size 16x256.
#       In order to construct m, we first construct matrices mE, mF, mK : \bar U_i(sl2) --> \bar U_i(sl2) of size 16x16
#       which represent the left multiplication by E, F, K respectively.
#   - A comultiplication tensor com : \bar U_i(sl2) --> \bar U_i(sl2) x \bar U_i(sl2), i.e. a matrix of size 256x16.
# We do not need the counit and the antipode but one can easily represent them as matrices of size 1x16 and 16x16 respectively.


## MATRICES FOR THE LEFT MULTIPLICATIONS BY E, F, K ##
# IMPORTANT REMARK about matrices in GAP: a matrix is a list of ROW vectors; for instance [ [1, 2], [3, 4] ] represents 
# the matrix  1  2
#             3  4
# Hence to construct the matrix Mat(f) of a linear map f in some basis, we first put the images the basis vectors in a list.
# This list represents the TRANSPOSE of the desired matrix (because the images of the basis vectors by f must be the COLUMNS
# of Mat(f)), thus we apply TransposedMat to it and we get Mat(f).
# For instance consider the linear map mE for the left multiplication by E in \bar U_i(sl2). We have:
# E*B[1] = E*1 = E = B[5],  E*B[2] = E*K = EK = B[6],  E*B[3] = E*(K^2) = EK^2 = B[7],  E*B[4] = E*(K^3) = EK^3 = B[8],
# E*B[5] = E*E = 0 = zero16,  E*B[6] = E*(EK) = 0 = zero16,  E*B[7] = E*(EK^2) = 0 = zero16,  E*B[8] = E*(EK^3) = 0 = zero16,
# E*B[9] = E*F = EF = B[13],  E*B[10] = E*(FK) = EFK = B[14],  E*B[11] = E*(FK^2) = EFK^2 = B[15],  E*B[12] = E*(FK^3) = EFK^3 = B[16],
# E*B[13] = E*(EF) = 0 = zero16,  E*B[14] = E*(EFK) = 0 = zero16, E*B[15] = E*(EFK^2) = 0 = zero16, E*B[16] = E*(EFK^3) = 0 = zero16.
# Hence the matrix of mE is given by:          
mE := TransposedMat(
      [ B[5], B[6], B[7], B[8], 
        zero16, zero16, zero16, zero16, 
        B[13], B[14], B[15], B[16], 
        zero16, zero16, zero16, zero16] );;

# One similarly gets:
mF := TransposedMat(
      [ B[9], B[10], B[11], B[12], 
        B[13] + (i/2)*B[2] - (i/2)*B[4], B[14] + (i/2)*B[3] - (i/2)*B[1], B[15] + (i/2)*B[4] - (i/2)*B[2], B[16] + (i/2)*B[1] - (i/2)*B[3],
        zero16, zero16, zero16, zero16,
        -(i/2)*B[10] + (i/2)*B[12], -(i/2)*B[11] + (i/2)*B[9], -(i/2)*B[12] + (i/2)*B[10], -(i/2)*B[9] + (i/2)*B[11] ] );;
mK := TransposedMat(
      [ B[2], B[3], B[4], B[1],
        -B[6], -B[7], -B[8], -B[5],
        -B[10], -B[11], -B[12], -B[9],
        B[14], B[15], B[16], B[13] ] );;


## MULTIPLICATION TENSOR ##
# A matrix m of size 16x256. Since we use the monomial basis, m is easily constructed from the matrices mE, mF and mK constructed before.
m := TransposedMat(
     Concatenation(
     I_16, TransposedMat(mK), TransposedMat(mK^2), TransposedMat(mK^3), 
     TransposedMat(mE), TransposedMat(mE*mK), TransposedMat(mE*mK^2), TransposedMat(mE*mK^3), 
     TransposedMat(mF), TransposedMat(mF*mK), TransposedMat(mF*mK^2), TransposedMat(mF*mK^3), 
     TransposedMat(mE*mF), TransposedMat(mE*mF*mK), TransposedMat(mE*mF*mK^2), TransposedMat(mE*mF*mK^3) ) );;


## MATRICES FOR THE RIGHT MULTIPLICATIONS BY E, F, K ##
# We will also need the matrices rE, rF, rK for the right multiplications by the generators E, F, K in order to construct the
# centralizers of the (iterated) coproducts, see section "Computation of the cochain spaces below".
# The matrices rE, rF, rK can be constructed by hand as the matrices mE, mF, mK but it is easier to define them thanks to the
# multiplication tensor m just constructed:
rE := m*KroneckerProduct(I_16,TransposedMat([B[5]]));;
rF := m*KroneckerProduct(I_16,TransposedMat([B[9]]));;
rK := m*KroneckerProduct(I_16,TransposedMat([B[2]]));;


## COMULTIPLICATION TENSOR
# A matrix com of size 256x16.

tens := function(u, v)  ## function for the tensor product of two vectors in \bar U_i ##
   return Flat( KroneckerProduct([u], [v]) ); # KroneckerProduct can be applied to matrices only.
                                              # So we regard u and v as 1x16 matrices [u] and [v] and we compute their Kronecker product.
                                              # This gives a 1x256 matrix and we apply Flat to transform this matrix into a vector.
end;;

# With the function tens just defined we easily construct the coproduct tensor:
com := TransposedMat(
        [ tens(B[1], B[1]), tens(B[2], B[2]), tens(B[3], B[3]), tens(B[4], B[4]),
          tens(B[1], B[5]) + tens(B[5], B[2]), tens(B[2], B[6]) + tens(B[6], B[3]), tens(B[3], B[7]) + tens(B[7], B[4]), tens(B[4], B[8]) + tens(B[8], B[1]),
          tens(B[9], B[1]) + tens(B[4], B[9]), tens(B[10], B[2]) + tens(B[1], B[10]), tens(B[11], B[3]) + tens(B[2], B[11]), tens(B[12], B[4]) + tens(B[3], B[12]),
          tens(B[9], B[5]) + tens(B[4], B[13]) + tens(B[13], B[2]) - tens(B[8], B[10]),
          tens(B[10], B[6]) + tens(B[1], B[14]) + tens(B[14], B[3]) - tens(B[5], B[11]),
          tens(B[11], B[7]) + tens(B[2], B[15]) + tens(B[15], B[4]) - tens(B[6], B[12]),
          tens(B[12], B[8]) + tens(B[3], B[16]) + tens(B[16], B[1]) - tens(B[7], B[9]) ] );;


## DISPLAYING FUNCTION ##
# Recall that we want to compute explicit cocycles in H3_DY. According to equation (52) of the paper, such cocycles can
# naturally be seen as elements of the third tensor power (\bar U_i)^{\otimes 3} of \bar U_i, which is our point of view
# in this program.
# However, elements in this third tensor power are vectors of length 16^3 = 4096, which is by far too long to be easily read.
# So we need to display elements of (\bar U_i)^{\otimes 3} in a more readable form, which is the purpose of the function
# disp3 below. It takes a vector of length 4096 as input (representing an element in (\bar U_i)^{\otimes 3}) and returns a
# string which describes this vector as a linear combination of tensor products of basis vectors in the monomial basis 
# of \bar U_i.

dispB := ["1", "K", "K2", "K3", "E", "EK", "EK2", "EK3", "F", "FK", "FK2", "FK3", "EF", "EFK", "EFK2", "EFK3"];
# The list dispB contains string which represent the basis vectors in the monomial basis of \bar U_i.

disp3 := function(ttt)
   local jj, kk, ll, res, first;
   res := "0";
   first := true;
   for jj in [1 .. 16] do
      for kk in [1 .. 16] do
         for ll in [1 .. 16] do
            if ttt[16^2*(jj-1)+16*(kk-1)+ll] <> 0 then
               if first then
                  first := false;
                  res := "";
               else
                  res := Concatenation(res, " + ");
               fi;
               res := Concatenation(res, "(", String(ttt[16^2*(jj-1)+16*(kk-1)+ll]), ")*[", dispB[jj], ", ", dispB[kk],
                                    ", ", dispB[ll],"]");
            fi;
         od;
      od;
   od;
   return res;
end;


## USEFUL FUNCTION ##
# Given a basis and a list of coordinates in this basis, there a GAP function LinearCombination which returns the vector
# represented by these coordinates in this basis. This function is unpleasantly slow in spaces with big dimension.
# Here we define a function LC which does the same job much faster; the idea is simply to ignore the coordinates which 
# are equal to 0.

LC := function(base, coeffs)
   local index, res;
   res := 0*base[1];
   for index in [1 .. Length(coeffs)] do
      if coeffs[index] <> 0 then
         res := res + coeffs[index]*base[index];
      fi;
   od;
   return res;
end;


##################################################
## COMPUTATION OF THE COCHAIN SPACE IN DEGREE 2 ##
##################################################
Print("Computing the DY cochain space in degree 2...\n\n");
# Recall that thanks to equation (52) of the paper, the cochain spaces are isomorphic to the centralizers of the images of 
# the iterated coproducts. In particular, the cochain space in degree 2 C2_DY is the centralizer of the image of the coproduct.
# Since \bar U_i(sl2) is generated by E, F and K, an element x in (\bar U_i(sl2)) x (\bar U_i(sl2)) commutes with the image
# of the iterated coproduct if and only if it commutes with \Delta(E), \Delta(F) and \Delta(K). Hence we will compute the
# centralizers Z2_E, Z2_F and Z2_K of these coproducts and then compute the intersection of Z2_E, Z2_F, Z2_K. Again these problems
# reduce to linear equations.


tComE2 := [ ]; # This list will be the transpose matrix of the commutator [ -, \Delta(E) ]
               # We construct this (transpose) matrix by computing the value of the commutator on the vectors of the tensor
               # product basis of (\bar U_i(sl2) x \bar U_i(sl2)) (which consists of tensor products of monomials).
for j in [1 .. 16] do
   for k in [1 .. 16] do
      tComE2[16*(j-1) + k] := tens(B[j],mE*B[k]) + tens(mE*B[j],mK*B[k])
                            - tens(B[j],rE*B[k]) - tens(rE*B[j],rK*B[k]);
   od;
od;

Z2_E := Subspace(Qi^(16^2), NullspaceMat(tComE2), "basis");; # Since NullspaceMat computes the right kernel we apply it
                                                              # to the transpose of the matrix of the commutator, i.e. tComE2.
Unbind(tComE2); # We delete the matrix tComE2 to decrease the memory space.

tComF2 := [ ]; # This list will be the transpose matrix of the commutator [ -, \Delta(F) ]
for j in [1 .. 16] do
   for k in [1 .. 16] do
      tComF2[16*(j-1) + k] := tens(mF*B[j],B[k]) + tens(mK^(-1)*B[j],mF*B[k])
                            - tens(rF*B[j],B[k]) - tens(rK^(-1)*B[j],rF*B[k]);
   od;
od;

Z2_F := Subspace(Qi^(16^2), NullspaceMat(tComF2), "basis");; # Same remark as before.
Unbind(tComF2);

tComK2 := [ ];
for j in [1 .. 16] do
   for k in [1 .. 16] do
      tComK2[16*(j-1) + k] := tens(mK*B[j],mK*B[k]) - tens(rK*B[j],rK*B[k]);
   od;
od;

Z2_K := Subspace(Qi^(16^2), NullspaceMat(tComK2), "basis");;
Unbind(tComK2);

C2_DY := Intersection(Z2_E, Z2_F, Z2_K);;
Unbind(Z2_E); Unbind(Z2_F); Unbind(Z2_K); # We will no longer use these subspaces thus we can delete them.
BC2_DY := BasisVectors(Basis(C2_DY));; # A list of basis vectors for C2_DY.


##################################################
## COMPUTATION OF THE COCHAIN SPACE IN DEGREE 3 ##
##################################################
Print("Computing the DY cochain space in degree 3...\n\n");
# The cochain space in degree 3 C3_DY is the centralizer of the image of the iterated coproduct (\Delta x id)\Delta.
# We compute it exactly as we did for C2_DY, as the intersection of Z3_E, Z3_F and Z3_K which are respectively the 
# centralizers of the iterated coproducts of E, F and K.

tComE3 := [ ];
for j in [1 .. 16] do
   for k in [1 .. 16] do
      for l in [1 .. 16] do
         tComE3[16^2*(j-1) + 16*(k-1) + l] := tens(tens(B[j],B[k]),mE*B[l]) + tens(tens(B[j],mE*B[k]),mK*B[l]) 
                                              + tens(tens(mE*B[j],mK*B[k]),mK*B[l])
                                              - tens(tens(B[j],B[k]),rE*B[l]) - tens(tens(B[j],rE*B[k]),rK*B[l]) 
                                              - tens(tens(rE*B[j],rK*B[k]),rK*B[l]);
      od;
   od;
od;

# N.B. : tComE3 is a matrix of size 4096x4096 so its construction and the computation of its kernel can take a long time, 
# depending on the computer.

Z3_E := Subspace(Qi^(16^3), NullspaceMat(tComE3), "basis");;
Unbind(tComE3);

tComF3 := [ ];
for j in [1 .. 16] do
   for k in [1 .. 16] do
      for l in [1 .. 16] do
         tComF3[16^2*(j-1) + 16*(k-1) + l] := tens(tens(mF*B[j],B[k]),B[l]) + tens(tens(mK^(-1)*B[j],mF*B[k]),B[l]) 
                                              + tens(tens(mK^(-1)*B[j],mK^(-1)*B[k]),mF*B[l])
                                              - tens(tens(rF*B[j],B[k]),B[l]) - tens(tens(rK^(-1)*B[j],rF*B[k]),B[l]) 
                                              - tens(tens(rK^(-1)*B[j],rK^(-1)*B[k]),rF*B[l]);
      od;
   od;
od;

Z3_F := Subspace(Qi^(16^3), NullspaceMat(tComF3), "basis");;
Unbind(tComF3);

tComK3 := [ ];
for j in [1 .. 16] do
   for k in [1 .. 16] do
      for l in [1 .. 16] do
         tComK3[16^2*(j-1) + 16*(k-1) + l] := tens(tens(mK*B[j],mK*B[k]),mK*B[l])
                                              - tens(tens(rK*B[j],rK*B[k]),rK*B[l]); 
      od;
   od;
od;

Z3_K := Subspace(Qi^(16^3), NullspaceMat(tComK3), "basis");;
Unbind(tComK3);

C3_DY := Intersection(Z3_E, Z3_F, Z3_K);
Unbind(Z3_E); Unbind(Z3_F); Unbind(Z3_K);
BC3_DY := BasisVectors(Basis(C3_DY));; # A list of basis vectors for C3_DY


##################################################################################
## CONSTRUCTION OF THE DEGREE 2 DIFFERENTIAL d2_DY AND COMPUTATION OF ITS IMAGE ##
##################################################################################
Print("Constructing the degree 2 DY differential and computing its image...\n\n");
# Recall that we use the description of the DY complex given in equation (52) of the paper. The DY differential is given in
# equation (53) and uses the linear maps \Delta x id and id \otimes \Delta (where "x" means tensor product \otimes).
# Thus we first define functions which for an element tt in (\bar U_i(sl2)) x (\bar U_i(sl2)) return (\Delta x id)(tt) and 
# (id \otimes \Delta)(tt).

DeltaId := function(tt) ## returns (\Delta \otimes id)(tt)
   local jj, kk, res;
   res := tens(zero16, zero16);
   for jj in [1 .. 16] do
      for kk in [1 .. 16] do
         res := res + tt[16*(jj-1)+kk]*tens(com*B[jj],B[kk]);
      od;
   od;
   return res;
end;

IdDelta := function(tt) ## returns (id \otimes \Delta)(tt)
   local jj, kk, res;
   res := tens(zero16, zero16);
   for jj in [1 .. 16] do
      for kk in [1 .. 16] do
         res := res + tt[16*(jj-1)+kk]*tens(B[jj],com*B[kk]);
      od;
   od;
   return res;
end;

td2_DY := [ ]; # This list will be the transpose matrix of the differential d2_DY : C2_DY --> C3_DY.
               # We construct this (transpose) matrix by computing the value of d2_DY on the vectors of the basis of C2_DY
               # computed before.
for j in [1 .. Length(BC2_DY)] do
   td2_DY[j] := tens(B[1], BC2_DY[j]) - DeltaId(BC2_DY[j]) + IdDelta(BC2_DY[j]) - tens(BC2_DY[j],B[1]); # see formula of the differential in equation (53) of the paper.
od;

Imd2_DY := Subspace(Qi^(16^3), td2_DY);; # Image of the d2_DY (the command Subspace computes the vector subspace generated 
                                         # by the vectors in the list td2_DY, or equivalently the subspace generated by the
                                         # ROWS of the matrix td2_DY.
BImd2_DY := BasisVectors(Basis(Imd2_DY));; # A list of basis vectors for Imd2_DY.

Unbind(td2_DY); # We do not longer need this matrix in the sequel so we can delete it.


###################################################################################
## CONSTRUCTION OF THE DEGREE 3 DIFFERENTIAL d2_DY AND COMPUTATION OF ITS KERNEL ##
###################################################################################
Print("Constructing the degree 3 DY differential and computing its kernel... (may take about 5 minutes)\n\n");
# We construct the differential exactly as before. This time we need three functions which for a vector ttt in 
# (\bar U_i(sl2))^{\otimes 3} (third tensor power) return (\Delta x id x id)(ttt), (id x \Delta x id)(ttt) and 
# (id x id x \Delta)(ttt).

DeltaIdId := function(ttt) ## returns (\Delta \otimes id \otimes id)(ttt)
   local jj, kk, ll, res;
   res := tens(tens(zero16, zero16),zero16);
   for jj in [1 .. 16] do
      for kk in [1 .. 16] do
         for ll in [1 .. 16] do
            if ttt[16^2*(jj-1)+16*(kk-1)+ll] <> 0 then
               res := res + ttt[16^2*(jj-1)+16*(kk-1)+ll]*tens(tens(com*B[jj],B[kk]),B[ll]);
            fi;
         od;
      od;
   od;
   return res;
end;

IdDeltaId := function(ttt) ## returns (id \otimes \Delta \otimes id)(ttt)
   local jj, kk, ll, res;
   res := tens(tens(zero16, zero16),zero16);
   for jj in [1 .. 16] do
      for kk in [1 .. 16] do
         for ll in [1 .. 16] do
            if ttt[16^2*(jj-1)+16*(kk-1)+ll] <> 0 then
               res := res + ttt[16^2*(jj-1)+16*(kk-1)+ll]*tens(tens(B[jj],com*B[kk]),B[ll]);
            fi;
         od;
      od;
   od;
   return res;
end;

IdIdDelta := function(ttt) ## returns (id \otimes id \otimes \Delta)(ttt)
   local jj, kk, ll, res;
   res := tens(tens(zero16, zero16),zero16);
   for jj in [1 .. 16] do
      for kk in [1 .. 16] do
         for ll in [1 .. 16] do
            if ttt[16^2*(jj-1)+16*(kk-1)+ll] <> 0 then
               res := res + ttt[16^2*(jj-1)+16*(kk-1)+ll]*tens(tens(B[jj],B[kk]),com*B[ll]);
            fi;
         od;
      od;
   od;
   return res;
end;

td3_DY := [ ];
for j in [1 .. Length(BC3_DY)] do
   td3_DY[j] := tens(B[1], BC3_DY[j]) - DeltaIdId(BC3_DY[j]) + IdDeltaId(BC3_DY[j]) - IdIdDelta(BC3_DY[j]) + tens(BC3_DY[j],B[1]);
od;


F := NullspaceMat(td3_DY);; # Contains a basis of vectors of the kernel of d3_DY. But we have to be careful that the list F
                            # actually contains the COORDINATES of the  kernel vectors in the basis BC3_DY. Thus we have to
                            # convert these coordinates into actual vectors in Qi^(16^3) thanks to the function LC defined above.
BKerd3_DY := [ ]; # This list will contain a basis of ker(d3_DY).
for j in [1 .. Length(F)] do
   BKerd3_DY[j] := LC(BC3_DY, F[j]);
od;
Unbind(F);

Kerd3_DY := Subspace(Qi^(16^3), BKerd3_DY, "basis");


###############################################
## COMPUTATION OF EXPLICIT COCYCLES IN H3_DY ##
###############################################
Print("Computing the three explicit cocycles in H3_DY...\n\n");

# The idea is simply to find a basis, denoted Bcohom, of a supplementary of Imd2_DY in Kerd3_DY:
Bcohom := [ ];
Supp := Imd2_DY;
for j in [1 .. Length(BKerd3_DY)] do
   if not(BKerd3_DY[j] in Supp) then
      Supp := Supp + Subspace(Qi^(16^3), [ BKerd3_DY[j] ]);
      Add(Bcohom, BKerd3_DY[j]);
   fi;
od;

# Bcohom is a list of three vectors in Qi^(16^3). These vectors are of length 4096 and it is not possible for a human to
# read them in their coordinate form. To read them we use the function disp3 defined in the first part above.


#############
## RESULTS ##
#############
Print("A basis of cocycles for H3_DY is given by :\n\n", disp3(Bcohom[1]), "\n\n", disp3(Bcohom[2]), "\n\n", disp3(Bcohom[3]));
