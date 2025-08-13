\ Model qcp
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  116.6666666666667 C0 + 50 C1 + 50 C2 + 50 C3 + 66.66666666666666 C4
   + 0 aux_var + [ 200 C1 * C2 - 100 C2 * C3 + 233.3333333333333 C3 * C4
   ] / 2 
Subject To
 R0: 4 C5 + 3 C6 - 7 C7 >= 0
 R1: 6 C5 - 3 C6 + 7 C7 >= 0
 R2: - 2 C5 + C6 + C7 >= 0
 qc0: 4 C5 + 3 C6 + 3 C7 + [ - 10 C0 * C6 ] = 0
 qc1: - 6 C5 + 3 C6 + 3 C7 + [ - 10 C1 * C6 ] = 0
 qc2: - 6 C5 + 3 C6 + 3 C7 + [ - 10 C2 * C7 ] = 0
 qc3: 4 C5 + 3 C6 - 7 C7 + [ - 10 C3 * C5 ] = 0
 qc4: 4 C5 + 3 C6 - 7 C7 + [ - 10 C4 * C6 ] = 0
Bounds
 C0 <= 1
 C1 <= 1
 C2 <= 1
 C3 <= 1
 C4 <= 1
 -1 <= C5 <= 1
 -1 <= C6 <= 1
 -1 <= C7 <= 1
 aux_var >= 0.5
General Constraints
 1-norm: aux_var = NORM ( 1 ) ( C5 , C6 , C7 )
End
