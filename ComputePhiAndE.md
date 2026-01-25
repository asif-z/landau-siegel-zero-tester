## Computation for $\phi$ and $E$ in paper

This file contains mathematica code that computes the constants in Corollary 6.4 of our paper. This is for debugging purposes; for a rigorous ball-arithmetic computation, see misc/rigorousBound.c

Define $\lambda$:

```wolfram
lambda = 1.1
```

Set up constants as in Proposition 5.2, as well as $\theta_i$, $A(n,m)$, and $B(n,m)$ as in the beginning of section 6

```wolfram
r = lambda/(10 Log[10])
CC = {2.97655, 5.13781, 9.49562}
theta = {Pi/2, ArcCos[(r + 1/8)/(r + 7/8)], 
  ArcCos[(r + 1/4)/(r + 7/8)], ArcCos[(r + 1/2)/(r + 7/8)], 
  ArcCos[(r + 3/4)/(r + 7/8)], 0}
A = Table[
  Integrate[
   Cos[x], {x, Pi - theta[[n + 1]], Pi - theta[[m + 1]]}], {n, 0, 
   5}, {m, 0, 5}]
B = Table[
  Integrate[(Cos[x])^2, {x, Pi - theta[[n + 1]], 
    Pi - theta[[m + 1]]}], {n, 0, 5}, {m, 0, 5}]
```

Compute $K_2,K_3$ as in Lemma 6.1

```wolfram
K = Table[(2^(6 - k) Log[CC[[4 - k]]] + 
      2 Log[15/8 + 2^(k - 4) + 
         2 r]) ((r + 2^(k - 5))/(Pi (r + 7/8)) A[[k, 
         k + 1]] + (B[[k, k + 1]])/Pi) - (2^(6 - k) Log[CC[[5 - k]]] +
       Log[15/8 + 2^(k - 4) + 
        2 r]) ((r + 2^(k - 4))/(Pi (r + 7/8)) A[[k, 
         k + 1]] + (B[[k, k + 1]])/Pi), {k, 2, 3}]
```

Compute $K_4,K_5$ as in Lemma 6.2

```wolfram
KK = Table[
  -(2^(k - 1) Log[CC[[k - 3]]] + 
       2 Log[7/8 + 2^(3 - k)]) ((r + 1 - 
          2^(2 - k))/(Pi (r + 7/8)) A[[k, k + 1]] + 
      B[[k, k + 1]]/Pi) + (2^(k - 1) Log[CC[[k - 2]]] + 
      Log[7/8 + 
        2^(3 - k)]) ((r + 1 - 2^(3 - k))/(Pi (r + 7/8)) A[[k, 
         k + 1]] + B[[k, k + 1]]/Pi), {k, 4, 5}
  ]
```

Compute the following constants as in Theorem 6.3, Corollary 6.4

```wolfram
K1 = (2 Log[CC[[3]]] + 1/8 Log[2 r + 2]) B[[1, 2]]/(Pi (r + 1/8))
KR = 2 (Log[Abs[2 + r - (r + 7/8) E^(I theta[[4]])]] - 
    Log[2 Pi]) ((r + 1/2)/(Pi (r + 7/8)) A[[4, 6]] + B[[4, 6]]/Pi)
rho = -(2 A[[1, 2]]/(Pi (r + 7/8)) + 2 B[[1, 2]]/(Pi (r + 1/8)))
K = 1.912/(Pi (r + 7/8)) + KR + K1 + K[[1]] + K[[2]] + KK[[1]] + 
  KK[[2]] - Log[2 Pi] + EulerGamma/2 + 1 + 1/2 PolyGamma[(3 + r)/2]
EE = K + rho Log[1 + 1/r]
```

Compute $\phi$ as in Theorem 6.3

```wolfram
phi = (r A[[2, 6]])/(Pi (r + 7/8)) + B[[1, 2]]/(Pi (8 r + 1)) + 
  B[[2, 6]]/Pi

```

