# Redox Reaction Balancer

MOLECULES are written in the form:

  `[COEFF] X[n]Yz[m] * [#HYDRATES]H2O {[CHARGE][+/-]}`

For example:
  
  `CH3COO{-}     -> Acetate`

  `Mg(ClO3)2     -> Magnesium Chlorate`

  `2FeSO4*7H2O   -> 2 x Iron(II) Sulfate Heptahydrate`

Parentheses can be used to group ions, but only one level is currently supported.

For example, this is valid:

  `(CH3)2`

But this is not:

  `(C(H)3)2`

REACTIONS must be in the following form:

  `nX + nY -> nZ`

The arrow Regex is `/→|=+|-*>/`, meaning these are all valid:

  `2CN{-} + Au ---> Au(CN)2{-}`

  `NH3 + ClO{-} > N2H4 + Cl{-}`

  `AlH4{-} + H2CO === Al{3+} + CH3OH`

  `Hg + ZnSO4 → HgO4S + Zn`
