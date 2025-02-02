# Installation

We assume that SINGULAR is installed in a recent stable version (4.1.1 Feb 2018  or 4.3.2.0 Jan_2022 ) .

Run Singular in the terminal,

```singular
 >Singular
                     SINGULAR                                 /  Development
 A Computer Algebra System for Polynomial Computations       /   version 4.1.1
                                                           0<
 by: W. Decker, G.-M. Greuel, G. Pfister, H. Schoenemann     \   Feb 2018
FB Mathematik der Universitaet, D-67653 Kaiserslautern        \  Debian 1:4.1.1-p2+ds-4build2
> 
```

Then install the library feynman.lib.
```singular
> LIB "/path_to/Singular_Feynman/feynman.lib";
// ** redefining setMat (LIB "/path_to/Feynman/Singular_Feynman/feynman.lib";)
// ** redefining setMat (LIB "/path_to/Feynman/Singular_Feynman/feynman.lib";)
// ** loaded /path_to/Feynman/Singular_Feynman/feynman.lib (4.3.2.0,Jan_2022)
// ** loaded /usr/bin/../share/singular/LIB/general.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/matrix.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/nctools.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/random.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/ring.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/primdec.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/absfact.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/triang.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/elim.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/poly.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/inout.lib (4.1.1.0,Dec_2017)
// ** loaded /usr/bin/../share/singular/LIB/linalg.lib (4.1.1.0,Dec_2017)
> 
```

