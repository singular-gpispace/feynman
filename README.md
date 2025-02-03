# Feynman
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://singular-gpispace.github.io/Feynman/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://singular-gpispace.github.io/Feynman/

[ga-img]: https://github.com/singular-gpispace/Feynman/actions/workflows/CI.yml/badge.svg?branch=main
[ga-url]: https://github.com/singular-gpispace/Feynman/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/singular-gpispace/Feynman/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/singular-gpispace/Feynman

| **Documentation**                                                         | **Build Status**                                      |
|:-------------------------------------------------------------------------:|:-----------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][ga-img]][ga-url] [![][codecov-img]][codecov-url] |

# Feynman

The package Feynman generate the Feynman integrand associated to a Feynman diagram (graph) and compute an IBP system using the powerful module-intersection integration-by-parts (IBP) method, suitable for multi-loop and multi-scale Feynman integral reduction. It will provide( soon! ) an application programming interface(API) in OSCAR to use packages NeatIBP, pfd-parallel to make this computation much faster and to solve the reduction problem associated to Feynman integrals completely.The package Feynman is based on the computer algebra system OSCAR and is provided as a package for the Julia programming language.


This package can generate the Feynman integrand associated to a Feynman graph $G$ if $\text{number of internal edges} \leq \frac{1}{2}(l+1)+el$.
Here $l$ is the loop number and $e=span<p_1,...,p_{n_{ext}}>$ is the number of linearly independent external momenta of $G$.

In the case $\text{number of internal edges} < \frac{1}{2}(l+1)+el$, the package generate the Feynman integrand associated to a larger Feynman graph $G'$ so that $G$ is a subgraph of  $G'$. User can set appropriately the denomiator powers to zero to obtain the Baikov representation of $G$.

# Installation

We assume that Julia is installed in a recent enough version to run OSCAR. Navigate in a terminal to the folder where you want to install the package and pull the package from Github:

```bash
git pull https://github.com/singular-gpispace/Feynman.git
```

In the same folder execute the following command:

```bash
julia --project
```

This will activate the environment for our package. In Julia install missing packages:

```bash
import Pkg; Pkg.instantiate()
```

and load our package. On the first run this may take some time.

```bash
using Feynman  
```