```@meta
CurrentModule = Feynman
```
# OVERVIEW:

We generate the Feynman integrand associated to a Feynman diagram and compute reduced IBP system to reduce given target integers to master integrals.

This package can generate the Feynman integrand associated to a Feynman graph $G$ if $\text{number of internal edges} \leq \frac{1}{2}(l+1)+el$.
Here $l$ is the loop number and $e=span<p_1,...,p_{n_{ext}}>$ is the number of linearly independent external momenta of $G$.

In the case $\text{number of internal edges} < \frac{1}{2}(l+1)+el$, the package generate the Feynman integrand associated to a larger Feynman graph $G'$ so that $G$ is a subgraph of  $G'$. User can set appropriately the denomiator powers to zero to obtained the Baikov representation of $G$.





