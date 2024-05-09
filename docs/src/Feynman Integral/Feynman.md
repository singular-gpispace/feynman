# Feynman Integral

```@meta
CurrentModule = feynman
DocTestSetup = quote
  using feynman
end
```

```@setup TropicalFeynman
using feynman
```

## Graph

A Feynman graph is a (non-metrized) graph Γ without ends with n vertices which are labeled $x_1, . . . , x_n$ and with labeled edges $q_1, . . . , q_r$.
The graph $G$ is represented as a collection of vertices $V$ and edges $E$. Each edge is a pair $(v,w)$ where both $v$ and $w$ are elements of the set of vertices $V$.

