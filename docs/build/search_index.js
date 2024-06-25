var documenterSearchIndex = {"docs":
[{"location":"Overview.html#Overview","page":"Functions","title":"Overview","text":"","category":"section"},{"location":"Overview.html","page":"Functions","title":"Functions","text":"Modules = [Feynman]","category":"page"},{"location":"Overview.html","page":"Functions","title":"Functions","text":"Modules = [Feynman]","category":"page"},{"location":"Overview.html#Feynman.Feynman","page":"Functions","title":"Feynman.Feynman","text":"Feynman is a package for computing integration-by-part identities (IBPs) of a Feynman Integral associated to Feynman graph using module intesecation method. \nThis package also provides an interface of Oscar to use the packages NeatIBP developed using Singular and GPI-Space.\n\n\n\n\n\n","category":"module"},{"location":"Overview.html#Feynman.ISP-Tuple{labeledgraph}","page":"Functions","title":"Feynman.ISP","text":"ISP(G::labeledgraph)\n\nUSAGE:  ISP(G); ASSUME: G is a labeled graph. RETURN: idal containing the irreducible scalar products(ISPs), that is, those scalar product which are not linearly dependent on the propagators.\n\n#Examples\n\njulia> G=simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);\n\njulia> G=labelGraph(G,0);\n\njulia> Gelim=eliminateVariables(G);\n\njulia> ISP(Gelim)\nIdeal generated by\n  p[3]*q[1]\n  p[1]*q[2]\n\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.balancingIdeal-Tuple{labeledgraph}","page":"Functions","title":"Feynman.balancingIdeal","text":"balancingIdeal(G::labeledgraph)\n\nUSAGE: balancingIdeal(G) ASSUME: G is a labeled graph. RETURN: Ideal of balancing condition of the graph. i.e Ideal generated by the relation of the momentums which are obtained by applying momentum conservation law to external mementa,         and at each vertex; This is an ideal of the ring G.over.\n\n#Examples\n\njulia> G=simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);\njulia> G=labelGraph(G,0);\njulia> balancingIdeal(G)\n\nIdeal generated by\np[1] + p[2] + p[3] + p[4]\np[1] + q[1] + q[2]\np[2] - q[2] + q[3]\np[3] - q[1] + q[4]\np[4] - q[3] - q[4]\n\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.computeBaikovMatrix-Tuple{Any}","page":"Functions","title":"Feynman.computeBaikovMatrix","text":"computeBaikovMatrix(G::simple_graphgraph) USAGE:  computeBaikovMatrix(G); ASSUME: G is a graph, or G is a labled graph where redundant variables have been eliminated by the procedure eliminateVariables, and deleted from the          ring by the procedure removeElimVars. RETURN: a labeled graph G, where the computed Baikov matrix and the polynomial ring where baikovmatrix is defined are stored in G.baikovmatrix and G.baikovover respectively.\n\n#Examples ```julia julia> G=simple_graph([1,2,3,4,5,6,7],[(6,1),(6,4),(1,2),(3,7),(4,3),(2,7),(5,6),(7,5),1,2,3,4,5]);\n\njulia> G=Feynman.labelGraph(G,0);\n\njulia> G=Feynman.eliminateVariables(G);\n\njulia> G=Feynman.removeElimVars(G); QQMPolyRingElem[p[1], p[2], p[3], p[4], p[5], q[1], q[2]] julia> G=Feynman.computeBaikovMatrix(G); labels used for Gram matrix of external loop momenta: [\"p[1]p[2] => 1//2t[1]\"] [\"p[1]p[3] => 1//2t[2]\"] [\"p[2]p[3] => 1//2t[3]\"] [\"p[1]p[4] => 1//2t[4]\"] [\"p[2]p[4] => 1//2t[5]\"] [\"p[3]p[4] => -1//2t[1] - 1//2t[2] - 1//2t[3] - 1//2t[4] - 1//2t[5]\"] Assignment of Baikov variables (Z_i) are: [\"z[1] => p[3]q[1]\"] [\"z[2] => p[4]q[1]\"] [\"z[3] => q[1]^2\"] [\"z[4] => -2p[1]q[1] + q[1]^2\"] [\"z[5] => 2p[1]p[2] - 2p[1]q[1] - 2p[2]q[1] + q[1]^2\"] [\"z[6] => p[1]q[2]\"] [\"z[7] => q[2]^2\"] [\"z[8] => -2p[1]p[2] - 2p[1]p[3] - 2p[1]p[4] - 2p[2]p[3] - 2p[2]p[4] - 2p[3]q[2] - 2p[4]q[2] + q[2]^2\"] [\"z[9] => -2p[4]q[2] + q[2]^2\"] [\"z[10] => q[1]^2 + 2q[1]q[2] + q[2]^2\"] [\"z[11] => -2p[1]q[1] - 2p[1]q[2] - 2p[2]q[1] - 2p[2]q[2] - 2p[3]q[1] - 2p[3]q[2] - 2p[4]q[1] - 2p[4]q[2] + q[1]^2 + 2q[1]*q[2] + q[2]^2\"]\n\njulia> G.baikovmatrix 6×6 Matrix{RingElem}:  0                      …  z[6]  1//2t[1]                 1//2t[2] + 1//2t[3] + 1//2t[4] + 1//2t[5] - z[1] - z[2] - 1//2z[3] + 1//2z[5] - z[6] - 1//2z[7] + 1//2z[8] + 1//2z[10] - 1//2z[11]  1//2t[2]                 -1//2t[1] - 1//2t[2] - 1//2t[3] - 1//2t[4] - 1//2t[5] - 1//2z[8] + 1//2z[9]  1//2t[3]                 1//2z[7] - 1//2z[9]  1//2z[3] - 1//2z[4]     -1//2z[3] - 1//2z[7] + 1//2*z[10]  z[6]                   …  z[7]```\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.eliminateVariables-Tuple{labeledgraph}","page":"Functions","title":"Feynman.eliminateVariables","text":"eliminateVariables(G::labeledgraph)\n\nUSAGE:  eliminateVariables(G) ASSUME: G is a labeled graph. RETURN: labeled graph with variables of the bounded edges eliminated according to balancing condition.\n\n#Examples\n\njulia> G=simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);\njulia> G=labelGraph(G,0);\njulia> G=eliminateVariables(G);\njulia> printLabeledGraph(G);\n\nGraph with 4 vertices and 4 bounded edges 4 unbounded edges\nEdge terms:\n[\"(1, 3)=>q[1]\", \"(1, 2)=>-p[1] - q[1]\", \"(2, 4)=>-p[1] - p[2] - q[1]\", \"(3, 4)=>-p[3] + q[1]\", \"1=>p[1]\", \"2=>p[2]\", \"3=>p[3]\", \"4=>-p[1] - p[2] - p[3]\"]\njulia> G.elimvar\n\n4-element Vector{QQMPolyRingElem}:\n p[4]\n q[2]\n q[3]\n q[4]\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.feynmanDenominators-Tuple{labeledgraph}","page":"Functions","title":"Feynman.feynmanDenominators","text":"feynmanDenominators(G::labeledgraph)\n\nUSAGE:  feynmanDenominators(G); ASSUME: G is a labeled graph with the variables of the bounded edges eliminated according to balancing condition. i.e. G is a labeled graph where          the function eliminatedVariables applied. RETURN: ideal containing the propagators in the Feynman integral\n\n#Examples\n\njulia> G=simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);\njulia> G=labelGraph(G,0);\njulia> Gelim=eliminateVariables(G);\njulia> feynmanDenominators(Gelim)\nIdeal generated by\n  q[1]^2\n  p[1]^2 + 2*p[1]*q[1] + q[1]^2\n  p[1]^2 + 2*p[1]*p[2] + 2*p[1]*q[1] + p[2]^2 + 2*p[2]*q[1] + q[1]^2\n  p[3]^2 - 2*p[3]*q[1] + q[1]^2\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.makePoly-Tuple{Int64, Int64}","page":"Functions","title":"Feynman.makePoly","text":"makePoly(n::Int,m::Int)\n\nUSAGE:  makePoly(m,n); ASSUME: m and n are positve integers. RETURN: A polynomial ring with vatiables t[1],...,t[n],z[1],...,z[m] over QQ. \n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.propagators-Tuple{labeledgraph}","page":"Functions","title":"Feynman.propagators","text":"propagators(G::labeledgraph)\n\nUSAGE:  propagators(G) ASSUME: G is a labeld graph. RETURN: ideal, containing the denominators in the Feynman integral.\n\n#Examples\n\njulia> G=simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);\n\njulia> G=labelGraph(G,0);\n\njulia> Gelim=removeElimVars(Gelim);\n\njulia> Gelim=removeElimVars(Gelim);\n\njulia> propagators(Gelim)\nIdeal generated by\n  q[1]^2\n  q[2]^2\n  2*p[1]*p[3] + 2*p[1]*q[1] - 2*p[1]*q[2] + 2*p[3]*q[1] - 2*p[3]*q[2] + q[1]^2 - 2*q[1]*q[2] + q[2]^2\n  2*p[1]*q[1] + q[1]^2\n  -2*p[2]*q[1] + q[1]^2\n  2*p[1]*q[1] - 2*p[1]*q[2] + q[1]^2 - 2*q[1]*q[2] + q[2]^2\n  -2*p[2]*q[1] + 2*p[2]*q[2] + q[1]^2 - 2*q[1]*q[2] + q[2]^2\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.removeElimVars-Tuple{labeledgraph}","page":"Functions","title":"Feynman.removeElimVars","text":"removeElimVars(G::labeledgraph)\n\nUSAGE:  removeElimVars(G); ASSUME: G is a labled graph. RETURN: Removes the variables from G.elimvars. This key is generated by the procedure eliminatedVariables.\n\n#Examples ```julia julia> G=simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);\n\njulia> G=labelGraph(G,0);\n\njulia> Gelim=eliminateVariables(G);\n\njulia> G=removeElimVars(Gelim); QQMPolyRingElem[p[1], p[2], p[3], p[4], q[1], q[2]] julia> G.elimvar Any[]```\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.removeParameter-Tuple{AbstractAlgebra.Ring, Vector}","page":"Functions","title":"Feynman.removeParameter","text":"removeParameter(P::Ring,l::Vector)\n\nUSAGE: removeParameter(R,l); ASSUME: R is a polynomial ring. RETURN: polynomial ring with the parameters at indeces in l removed.\n\n#Examples\n\njulia> R,v,w=polynomial_ring(QQ,\"p\"=>(1:5),\"q\"=>(1:6));\njulia> I=ideal(R,[v[1],v[2],v[3],v[4],v[5],w[1]]);\njulia> u=complement_of_prime_ideal(I);\njulia> S,iso=localization(R,u);\njulia> S\n\njulia> S\nLocalization\n  of multivariate polynomial ring in 11 variables p[1], p[2], p[3], p[4], ..., q[6]\n    over rational field\n  at complement of prime ideal (p[1], p[2], p[3], p[4], p[5], q[1])\n\n\n  julia> removeParameter(S,[2]) \nLocalization\n  of multivariate polynomial ring in 10 variables p[1], p[3], p[4], p[5], ..., q[6]\n    over rational field\n  at complement of prime ideal (p[1], p[3], p[4], p[5], q[1])\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.removeVariable-Tuple{AbstractAlgebra.Ring, Vector}","page":"Functions","title":"Feynman.removeVariable","text":"removeVariable(R::Ring,l::Vector)\n\nUSAGE:  removeVariable(G,l) ASSUME: R is a polynomial ring. RERUTN: polynomial ring with the vaiables at indeces given in l removed.\n\n#Examples\n\njulia> R,v=polynomial_ring(QQ,\"p\"=>(1:3));\njulia> Feynman.removeVariable(R,[2])\nMultivariate polynomial ring in 2 variables p[1], p[3]\n  over rational field\n\n\n\n\n\n","category":"method"},{"location":"Overview.html#Feynman.removeVariableLocal-Tuple{AbstractAlgebra.Ring, Vector}","page":"Functions","title":"Feynman.removeVariableLocal","text":"removeVariableLocal(P::Ring,l::Vector)\n\nUSAGE: removeVariableLocal(P,l); ASSUME: P is a local ring locaized by the maximal ideal generated by parameters. RETURN: local ring where the variables at indeces in l removed. \n\n\n\n\n\n","category":"method"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"CurrentModule=Feynman","category":"page"},{"location":"Example.html#Example-:-Fully-massless-nonplanar-double-pentagon","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"","category":"section"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"(Image: alt text)","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"To provide an example on how to use our package, we calculate the Baikov matrix of the fully massless nonplanar double pentagon. ","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"We define the graph G from the list of vertices and list of edges. The direction of momenta are taken from the direction  of edges. All external momenta are taken to be outgoing.","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=simple_graph([1,2,3,4,5,6,7],[(6,1),(6,4),(1,2),(3,7),(4,3),(2,7),(5,6),(7,5),1,2,3,4,5]);","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"We then assign polynomial variables q[i] at bounded edges and function field variables P[i] at the unbounded edges over a prime filed of characteristic 0.","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"\njulia> G=labelGraph(G,0);\n\njulia> printLabeledGraph(G);\nGraph with 7 vertices and 8 bounded edges 5 unbounded edges\nEdge terms:\n[\"(6, 1)=>q[1]\", \"(6, 4)=>q[2]\", \"(1, 2)=>q[3]\", \"(3, 7)=>q[4]\", \"(4, 3)=>q[5]\", \"(2, 7)=>q[6]\", \"(5, 6)=>q[7]\", \"(7, 5)=>q[8]\", \"1=>p[1]\", \"2=>p[2]\", \"3=>p[3]\", \"4=>p[4]\", \"5=>p[5]\"]\n","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"Then we use balancing condition of the graph (realtions of momenta which are obtained by applying momentum conservation law at each vertex of the graph and to the whole graph) to rewrite each dependent momenta in terms of the eliments in the ordered set V of external momenta and loop momenta. Here we use invlex ordering on p[1],...,p[E],q[1],...,q[L]$ to choose independent external momenta and independent loop momenta. G.elimVars will store the eliminated variables.","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=eliminateVariables(G);\njulia> printLabeledGraph(G);\nGraph with 7 vertices and 8 bounded edges 5 unbounded edges\nEdge terms:\n[\"(6, 1)=>q[1]\", \"(6, 4)=>q[2]\", \"(1, 2)=>-p[1] + q[1]\", \"(3, 7)=>-p[3] - p[4] + q[2]\", \"(4, 3)=>-p[4] + q[2]\", \"(2, 7)=>-p[1] - p[2] + q[1]\", \"(5, 6)=>q[1] + q[2]\", \"(7, 5)=>-p[1] - p[2] - p[3] - p[4] + q[1] + q[2]\", \"1=>p[1]\", \"2=>p[2]\", \"3=>p[3]\", \"4=>p[4]\", \"5=>-p[1] - p[2] - p[3] - p[4]\"]\njulia> G.elimvar\n7-element Vector{QQMPolyRingElem}:\n p[5]\n q[3]\n q[4]\n q[5]\n q[6]\n q[7]\n q[8]","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"Then the irreducible scalar products associated to G can be printed as follows:","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia>ISP(G)\nIdeal generated by\n  p[3]*q[1]\n  p[4]*q[1]\n  p[1]*q[2]","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"We should remove the eliminated variables from G, in order to compute the Baikov matrix of G.","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=removeElimVars(G);\nQQMPolyRingElem[p[1], p[2], p[3], p[4], p[5], q[1], q[2]]","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"We then calculate the Baikov matrix associated to Feynman integral of G. It will also print the assignment of Baikov variables z[i] to each inverse propagators and irreducible scalar products of G.","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=computeBaikovMatrix(G);\nlabels used for Gram matrix of external loop momenta:\n[\"p[1]*p[2] => 1//2*t[1]\"]\n[\"p[1]*p[3] => 1//2*t[2]\"]\n[\"p[2]*p[3] => 1//2*t[3]\"]\n[\"p[1]*p[4] => 1//2*t[4]\"]\n[\"p[2]*p[4] => 1//2*t[5]\"]\n[\"p[3]*p[4] => -1//2*t[1] - 1//2*t[2] - 1//2*t[3] - 1//2*t[4] - 1//2*t[5]\"]\nAssignment of Baikov variables (Z_i) are:\n[\"z[1] => p[3]*q[1]\"]\n[\"z[2] => p[4]*q[1]\"]\n[\"z[3] => q[1]^2\"]\n[\"z[4] => -2*p[1]*q[1] + q[1]^2\"]\n[\"z[5] => 2*p[1]*p[2] - 2*p[1]*q[1] - 2*p[2]*q[1] + q[1]^2\"]\n[\"z[6] => p[1]*q[2]\"]\n[\"z[7] => q[2]^2\"]\n[\"z[8] => -2*p[1]*p[2] - 2*p[1]*p[3] - 2*p[1]*p[4] - 2*p[2]*p[3] - 2*p[2]*p[4] - 2*p[3]*q[2] - 2*p[4]*q[2] + q[2]^2\"]\n[\"z[9] => -2*p[4]*q[2] + q[2]^2\"]\n[\"z[10] => q[1]^2 + 2*q[1]*q[2] + q[2]^2\"]\n[\"z[11] => -2*p[1]*q[1] - 2*p[1]*q[2] - 2*p[2]*q[1] - 2*p[2]*q[2] - 2*p[3]*q[1] - 2*p[3]*q[2] - 2*p[4]*q[1] - 2*p[4]*q[2] + q[1]^2 + 2*q[1]*q[2] + q[2]^2\"]\n\njulia> G.baikovmatrix\n6×6 Matrix{RingElem}:\n 0                      …  z[6]\n 1//2*t[1]                 1//2*t[2] + 1//2*t[3] + 1//2*t[4] + 1//2*t[5] - z[1] - z[2] - 1//2*z[3] + 1//2*z[5] - z[6] - 1//2*z[7] + 1//2*z[8] + 1//2*z[10] - 1//2*z[11]\n 1//2*t[2]                 -1//2*t[1] - 1//2*t[2] - 1//2*t[3] - 1//2*t[4] - 1//2*t[5] - 1//2*z[8] + 1//2*z[9]\n 1//2*t[3]                 1//2*z[7] - 1//2*z[9]\n 1//2*z[3] - 1//2*z[4]     -1//2*z[3] - 1//2*z[7] + 1//2*z[10]\n z[6]                   …  z[7]","category":"page"},{"location":"Example.html#Example:-Fully-massless-planar-double-box","page":"Example : Fully massless nonplanar double pentagon","title":"Example: Fully massless planar double box","text":"","category":"section"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"(Image: alt text)","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=Feynman.simple_graph([1,2,3,4,5,6],[(6,1),(4,6),(1,2),(3,5),(4,3),(2,5),(5,6),1,2,3,4]);","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=Feynman.labelGraph(G,0);","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"G=eliminateVariables(G);","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> printLabeledGraph(G);\nGraph with 6 vertices and 7 bounded edges 4 unbounded edges\nEdge terms:\n[\"(6, 1)=>q[1]\", \"(4, 6)=>q[2]\", \"(1, 2)=>-p[1] + q[1]\", \"(3, 5)=>p[1] + p[2] - q[2]\", \"(4, 3)=>p[1] + p[2] + p[3] - q[2]\", \"(2, 5)=>-p[1] - p[2] + q[1]\", \"(5, 6)=>q[1] - q[2]\", \"1=>p[1]\", \"2=>p[2]\", \"3=>p[3]\", \"4=>-p[1] - p[2] - p[3]\"]","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G.elimvar\n6-element Vector{QQMPolyRingElem}:\n p[4]\n q[3]\n q[4]\n q[5]\n q[6]\n q[7]","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=removeElimVars(G);\nQQMPolyRingElem[p[1], p[2], p[3], p[4], q[1], q[2]]","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G=computeBaikovMatrix(G);\nlabels used for Gram matrix of external loop momenta:\n[\"p[1]*p[2] => 1//2*t[1]\"]\n[\"p[1]*p[3] => 1//2*t[2]\"]\n[\"p[2]*p[3] => -1//2*t[1] - 1//2*t[2]\"]\nAssignment of Baikov variables (Z_i) are:\n[\"z[1] => p[3]*q[1]\"]\n[\"z[2] => q[1]^2\"]\n[\"z[3] => -2*p[1]*q[1] + q[1]^2\"]\n[\"z[4] => 2*p[1]*p[2] - 2*p[1]*q[1] - 2*p[2]*q[1] + q[1]^2\"]\n[\"z[5] => p[1]*q[2]\"]\n[\"z[6] => q[2]^2\"]\n[\"z[7] => 2*p[1]*p[2] - 2*p[1]*q[2] - 2*p[2]*q[2] + q[2]^2\"]\n[\"z[8] => -2*p[1]*q[2] - 2*p[2]*q[2] - 2*p[3]*q[2] + q[2]^2\"]\n[\"z[9] => q[1]^2 - 2*q[1]*q[2] + q[2]^2\"]","category":"page"},{"location":"Example.html","page":"Example : Fully massless nonplanar double pentagon","title":"Example : Fully massless nonplanar double pentagon","text":"julia> G.baikovmatrix\n5×5 Matrix{RingElem}:\n 0                      1//2*t[1]                                 …  z[5]\n 1//2*t[1]              0                                            1//2*t[1] - z[5] + 1//2*z[6] - 1//2*z[7]\n 1//2*t[2]              -1//2*t[1] - 1//2*t[2]                       -1//2*t[1] + 1//2*z[7] - 1//2*z[8]\n 1//2*z[2] - 1//2*z[3]  1//2*t[1] + 1//2*z[3] - 1//2*z[4]            1//2*z[2] + 1//2*z[6] - 1//2*z[9]\n z[5]                   1//2*t[1] - z[5] + 1//2*z[6] - 1//2*z[7]     z[6]\n","category":"page"},{"location":"FeynmanIBP.html","page":"-","title":"-","text":"CurrentModule = Feynman","category":"page"},{"location":"index.html","page":"Feynman","title":"Feynman","text":"CurrentModule = Feynman","category":"page"},{"location":"index.html#Feynman","page":"Feynman","title":"Feynman","text":"","category":"section"},{"location":"index.html","page":"Feynman","title":"Feynman","text":"Documentation for Feynman.","category":"page"},{"location":"index.html","page":"Feynman","title":"Feynman","text":"The package Feynman computes complete set of IBP identities of the Feynman integral associated to a given Feynman graph using the powerful module-intersection integration-by-parts (IBP) method, suitable for multi-loop and multi-scale Feynman integral reduction. It will provide an application programming interface(API) in OSCAR to use packages NeatIBP, pfd-parallel to make this computation much faster.The package Feynman is based on the computer algebra system OSCAR and is provided as a package for the Julia programming language.","category":"page"},{"location":"Installation.html#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"We assume that Julia is installed in a recent enough version to run OSCAR. Navigate in a terminal to the folder where you want to install the package and pull the package from Github:","category":"page"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"git pull https://github.com/singular-gpispace/Feynman.git","category":"page"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"In the same folder execute the following command:","category":"page"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"julia --project","category":"page"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"This will activate the environment for our package. In Julia install missing packages:","category":"page"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"import Pkg; Pkg.instantiate()","category":"page"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"and load our package. On the first run this may take some time.","category":"page"},{"location":"Installation.html","page":"Installation","title":"Installation","text":"using Feynman  ","category":"page"}]
}
