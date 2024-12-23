using Oscar
#include("DataTypes.jl")

export greet
export printNet
export printGraph
export printLabeledGraph
export labelGraph
export balancingIdeal
export substituteGraph
export eliminateVariables
export removeVariable
export removeParameter
export feynmanDenominators
export propagators
export ISP
export removeElimVars
export computeBaikovMatrix
export makePoly
export removeVariableLocal
export computeIBP
export printIBP
export computeM1
export enumerateTerms
export Heviside
export power_monomials
export generating_probes
export oneRR
export manyRR
export pickDen
export evaluateBB
export pickWeights
export numeratorAnsatz
export interpolateAnsatz
export getWeights
export setBB
export getCandidateDenominators
export getIndDen
export interpolateAnsatzTemp
export computeCoef
#----------------------------------
export sample1
export Net
export simple_graph
export labeledgraph
export IBP
mutable struct sample1
    a::Int
end
mutable struct simple_graph
    vertices :: Vector{Int64}
    edges :: Vector
end

mutable struct Net
    rows :: Vector
end



mutable struct labeledgraph
    vertices :: Vector{Int64}
    edges :: Vector
    over :: Ring
    labels :: Vector
    overpoly ::Ring
    elimvar ::Vector
    baikovover :: Ring
    baikovmatrix :: Matrix{RingElem}
end

mutable struct IBP
    baikovind::Vector
    baikovover::Ring
    paraind::Int64
    setIBP::Vector
end

#----------------------------------






function greet()
    return "Hello"
end



function printNet(N::Vector)
    println('[',join(repr.(N),", "),']');
end

## This is the print function to print a simple graph
function printGraph(G::simple_graph)
    num_ver=length(G.vertices);
    num_edg=length(G.edges);
   
    ct=0;
        for i in 1:num_edg
            if length(G.edges[i])==1
                ct=ct+1;
            end    
        end

    if ct!=0
        println("Graph with ", num_ver ," vertices and ",num_edg-ct ," bounded edges ",ct," unbounded edges")
    else
        println("Graph with ", num_ver ," vertices and ",num_edg ," edges")
    end
end
#Make simple graph from list of vertices and edges

@doc raw"""
printLabeledGraph(G::labeledgraph)

USAGE: printLabeledGraph(G);
ASSUME: G is a labeled graph.
Theory: This is the print function used in julia to print a labeled graph.

#Examples
```julia
julia> var=["x","y","z","p","q","r"];
julia> R, (x,y,z,p,q,r)=polynomial_ring(QQ,var);
julia> G=labeledgraph([1,2,3,4],[(1,3),(1,2),(1,2),(2,4),(3,4),(3,4)],R,var,R,[3,4],R,[[y+1,R(2)] [R(3),r+3]]);
julia> printLabeledGraph(G);

Graph with 4 vertices and 6 edges
Edge terms:
["(1, 3)=>x", "(1, 2)=>y", "(1, 2)=>z", "(2, 4)=>p", "(3, 4)=>q", "(3, 4)=>r"]

```
"""

function printLabeledGraph(G::labeledgraph)
   
    printGraph(simple_graph(G.vertices,G.edges));
ev=String[];
println("Edge terms:");
for i in 1:length(G.edges)
    t=string(string(G.edges[i]),"=>",string(G.labels[i]));
    push!(ev,t) ;    
end
printNet(ev);
end

@doc raw"""
labelGraph(G::simple_graph,ch::Int)

USAGE: labelGraph(G,ch)
ASSUME: G is a graph and ch is either zero or prime.
RETURN: labeled graph with polynomialvariables qi at the bounded edges and functin filed variables pi at the unbounded edges over a prime filed of characteristic ch
Initially we it sets the fields Baikovmatrix and elimvar empty.
 
#Examples
```julia
julia> G3=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
julia> G4=Feynman.labelGraph(G3,0);
julia> Feynman.printLabeledGraph(G4);
Graph with 4 vertices and 4 bounded edges 4 unbounded edges
Edge terms:
["(1, 3)=>q[1]", "(1, 2)=>q[2]", "(2, 4)=>q[3]", "(3, 4)=>q[4]", "1=>p[1]", "2=>p[2]", "3=>p[3]", "4=>p[4]"]
```

"""

function labelGraph(G::simple_graph,ch::Int)
    if length(G.edges)==0
        error("Graph must have vertices")
    elseif ch<0

        error("ch must be nonnegative integer")
    
    end
    ct=0;
    for i in 1:length(G.edges)
        if length(G.edges[i])==1
            ct=ct+1;
        end    
    end
    anzq=length(G.edges)-ct;
    
    ##P=(Z/pZ)[p(1),...,p(ct),q(1),...q(anzq)]
    if ch==0
        P,p,q=polynomial_ring(QQ,"p"=>(1:ct),"q"=>(1:anzq));
        #P,q,p=polynomial_ring(QQ,"q"=>(1:anzq),"p"=>(1:ct));
    else
        F=GF(ch);
        P,p,q=polynomial_ring(F,"p"=>(1:ct),"q"=>(1:anzq));
        #P,q,p=polynomial_ring(F,"q"=>(1:anzq),"p"=>(1:ct));
    end
    
    
   
    ##Q=complement_of_prime_ideal(ideal(P,p));
    Q=complement_of_prime_ideal(ideal(P,p));
    R,iso=localization(P,Q);

    ##Making list of labels
    pidx=1;
    qidx=1;
    lab=[];
    p=gens(R);
    for i in 1:length(G.edges)
        if length(G.edges[i])==1
            #lab[i]=p[pidx];
            push!(lab,p[pidx]);
            pidx=pidx+1;
        else
            #lab[i]=q[qidx];
            push!(lab,q[qidx]);
            qidx=qidx+1;
        end 
        
    end
    lG=labeledgraph(G.vertices,G.edges,R,lab,P,[],P,[;;]);
    return lG
end

@doc raw"""
balancingIdeal(G::labeledgraph)

USAGE: balancingIdeal(G)
ASSUME: G is a labeled graph.
RETURN: Ideal of balancing condition of the graph. i.e Ideal generated by the relation of the momentums which are obtained by applying momentum conservation law to external mementa,
        and at each vertex; This is an ideal of the ring G.over.

#Examples
```julia
julia> G=simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
julia> G=labelGraph(G,0);
julia> balancingIdeal(G)

Ideal generated by
p[1] + p[2] + p[3] + p[4]
p[1] + q[1] + q[2]
p[2] - q[2] + q[3]
p[3] - q[1] + q[4]
p[4] - q[3] - q[4]

```

"""
function balancingIdeal(G::labeledgraph)
 v=G.vertices;
 e=G.edges;
 lab=G.labels;
 R=G.over;
 edg=R(0);
 rel=R(0);
 w=Vector{typeof(R(1))}(undef,0);
 for i in 1:length(v)
    edg=R(0);#changed
    for j in 1:length(e)
        if length(e[j])==2
            if e[j][1]==v[i]
                edg=edg+R(lab[j]);
                
            end
            if e[j][2]==v[i]
                edg=edg-R(lab[j]);
            end
        else
            if e[j][1]==v[i]
                edg=edg+R(lab[j]);
            end
        end 

        
    end 
   push!(w,edg); #changed
 end   
 for j in 1:length(e)
    if length(e[j])==1
        rel=rel+R(lab[j]);
        
    end 
    
 end
 I=ideal(R,rel);#changed
 J=ideal(R,w);
 return I+J;
end

@doc raw""" 
substituteGraph(G::labeledgraph,a::RingElement,b::RingElement)
USAGE:substituteGraph(G,a,b)
ASSUME: G is a labeled graph
RETURN:a labelled graph with labelling where each 'a' is replaced 'b'
"""

function substituteGraph(G::labeledgraph,a::RingElement,b::RingElement)
    L=G.labels;
    for i in 1:length(L)
       if L[i]==a
        L[i]=b;       
       end         
    end
   G.labels=L;
   return G; 
    
end


@doc raw"""
eliminateVariables(G::labeledgraph)

USAGE:  eliminateVariables(G)
ASSUME: G is a labeled graph.
RETURN: labeled graph with variables of the bounded edges eliminated according to balancing condition.

#Examples
```julia
julia> G=simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
julia> G=labelGraph(G,0);
julia> G=eliminateVariables(G);
julia> printLabeledGraph(G);

Graph with 4 vertices and 4 bounded edges 4 unbounded edges
Edge terms:
["(1, 3)=>q[1]", "(1, 2)=>-p[1] - q[1]", "(2, 4)=>-p[1] - p[2] - q[1]", "(3, 4)=>-p[3] + q[1]", "1=>p[1]", "2=>p[2]", "3=>p[3]", "4=>-p[1] - p[2] - p[3]"]
julia> G.elimvar

4-element Vector{QQMPolyRingElem}:
 p[4]
 q[2]
 q[3]
 q[4]
```
"""
function eliminateVariables(G::labeledgraph)
  RP=G.overpoly;
  R=G.over;
  I=balancingIdeal(G); #Ideal in G.over
  n=ngens(I);
  #Need ring homomorphism from R to RP to get ideal 
  H=hom(R,RP,gens(RP));
  T=hom(RP,R,gens(R));
  v=[];
  for i in 1:n
    push!(v,H(I[i])); 
  end
  I=ideal(RP,v);

  #------Getting variables and parameters
  p=gens(base_ring(R));
  para=Vector{typeof(RP(1))}(undef,0);
  var=Vector{typeof(RP(1))}(undef,0);
  
  for i in 1:length(p) 
  
      if p[i] in inverted_set(R)
        push!(var,p[i]);
    
      else
          push!(para,p[i]); 
          
      end
        
  end

#--------------------------------------
#od=invlex(var)*invlex(para);
#od=lex(para,var);
od=invlex(RP);
#--------------------------------------
  I=standard_basis(I,ordering=od,complete_reduction=true);
#  I=standard_basis(I,ordering=od);
  #G1=G;
  eliminatedVariables=Vector{typeof(RP(1))}(undef,0);
  G.elimvar=Vector{typeof(RP(1))}(undef,0);
  n=length(I);
  for i in 1:n
    ld=leading_term(I[i],ordering=od);
    ta=ld-I[i];
    #ld=RP(ld);

    if length(G.elimvar)==0
        push!(G.elimvar,H(T(ld)));
    else
        k=0;
        for j in 1:length(G.elimvar)
            if G.elimvar[j]==ld || G.elimvar[j]==-ld 
                k=k+1;
                #print(eliminatedVariables);
            end 
            
        end
        if k==0
            push!(G.elimvar,H(T(ld)));            
        end
        
    end
    #ta=RP(ta); 
    G=substituteGraph(G,ld,ta);  
    ld=0;
    ta=0;
  end  
  #G.elimvar=Vector{typeof(RP(1))}(undef,0);
  #G.elimvar=eliminatedVariables;
  return G;
end


@doc raw"""
removeVariable(R::Ring,l::Vector)

USAGE:  removeVariable(G,l)
ASSUME: R is a polynomial ring.
RERUTN: polynomial ring with the vaiables at indeces given in l removed.

#Examples
```julia
julia> R,v=polynomial_ring(QQ,"p"=>(1:3));
julia> Feynman.removeVariable(R,[2])
Multivariate polynomial ring in 2 variables p[1], p[3]
  over rational field
```

"""
function removeVariable(R::Ring,l::Vector)
    nv=0;
    for j in 1:length(l) 
        if l[j]<0|| l[j]>ngens(R)
            throw(error("Index out of range"));
        end
    end
    
    genR=deleteat!(gens(R),l);
    v=String[];
    for i in 1:length(genR)
        push!(v,string(genR[i]))       
    end
    R,c=polynomial_ring(coefficient_ring(R),v);
    return R
end


@doc raw"""
removeVariableLocal(P::Ring,l::Vector)

USAGE: removeVariableLocal(P,l);
ASSUME: P is a local ring locaized by the maximal ideal generated by parameters.
RETURN: local ring where the variables at indeces in l removed. 

"""
function removeVariableLocal(P::Ring,l::Vector)

    R=base_ring(P);
    nv=0;
    genR=gens(R);
    for j in 1:length(l)
        if l[j]<0|| l[j]>ngens(R) 
            throw(error("Index out of range"));
        end 
        
    end

    genR=deleteat!(gens(R),l);
        #print(genR);

    v=String[];
    for i in 1:length(genR)
        push!(v,string(genR[i]))       
    end
    #delete jth variable from base ring
    S,c=polynomial_ring(coefficient_ring(R),v);
    #-------
    u=[];
    k=0;
    genR=gens(R);
    for i in 1:length(genR) 
        #J=ideal(R,u);
        #L=complement_of_prime_ideal(J);
        if genR[i] in inverted_set(P) 
            
            continue;
        else
            push!(u,c[i]); 
            k=k+1;
        end
       
    end
    #-------
    #x=[];
   # for i in 1:k 
    #    push!(x,c[i]);
    #end
    #-------
    #w=deleteat!(x,j);
    J=ideal(S,u);
    U=complement_of_prime_ideal(J);
    T,iso=localization(S,U);
    return T;
    
end


@doc raw"""
removeParameter(P::Ring,l::Vector)

USAGE: removeParameter(R,l);
ASSUME: R is a polynomial ring.
RETURN: polynomial ring with the parameters at indeces in l removed.

#Examples
```julia
julia> R,v,w=polynomial_ring(QQ,"p"=>(1:5),"q"=>(1:6));
julia> I=ideal(R,[v[1],v[2],v[3],v[4],v[5],w[1]]);
julia> u=complement_of_prime_ideal(I);
julia> S,iso=localization(R,u);
julia> S

julia> S
Localization
  of multivariate polynomial ring in 11 variables p[1], p[2], p[3], p[4], ..., q[6]
    over rational field
  at complement of prime ideal (p[1], p[2], p[3], p[4], p[5], q[1])


  julia> removeParameter(S,[2]) 
Localization
  of multivariate polynomial ring in 10 variables p[1], p[3], p[4], p[5], ..., q[6]
    over rational field
  at complement of prime ideal (p[1], p[3], p[4], p[5], q[1])
```
"""
function removeParameter(P::Ring,l::Vector)
    R=base_ring(P);
    genR=gens(R);
    nv=0;
    for j in 1:length(l)
        if l[j]<0|| l[j]>ngens(R)
            throw(error("Index out of range"));
        end
    end
    
    genR=deleteat!(gens(R),l);
    v=String[];
    for i in 1:length(genR)
        push!(v,string(genR[i]))       
    end
    #delete jth variable from base ring
    S,c=polynomial_ring(coefficient_ring(R),v);
    #-------
    u=[];
    k=0;
    for i in 1:length(genR) 
        #J=ideal(R,u);
        #L=complement_of_prime_ideal(J);
        if genR[i] in inverted_set(P) 
            continue;
        else
            push!(u,c[i]); 
            k=k+1; 
        end
         
    end
    #-------
    #x=[];
   # for i in 1:k 
    #    push!(x,c[i]);
    #end
    #-------
    #w=deleteat!(x,j);
    #w=deleteat!(u,j);
    J=ideal(S,u);
    U=complement_of_prime_ideal(J);
    T,iso=localization(S,U);
    return T;
end

@doc raw"""
feynmanDenominators(G::labeledgraph)

USAGE:  feynmanDenominators(G);
ASSUME: G is a labeled graph with the variables of the bounded edges eliminated according to balancing condition. i.e. G is a labeled graph where 
        the function eliminatedVariables applied.
RETURN: ideal containing the propagators in the Feynman integral

#Examples
```julia
julia> G=simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
julia> G=labelGraph(G,0);
julia> Gelim=eliminateVariables(G);
julia> feynmanDenominators(Gelim)
Ideal generated by
  q[1]^2
  p[1]^2 + 2*p[1]*q[1] + q[1]^2
  p[1]^2 + 2*p[1]*p[2] + 2*p[1]*q[1] + p[2]^2 + 2*p[2]*q[1] + q[1]^2
  p[3]^2 - 2*p[3]*q[1] + q[1]^2
```

"""
function feynmanDenominators(G::labeledgraph)
    L=G.labels;
    S=G.over;
    v=[];
    for i in 1:length(L)
        if length(G.edges[i])==2
            push!(v,L[i]^2);    
        end
                
    end
    return ideal(S,v);
end


@doc raw"""
propagators(G::labeledgraph)

USAGE:  propagators(G)
ASSUME: G is a labeld graph.
RETURN: ideal, containing the denominators in the Feynman integral.

#Examples
```julia
julia> G=simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);

julia> G=labelGraph(G,0);

julia> Gelim=removeElimVars(Gelim);

julia> Gelim=removeElimVars(Gelim);

julia> propagators(Gelim)
Ideal generated by
  q[1]^2
  q[2]^2
  2*p[1]*p[3] + 2*p[1]*q[1] - 2*p[1]*q[2] + 2*p[3]*q[1] - 2*p[3]*q[2] + q[1]^2 - 2*q[1]*q[2] + q[2]^2
  2*p[1]*q[1] + q[1]^2
  -2*p[2]*q[1] + q[1]^2
  2*p[1]*q[1] - 2*p[1]*q[2] + q[1]^2 - 2*q[1]*q[2] + q[2]^2
  -2*p[2]*q[1] + 2*p[2]*q[2] + q[1]^2 - 2*q[1]*q[2] + q[2]^2
```
"""
function propagators(G::labeledgraph)
    L=G.labels;
    S=G.over;
    RP=G.overpoly;
   
    u=[];
    for i in 1:length(L)
        if length(G.edges[i])==2
                    push!(u,L[i]^2);                 
        end
                
    end
    J=ideal(S,u);
    
    v=[];
    for i in 1:length(L)
        if length(G.edges[i])==1          
                    push!(v,L[i]^2);  
        end             
    end
    infedges=ideal(S,v);

#ring Map from RP to S
    H=hom(S,RP,gens(RP));
#Ideal I an an Ideal of RP
  v=[];
  u=[];
  for i in 1:ngens(J)
    push!(v,H(J[i]));
    push!(u,H(J[i])); 
  end
  J=ideal(RP,v);

#Ideal J as an ideal of RP
  v=[];
  
  for i in 1:ngens(infedges)
    push!(v,H(infedges[i])); 
  end
  infedges=ideal(RP,v);

#=-------------------------------------------
Code chunk for required sorting of generators in ideal in order to align with the example on the paper
=#
#=
npars=0;
    p=gens(base_ring(S));
    para=[];
    
    for i in 1:length(p) 
    
        if p[i] in inverted_set(S)
            continue;
        else
            push!(para,p[i]); 
            npars=npars+1;
        end
          
    end
    startvars=npars+1;
    varRP=gens(RP);
#-----------------------------------
b=Vector{typeof(RP(1))}(undef,0);
c=Vector{typeof(RP(1))}(undef,0);
for i in 1:length(varRP) 
    if i>=startvars
        push!(b,gens(RP)[i]);
    else
        push!(c,gens(RP)[i]);
    end

end


    o=invlex(b)*lex(c);
    w=Vector{typeof(RP(1))}(undef,0);
    while length(u)!=0 
        k=1;
        

        for j in 2:length(u) 
            f=u[k];
            f1=leading_monomial(f,ordering=o);
            g=leading_monomial(u[j],ordering=o);
                    if cmp(o,f1,g)==-1 || cmp(o,f1,g)==0  
                        continue;
                    else
                        k=j;
                    
                        
                    end
            
        end
            push!(w,u[k]);
            deleteat!(u,[k]);       
    end
    w
    b=w[6];
    w[6]=w[4];
    w[4]=w[5];
    w[5]=b;
    J=ideal(RP,w);

#-------------------------------------------=#
  #=
Code chunk ended
=#
#reduce J w.r.t std basis of infedges
SB_infedges=standard_basis(infedges,ordering=invlex(RP));

N=Vector{typeof(RP(1))}(undef,0);
for i in 1:ngens(J)
    push!(N,RP(J[i]));     
end

M=Vector{typeof(RP(1))}(undef,0);
for i in 1:length(SB_infedges)
    push!(M,RP(SB_infedges[i]));     
end

N=reduce(N,M,ordering=invlex(RP),complete_reduction=true);

#rewite the ideal J as an ideal of S
T=hom(RP,S,gens(S));

v=[];
for i in 1:length(N)
push!(v,T(N[i])); 
end

J=ideal(S,v);
return J;
end


@doc raw"""
ISP(G::labeledgraph)

USAGE:  ISP(G);
ASSUME: G is a labeled graph.
RETURN: idal containing the irreducible scalar products(ISPs), that is, those scalar product which are not linearly dependent on the propagators.

#Examples
```julia
julia> G=simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);

julia> G=labelGraph(G,0);

julia> Gelim=eliminateVariables(G);

julia> ISP(Gelim)
Ideal generated by
  p[3]*q[1]
  p[1]*q[2]

```
"""
function ISP(G::labeledgraph)
    S=G.over;
    RP=G.overpoly;
    J=Feynman.propagators(G);
    el=G.elimvar;
    x=gens(S);
    v=[];
    for i in 1:length(G.edges)
        if length(G.edges[i])==1
            push!(v,G.labels[i]);
            
        end   
    end
    infedges=ideal(S,v);
    
    T=hom(S,RP,gens(RP));
    
    v=[];
    for i in 1:ngens(J)
    push!(v,T(J[i])); 
    end
    J=ideal(RP,v);
    
    v=[];
    for i in 1:ngens(infedges)
    push!(v,T(infedges[i])); 
    end
    infedges=ideal(RP,v);
    
    J=J+infedges^2;

    v=[];
    for i in 1:length(el)
    push!(v,T(el[i])); 
    end
    el=v;
    
    for i in 1:length(el)
    J=J+ideal(RP,el[i]);    
    end
    
    J=groebner_basis(J,ordering=invlex(RP));

    v=Vector{typeof(RP(1))}(undef,0);
    for i in 1:length(J)
    push!(v,T(J[i])); 
    end
    J=ideal(RP,v);

    w=Vector{typeof(RP(1))}(undef,0);   
    gens_RP=gens(RP);
    for i in 1:length(gens_RP) 
        Q,h=reduce_with_quotients(gens_RP[i]^2,v,ordering=invlex(RP),complete_reduction=true);
        if h!=0
            push!(w,h);
        end
        for j in i+1:length(gens_RP) 
        Q,h=reduce_with_quotients(gens_RP[i]*gens_RP[j],v,ordering=invlex(RP),complete_reduction=true);
        if h!=0
            push!(w,h);
        end
        end
    end
    
    K=ideal(RP,w);
    K=standard_basis(K,ordering=invlex(RP),complete_reduction=true);
    H=hom(RP,S,gens(S));
    u=Vector{typeof(S(1))}(undef,0);
    for i in 1:length(K) 
        push!(u,H(K[i]));
    end
    K=ideal(S,u);
    return K;
     
end


@doc raw"""
removeElimVars(G::labeledgraph)

USAGE:  removeElimVars(G);
ASSUME: G is a labled graph.
RETURN: Removes the variables from G.elimvars. This key is generated by the procedure eliminatedVariables.

#Examples
```julia
julia> G=simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);

julia> G=labelGraph(G,0);

julia> Gelim=eliminateVariables(G);

julia> G=removeElimVars(Gelim);
QQMPolyRingElem[p[1], p[2], p[3], p[4], q[1], q[2]]
julia> G.elimvar
Any[]```
"""
function removeElimVars(G::labeledgraph)
    R=G.over;
    RP=G.overpoly;
    el=G.elimvar;
    lb=G.labels;
    
    p=gens(base_ring(R));
    para=[];
    
    for i in 1:length(p) 
    
        if p[i] in inverted_set(R)
            continue;
        else
            push!(para,p[i]); 
        end          
    end
    
    iv=[];
    ip=[];
    ia=[];
        for i in 1:length(el)
            t=findall(x->x==el[i],p);
            push!(ia,t[1]);
            if el[i] in para
                push!(ip,t[1]);
            else
                push!(iv,t[1]);
            end
                   
        end
    
    R1=R;
    RP1=RP;

    if length(iv)!=0
      R1=Feynman.removeVariableLocal(R1,iv);
    end

    if length(ip)!=0
      R1=Feynman.removeParameter(R1,ip); 
    end

    if length(ia)!=0
      RP1=Feynman.removeVariable(RP1,ia);
        
    end
    #G.labels=lb;
    

    gR1=gens(R1);
  
    u=Vector{typeof(R1(1))}(undef,ngens(R));
    for i in 1:ngens(R)
      u[i]=R1(0); 
    end

    j=1;
    for i in 1:length(u) 
      if i in ia
      continue;
      else
        u[i]=gR1[j];
        j=j+1;
      end
    end
    T=hom(R,R1,u);
    u=Vector{typeof(R1(1))}(undef,0);
    for i in 1:length(lb) 
      push!(u,T(lb[i]));
    end
    G.labels=u;
    G.over=R1;
    G.overpoly=RP1;
    G.elimvar=[];
     return G; 
end


@doc raw"""
computeBaikovMatrix(G::simple_graphgraph)
USAGE:  computeBaikovMatrix(G);
ASSUME: G is a graph, or G is a labled graph where redundant variables have been eliminated by the procedure eliminateVariables, and deleted from the 
        ring by the procedure removeElimVars.
RETURN: a labeled graph G, where the computed Baikov matrix and the polynomial ring where baikovmatrix is defined are stored in G.baikovmatrix and G.baikovover respectively.

#Examples
```julia
julia> G=simple_graph([1,2,3,4,5,6,7],[(6,1),(6,4),(1,2),(3,7),(4,3),(2,7),(5,6),(7,5),1,2,3,4,5]);

julia> G=Feynman.labelGraph(G,0);

julia> G=Feynman.eliminateVariables(G);

julia> G=Feynman.removeElimVars(G);
QQMPolyRingElem[p[1], p[2], p[3], p[4], p[5], q[1], q[2]]
julia> G=Feynman.computeBaikovMatrix(G);
labels used for Gram matrix of external loop momenta:
["p[1]*p[2] => 1//2*t[1]"]
["p[1]*p[3] => 1//2*t[2]"]
["p[2]*p[3] => 1//2*t[3]"]
["p[1]*p[4] => 1//2*t[4]"]
["p[2]*p[4] => 1//2*t[5]"]
["p[3]*p[4] => -1//2*t[1] - 1//2*t[2] - 1//2*t[3] - 1//2*t[4] - 1//2*t[5]"]
Assignment of Baikov variables (Z_i) are:
["z[1] => p[3]*q[1]"]
["z[2] => p[4]*q[1]"]
["z[3] => q[1]^2"]
["z[4] => -2*p[1]*q[1] + q[1]^2"]
["z[5] => 2*p[1]*p[2] - 2*p[1]*q[1] - 2*p[2]*q[1] + q[1]^2"]
["z[6] => p[1]*q[2]"]
["z[7] => q[2]^2"]
["z[8] => -2*p[1]*p[2] - 2*p[1]*p[3] - 2*p[1]*p[4] - 2*p[2]*p[3] - 2*p[2]*p[4] - 2*p[3]*q[2] - 2*p[4]*q[2] + q[2]^2"]
["z[9] => -2*p[4]*q[2] + q[2]^2"]
["z[10] => q[1]^2 + 2*q[1]*q[2] + q[2]^2"]
["z[11] => -2*p[1]*q[1] - 2*p[1]*q[2] - 2*p[2]*q[1] - 2*p[2]*q[2] - 2*p[3]*q[1] - 2*p[3]*q[2] - 2*p[4]*q[1] - 2*p[4]*q[2] + q[1]^2 + 2*q[1]*q[2] + q[2]^2"]

julia> G.baikovmatrix
6×6 Matrix{RingElem}:
 0                      …  z[6]
 1//2*t[1]                 1//2*t[2] + 1//2*t[3] + 1//2*t[4] + 1//2*t[5] - z[1] - z[2] - 1//2*z[3] + 1//2*z[5] - z[6] - 1//2*z[7] + 1//2*z[8] + 1//2*z[10] - 1//2*z[11]
 1//2*t[2]                 -1//2*t[1] - 1//2*t[2] - 1//2*t[3] - 1//2*t[4] - 1//2*t[5] - 1//2*z[8] + 1//2*z[9]
 1//2*t[3]                 1//2*z[7] - 1//2*z[9]
 1//2*z[3] - 1//2*z[4]     -1//2*z[3] - 1//2*z[7] + 1//2*z[10]
 z[6]                   …  z[7]```
"""

function computeBaikovMatrix(G::simple_graph)
    G=labelGraph(G,0);
    G=eliminateVariables(G);
    G=removeElimVars(G);
return computeBaikovMatrix(G);
end

function computeBaikovMatrix(G::labeledgraph)
    if typeof(G)=="simple_graph"
        lG=labelGraph(G,0);
        G1=eliminateVariables(lG);
        G2=removeElimVars(G1);
        return computeBaikovMatrix(G2);
        
    end
    
    #if typeof(G)!="labeledgraph"
    ##    throw(error("expected a graph or labeledgraph"));
    # end

    R=G.over;
    RP=G.overpoly;
    P=Feynman.propagators(G);
    ngens(P)
    I=Feynman.ISP(G);
    PI=P+I;
    idx=0;
#-------------------------------------------------------------
    H=hom(R,RP,gens(RP));
    T=hom(RP,R,gens(R));
#--------------------------------------------------------------
    a=Vector{typeof(RP(1))}(undef,0);
    for i in 1:ngens(P) 
        push!(a,H(P[i]));
    end
    P=ideal(RP,a);
    v=Vector{typeof(RP(1))}(undef,0);
   
    for i in 1:ngens(P) 
        push!(v,P[i]);
       
    end
    
#-------------------------------------ordering elements
    npars=0;
    p=gens(base_ring(R));
    para=[];
    
    for i in 1:length(p) 
    
        if p[i] in inverted_set(R)
            continue;
        else
            push!(para,p[i]); 
            npars=npars+1;
        end
          
    end
    startvars=npars+1;
    varRP=gens(RP);
    #-----------------------------------------

    gram1=Vector{typeof(varRP[1])}(undef,0);

    #H=hom(R,RP,gens(RP));
    idx=0;
    for i in 1:length(varRP)
        for j in 1:length(varRP) 
        if (i>=startvars) || (j>=startvars)
            push!(gram1,H(varRP[i]*varRP[j]));
            idx=idx+1;
        else
            push!(gram1,RP(0));
        end 
        end
    
    end
    #gram1

#-----print gram vector as a matrix
    W=Matrix{typeof(RP(1))}(undef,ngens(RP),ngens(RP));    
    for j in 1:length(gram1)
            i=div(j,ngens(RP))+1;
            k=mod(j,ngens(RP));
            #print(gram1[j])
            if k==0
                W[i-1,ngens(RP)]=gram1[j];  
                #print("(",i-1,",",ngens(RP),")=",gram1[j])
            else
                W[i,k]=gram1[j];   
               # print("(",i,",",k,")=",gram1[j])
            end     
    end

    
#---Adding terms independent of loop momenta
    for i in 1:startvars-1 
        for j in i+1:startvars-1 
            PI=PI+ideal(R,[varRP[i]*varRP[j]]);
        end
    
    end
    #ngens(PI)
    #PI
#---Writing PI as an idal of RP

    a=Vector{typeof(RP(1))}(undef,0);
    for i in 1:ngens(PI) 
        push!(a,H(PI[i]));
    end
    PI=ideal(RP,a);
   
#---------------------------------

    v=Vector{typeof(RP(1))}(undef,0);
    for i in 1:ngens(PI) 
        push!(v,H(PI[i]));
    end

#---------------Ordering terms in IP appropreately for further calculations--------------------------
    o=invlex(RP);
    w=Vector{typeof(RP(1))}(undef,0);
    while length(v)!=0 
        k=1;
        for j in 2:length(v) 
            f=v[k];
            f1=leading_monomial(f,ordering=o);
            g=leading_monomial(v[j],ordering=o);
                    if cmp(o,f1,g)==1   
                        k=j;                            
                    end
            
        end
            push!(w,v[k]);
            deleteat!(v,[k]);       
    end
    v=w;
    PI=ideal(RP,v);

#----------------------------------
    m=length(para)
    m2=Int(m*(m-1)/2);  
    if m2==0
        mt=1;
    else
        mt=m2-1;
    end
    n=ngens(PI)-m2;

    #ngens(PI)
    Z=Feynman.makePoly(mt,n);
    t=gens(Z);
    t=deleteat!(t,mt+1);
    z=Vector{typeof(t[1])}(undef,0);
    for i in 1:n
        push!(z,t[mt+i]);         
    end
    
    B=zeros(Z,length(varRP),length(varRP));
    pq=Vector{typeof(t[1])}(undef,0);
    sumt=Z(0);

    idx=1;
    for i in 1:m
        for j in i+1:m 
            if idx<=mt
                B[i,j]=Z(1//2*t[idx]);
                sumt=sumt+Z(1//2*t[idx]);
                push!(pq,Z(1//2*t[idx]));
            else
                B[i,j]=-sumt;
                push!(pq,-sumt);
            end 
            B[j,i]=B[i,j];
            idx=idx+1;
        end 
    
    end


    zvar=Matrix{typeof(t[1])}(undef,1,n+m2);
    
    for i in 1:m2
        zvar[1,i]=pq[i];     
    end
    for i in 1:n    
        zvar[1,i+m2]=z[i]; 
    end

    #-----------------Printing of Baikov variable assigment
    ev=[];
    println("labels used for Gram matrix of external loop momenta:");
    for i in 1:m2 
        t=string(string(gens(PI)[i])," => ",string(zvar[i]));
        push!(ev,t);
        printNet(ev);
        ev=[];
    end
  
    println("Assignment of Baikov variables (Z_i) are:");
    ev=[];
    for i in 1:n 
        t=string(string(zvar[i+m2])," => ",string(gens(PI)[i+m2]));
        push!(ev,t);
        printNet(ev);
        ev=[];
    end
  

#---------------------------------------------------------
    X=ideal(RP,gram1);
#    is_subset(X,PI)
    H,I=groebner_basis_with_transformation_matrix(PI,ordering=invlex(RP));
    T,J=groebner_basis_with_transformation_matrix(X,ordering=invlex(RP));
    I1=Matrix{typeof(Z(1))}(undef,ngens(PI),length(H));
    I2=Matrix{Int}(undef,ngens(PI),length(H));
    I3=Matrix{typeof(RP(1))}(undef,ngens(PI),length(H));
    J1=Matrix{typeof(Z(1))}(undef,ngens(X),length(T));
#   gens(PI)*I==gens(H)
    #gens(PI)
    #gens(PI)
    y=Matrix{typeof(RP(1))}(undef,1,ngens(PI));
    for i in 1:ngens(PI) 
        y[1,i]=PI[i];
    end
    
    for i in 1:ngens(PI) 
        for j in 1:length(H)
            I1[i,j]=Z(constant_coefficient(I[i,j]));
            I2[i,j]=Int(constant_coefficient(I[i,j]));
        end
    end

    for i in 1:ngens(X)
        for j in 1:length(T) 
           J1[i,j]=Z(constant_coefficient(RP(J[i,j])));
           
        end
        
    end

    v=Vector{typeof(RP(1))}(undef,0);
    w=Vector{typeof(RP(1))}(undef,0);
    V=Matrix{typeof(RP(1))}(undef,1,length(H));
    for i in 1:length(H) 
        push!(v,H[i]);
        V[1,i]=H[i];
    end
    for i in 1:length(T) 
        push!(w,T[i]);
    end

   # is_subset(ideal(RP,w),ideal(RP,v))

    A=Matrix{typeof(Z(1))}(undef,length(v),length(gram1));
    D=Matrix{typeof(Z(1))}(undef,length(v),length(gram1));

    for i in 1:length(gram1)
         u=reduce_with_quotients(gram1[i],v,ordering=invlex(RP));
         #print(u)
            for j in 1:length(u[1]) 
                if u[1][j]==0
                    A[j,i]=Z(0);
                    D[j,i]=RP(0)
                else
                    A[j,i]=Z(constant_coefficient(u[1][j]));
                    D[j,i]=RP(u[1][j]);

                end
            end 
    end

    #transpose(V*D)
    f=Matrix{typeof(RP(1))}(undef,nrows(I),ncols(I));
    for i in 1:nrows(I)
        for j in 1:ncols(I)
            f[i,j]=RP(I[i,j]);
        end 
        
    end
    #transpose(y*f*D)

    #Bentries=zvar*transpose(I1)*A
    Bentries=zvar*I1*A;
    B1=Matrix{typeof(Z(1))}(undef,ngens(RP),ngens(RP));    
    for j in 1:length(Bentries)
            i=div(j,ngens(RP))+1;
            k=mod(j,ngens(RP));
            if k==0
                B1[i-1,ngens(RP)]=Bentries[1,j];  
            else
                B1[i,k]=Bentries[1,j];   
            end     
    end
    
    B=B1+B;
    G.baikovover=Z;
    G.baikovmatrix=B;
    return G;


    
end

@doc raw"""
makePoly(n::Int,m::Int)

USAGE:  makePoly(m,n);
ASSUME: m and n are positve integers.
RETURN: A polynomial ring with vatiables t[1],...,t[n],z[1],...,z[m] over QQ. 
"""
function makePoly(n::Int,m::Int)
    v=Vector{String}(undef,0);
    for i in 1:n 
        push!(v,"t[$i]");  
    end
        push!(v,"D");
    for i in 1:m 
        push!(v,"z[$i]");  
    end
    #Z,t,z=polynomial_ring(QQ,"t"=>(1:n),"z"=>(1:m));
    Z,t=polynomial_ring(QQ,v);
    return Z;
end

@doc raw"""
computeIBP(G::labeledgraph,Nu::Vector{Int64},cutDeg::Int)

USAGE:  computeIBP(G,ν,d);
ASSUME: G is a labeled graph, d is a positive integer and ν is vector of integers correspond to the parent diagram of the integral.
RETURN: A set of simplified IBP identities without double propagators (without performing trimming) . 
"""
function computeIBP(G::simple_graph,Nu::Vector{Int64},cutDeg::Int,showGens::Bool)
    G=Feynman.labelGraph(G,0);
    G=Feynman.eliminateVariables(G);
    G=Feynman.removeElimVars(G);
    G=Feynman.computeBaikovMatrix(G);
    return computeIBP(G,Nu,cutDeg,showGens);
end

function computeIBP(G::simple_graph,Nu::Vector{Int64},cutDeg::Int)
    G=Feynman.labelGraph(G,0);
    G=Feynman.eliminateVariables(G);
    G=Feynman.removeElimVars(G);
    G=Feynman.computeBaikovMatrix(G);
    return computeIBP(G,Nu,cutDeg,true);
end

function computeIBP(G::labeledgraph,Nu::Vector{Int64},cutDeg::Int)
    return computeIBP(G,Nu,cutDeg,true);
end

function computeIBP(G::labeledgraph,Nu::Vector{Int64},cutDeg::Int,showGens::Bool)

    RZ=G.baikovover;
    R=G.over;
    gens_RZ=gens(RZ);
    B=G.baikovmatrix;
#-----------------------------------------------------------------------------------#
    S=matrix_space(RZ,size(G.baikovmatrix)[1],size(G.baikovmatrix)[2]);
    C=S(G.baikovmatrix);
    f=det(C);  
#-----------------------------------------------------------------------------------#
#Compute E and Count L
   npars=0;
   p=gens(base_ring(R));
   para=[];
   for i in 1:length(p) 
       if p[i] in inverted_set(R)
           continue;
       else
           push!(para,p[i]); 
           npars=npars+1;
       end
         
   end
   m=length(para);
   E=m;
   m2=Int(m*(m-1)/2);
   mt=0;  
   if m2==0
       mt=1;
   else
       mt=m2-1;
   end
   L=length(p)-E;        
   #----------Getting t vector and z vector-------------
   var_t=Vector{typeof(RZ(1))}(undef,0);
   var_z=Vector{typeof(RZ(1))}(undef,0);
   for i in 1:mt
       push!(var_t,gens_RZ[i]);
   end
   for i in mt+2:length(gens_RZ) #changed 
       push!(var_z,gens_RZ[i]);
   end

    if length(Nu)!=length(var_z)
        error("length of the vector nu must equal to number of Baikov variables");
    end  

#-------------Computing generators of M1--------------------------
   t=computeM1(G);
   n=L*E+Int(L*(L+1)/2);
#=Test for whether the computed generators are correct

    for j in 1:length(t) 
        g=RZ(0);
    for i in 1:(length(t[1])-1) 
        g=g+t[j][i]*derivative(f,length(var_t)+1+i);
    end
        println(g+t[j][length(t[1])]*f)
    end
=#
    m=length(var_z);


#-----------------------------------------------------------------------------------#
# --------------compute module intersection and  Groebner basis---------------------#
#-----------------------------------------------------------------------------------#

##------Defining the polynomial ring RZ in Singular----------------------------------#
    v1=Vector{String}(undef,0);
    v2=Vector{String}(undef,0);

    for i in 1:length(var_t) 
        push!(v1,"t[$i]");  
    end
    push!(v1,"D");

    for i in 1:length(var_z) 
        push!(v2,String("z[$i]"));  
    end

    v=convert(AbstractArray,vcat(v1,v2));
    T,var=Singular.polynomial_ring(Singular.QQ,v,ordering=Singular.ordering_ls(length(v1))*Singular.ordering_dp(),degree_bound=cutDeg);
##-------Build the module M1----------------------------------------------------------#

    gens_M1=[];
    for i in 1: length(t)#size(A)[2]
        v=Vector{typeof(T(1))}(undef,0);
   #u=transpose(A[:,i]);
        u=t[i];
    for j in 1:length(u)-1 
        push!(v,T(u[j]))  
    end
        push!(gens_M1,Singular.vector(T,v...))
   
    end
    M1=Singular.Module(T,gens_M1...);

##--------Build the module M2-----------------------------------------------------------#
    gens_M2=[];
    for j in 1:length(t[1])-1
    v=Vector{typeof(T(1))}(undef,0);;
    for i in 1:length(t[1])-1 
        if i==j && j<=length(t[1])-1
             if Nu[j]>0 
                 push!(v,T(var_z[j]));       
            else #j<=length(t[1]) #j<=length(t[1])-1
                push!(v,T(1));
            end
    
        else
                push!(v,T(0));
        end          
    end 
    push!(gens_M2,Singular.vector(T,v...));
    end

    M2=Singular.Module(T,gens_M2...);
  
 #=   if showGens
 #   println("Generators of M1 are:");
    printNet(gens_M1);
    println("Generators of M2 are:");
    printNet(gens_M2);
   end
   =# 
##-------------------Compute Module intersection------------------------------------------# 
    G=Singular.intersection(M1,M2);
    G=Singular.std(G,complete_reduction=true);

#-----------------------------------------------------------------------------------#
# ------compute IBP identities correspond to generators of the Groebner basis--------#
#------------------------------------------------------------------------------------#

##------Convert generators of Groebner basis (in Singular) to generators in Oscar-----
    vecG=[];
    for i in 1:Singular.number_of_generators(G) 
        u=Singular.Base.Array(G[i]);
        v=[];
        for j in 1:length(u) 
            push!(v,RZ(u[j]));
        end
   
       push!(vecG,v);
    end

#--------------------------------------------------
#--------------------------------------------------
    for j in 1:length(vecG) 
        w=vecG[j]
        g=RZ(0);
        for i in 1:length(w) 
            g=g+w[i]*derivative(f,length(var_t)+1+i);
        end

        if g==0
            bj=0;
            hj=0;
            push!(vecG[j],bj[1]);
        else
       #e=j-invP;
            bj, hj = reduce_with_quotients(g, [f], ordering = invlex(RZ));
            if hj==0
                push!(vecG[j],bj[1]);
            else
            deleteat!(VecG,j);
            end
        end

    end

#=------------------test again for correct calculation
    for i in 1:length(vecG) 
        w=vecG[i]
        g=RZ(0);
        for i in 1:length(w)-1 
            g=g+w[i]*derivative(f,length(var_t)+1+i);
        end
        println(g-w[length(w)]*f)

    end  
=#
##---------------Computation of IBP identities------------------------------------------------
    set_IBP=[];
    for i in 1:length(vecG) 
        u=vecG[i];
        v=[];
        proz=RZ(1);        
        for j in 1:m
            proz=proz*var_z[j];
        end

        for j in 1:length(u)-1 
       
               if u[j]==0
                   bj=0;
                     hj=0;
                else
              #e=j-invP;
                    g=proz*u[j];
                    bj, hj = reduce_with_quotients(g, [var_z[j]], ordering = invlex(RZ));
                end
                push!(v,bj[1]);      
        end
   #Compute Baikov identity associated to generator vecG[i]
   ##---1---compute m1= (\sum_{i=1}^{m}\frac{\partial a_i}{\partial z_i}
        m1=RZ(0);
        for j in 1:length(u)-1 
            m1=m1+derivative(u[j],length(var_t)+1+j);
        end
   ##---2---compute m3=\sum_{i=1}^{m} \frac{\nu_i a_i}{z_i}
        m2=RZ(0);
        for j in 1:length(u)-1 
            m2=m2+RZ(Nu[j])*v[j];
        end
   M=RZ((m1-(gens_RZ[length(var_t)+1]-L-E-1)*RZ(1//2)*u[m+1])*proz-m2);
        single_IBP=[];
        single_IBP_i=[];
        single_IBP_c=[];
        for j in 1:length(M) 
            degT=degrees(monomial(M,j)); #degree sequence of (t1,...,t_k,D,z_1,...,z_m)
            deg_t=first(degT,length(var_t)+1);
            deg_z=last(degT,length(var_z));
            cTerm=RZ(1);
            for l in 1:length(var_t)+1
          #cTerm=cTerm*var_t[l]^deg_t[l];
                cTerm=cTerm*gens(RZ)[l]^deg_t[l];
            end
            cTerm=cTerm*coeff(M,j);
            iTerm=deg_z-Nu;
       #Since we multiplied the expression by prod_z, 
            for l in 1:length(iTerm) 
                iTerm[l]=iTerm[l]-1;
            end
            push!(single_IBP,[cTerm,iTerm]);
            push!(single_IBP_c,cTerm);
            push!(single_IBP_i,iTerm);
        end
        if single_IBP!=[]
            #Need to simplify the IBP
            IBP=[];
            while length(single_IBP_i)!=0
                w=single_IBP_i[1];
                ind=findall(x->x==w,single_IBP_i);
                coef=RZ(0);
                for l in 1:length(ind) 
                    coef=coef+single_IBP_c[ind[l]];
                end
                push!(IBP,[coef,w]);
                deleteat!(single_IBP_i,ind);
                deleteat!(single_IBP_c,ind);
            end
            push!(set_IBP,IBP);   
        end
    end
 
return IBP(Nu,RZ,length(v1),set_IBP);
end

@doc raw"""
printIBP(set_IBP::Vector{Vector{}},n::Int64)
USAGE:  printIBP(single_IBP,n); 
"""
function printIBP(set_IBP::Vector,n::Int64)
    
    if n>length(set_IBP)
        error("N must be less or equal to number of IBP identities");
    end
    println("First ",n," IBP identities associated to G  (Total number of relations=",length(set_IBP),"):")
   
    for j in 1:n
        ve=set_IBP[j];
        ev=[];
    for i in 1:length(ve) 
        t1=join(repr.(ve[i][2]),",");
        t=string(string("("),string(ve[i][1]),string(")"),string("I("),t1,string(")"));
        push!(ev,t);
    end
    println('0','=',join(ev,"+")); 
     println();   
    end
end

function computeM1(G::labeledgraph)
    RZ=G.baikovover;
    R=G.over;
    gens_RZ=gens(RZ);
    B=G.baikovmatrix;
    
    S=matrix_space(RZ,size(G.baikovmatrix)[1],size(G.baikovmatrix)[2]);
    C=S(G.baikovmatrix);
    f=det(C);

    npars=0;
    p=gens(base_ring(R));
    para=[];
    for i in 1:length(p) 
    
        if p[i] in inverted_set(R)
            continue;
        else
            push!(para,p[i]); 
            npars=npars+1;
        end
          
    end
    m=length(para);
    E=m;
    m2=Int(m*(m-1)/2);
    mt=0;  
    if m2==0
        mt=1;
    else
        mt=m2-1;
    end
    L=length(p)-E;        
    #----------Getting t vector and z vector-------------
    var_t=Vector{typeof(RZ(1))}(undef,0);
    var_z=Vector{typeof(RZ(1))}(undef,0);
    for i in 1:mt
        push!(var_t,gens_RZ[i]);
    end
    for i in mt+2:length(gens_RZ) #changed 
        push!(var_z,gens_RZ[i]);
    end

    t=[];
   n=L*E+Int(L*(L+1)/2);


##------------Compute the partial derivative correctly z_alpha/x(i,j)-------------------------------
    W=Matrix{typeof(RZ(1))}(undef,length(var_z),length(var_z));    
    r=0;
    ij=[];
    for i in E+1:E+L 
        for j in 1:E+L-1 
            r=r+1;
            push!(ij,[i,j]);
            for l in 1:length(var_z) 
                W[r,l]=derivative(B[i,j],var_z[l]);
            end
        end
    end
    push!(ij,[E+L,E+L])
    for l in 1:length(var_z) 
        W[length(var_z),l]=derivative(B[E+L,E+L],var_z[l]);

    end
    S=matrix_space(RZ,size(W)[1],size(W)[2]);
    C=S(W);
    C=inv(C);
    for i in E+1:E+L
        for j in 1:E+L
            v=Vector{typeof(RZ(1))}(undef,0);
            for l in 1:n #here l is alpha
                a=RZ(0);
                for k in 1:E+L 
                if i==k
                    d=RZ(2);
                else
                     d=RZ(1);
                end
           
                if findall(x->x==[i,k],ij)!=[]
                    c=C[l,findall(x->x==[i,k],ij)[1]];
                elseif findall(x->x==[k,i],ij)!=[]
                        c=C[l,findall(x->x==[k,i],ij)[1]];
                else
                        c=RZ(0);   
                end
           
                a=a+d*c*B[j,k];
                end   
                push!(v,a);
        
            end
                if i==j
                    b=RZ(-2);
                else
                    b=RZ(0);
                end
                push!(v,b);
                push!(t,v);
        #   v=Vector{typeof(RZ(1))}(undef,0);
        end 
   
    end
    return t;
end

function enumerateTerms(sets)
    if isempty(sets)
        return [()]
    end
    first_set=sets[1];
    rest_sets=sets[2:end];
    rest_product=enumerateTerms(rest_sets);
    result=[];
    for elim in first_set
        for combo in rest_product 
            push!(result,(elim,combo...));
        end
    end
        return result;

end
function Heviside(x)
    if x>=0
        return 1;
    else
        return 0;
        
    end
end

function power_monomials(nvar::Int64,deg::Int64)
    if nvar==1
        return [[deg]];        
    end
    vecPower=[];
    for i in 0:deg 
        for j in power_monomials(nvar-1,deg-i) 
            push!(vecPower,[i;j])
        end
    end
    return vecPower
end

function getIndDen(R::Ring,common_factors::Vector)
    n=nvars(R);
    vars=gens(R);
  #  common_factors=Vector{typeof(R(1))}(undef,0);
  #  for i in 1:length(common_den) 
  #      push!(common_factors,R(common_den[i][1]));
  #  end

    X=matrix_space(QQ,length(common_factors),length(vars));
    x=Matrix(zero_matrix(QQ,length(common_factors),length(vars)));
    for i in 1:length(common_factors) 
        f=common_factors[i];
        while f!=0
            for j in 1:length(vars) 
                if leading_monomial(f)==vars[j]
                    x[i,j]=leading_coefficient(f);#R1(leading_coefficient(f));
                end
            end   
            f=tail(f); 
        end
    end
    x=X(x);
    r, P, L, U = lu(x);
    ind_den=[];
    ind_w=[];
    if r<=length(vars)
        for i in 1:r 
            push!(ind_den,common_factors[P[i]]);
            push!(ind_w,P[i]);
        end
    end
    if r<n 
        x=x[1:r,:];
        nl, N = nullspace(x);
        for j in 1:(n-r)#(size(N))[2]   #applying rank-nullity theorem 
            f=0;
            for k in 1:length(vars) 
               f=f+N[k,j]*vars[k]; 
            end
            push!(ind_den,f);
        end
    end
 
    if r>n
        error("Something is wrong in the denominators provided");        
    end
    return(ind_den,ind_w);
end

function generating_probes(R::Ring,P::Int64,common_factors::Vector,w1::Vector,t::Int64)
 CF=common_factors;
    p=P;
    R1,_= QadicField(ZZ(p), 1,30);
    Qx, x = QQ["x"];
F=GF(P);
Qx, x = QQ["x"];
vars=gens(R);
common_factors,i_w=getIndDen(R,common_factors);
#---------------------------------------------------------
    w=Vector{typeof(1)}(undef,length(vars));
    for i in 1: length(vars)
        w[i]=w1[i]#rand(0:8)
    end
    for i in 1:length(i_w)
        w[i]= w1[i_w[i]];
    end
    #undef_indices=findall(x->x===nothing || x ===missing,w);
    #for i in 1:length(undef_indices) 
    #    w[undef_indices[i]]=rand(0:4);
    #end
    X=matrix_space(QQ,length(common_factors),length(vars));
    x=Matrix(zero_matrix(QQ,length(common_factors),length(vars)))
    
    for i in 1:length(common_factors) 
        f=common_factors[i];
        while f!=0
            for j in 1:length(vars) 
                if leading_monomial(f)==vars[j]
                    x[i,j]=leading_coefficient(f)#R1(leading_coefficient(f));                   
                end
            end   
            f=tail(f); 
        end
    end
    x=X(x);
    C=matrix_space(QQ,length(common_factors),1);
    c=Matrix{typeof(QQ(1))}(undef,length(common_factors),1)
    d=Matrix{typeof(QQ(1))}(undef,length(common_factors),1)
    #P=Int64(101)

    for i in 1:length(common_factors) 
        l1=rand(1:(t*1000));
        l=rand(1:t)//l1;
        x1=R1(l)+O(R1,P^w[i]); 
        p= l-lift(Qx, R1(x1));
        l2=coeff(p,0);
        l3=coeff(lift(Qx, R1(x1)),0);
        while denominator(l2)%P==0 #&& denominator(l3)%P==0 
            l1=rand(1:(t*1000));
            l=rand(1:t)//l1;
            x1=R1(l)+O(R1,P^w[i]); 
            p= l-lift(Qx, R1(x1));
            l2=coeff(p,0);
            l3=coeff(lift(Qx, R1(x1)),0);
        end
    #l=rand(1:t)//l1;
    #x1=R1(l)+O(R1,P^w[i]); 
    #p= l-lift(Qx, R1(x1));
    #c[i,1]=coeff(p,0);
        c[i,1]=l2;
        d[i,1]=coeff(lift(Qx, R1(x1)),0);
    end
    c=C(c)
    d=C(d)
    prob1=Matrix(inv(x)*c)
    prob2=Matrix(inv(x)*d)
    e1=prob1[1,1]
    e2=prob1[2,1]
    e3=prob1[3,1]

#=    CF=common_factors;
    p=P;
    R1,_= QadicField(ZZ(p), 1,30);
    Qx, x = QQ["x"];
    common_factors,i_w=getIndDen(R,common_factors);
    vars=gens(R);
    w=Vector{typeof(1)}(undef,length(vars));
    for i in 1:length(i_w)
        w[i]= w1[i_w[i]];
        #if i<=length(i_w)
         #   push!(w,w1[i_w[i]]);
        #else
         #   push!(w,rand(0:4));
        #end
    end
    
    undef_indices=findall(x->x===nothing || x ===missing,w);
    
    for i in 1:length(undef_indices) 
        w[undef_indices[i]]=rand(0:4);
    end

    X=matrix_space(QQ,length(common_factors),length(vars));
    x=Matrix(zero_matrix(QQ,length(common_factors),length(vars)))
    
    for i in 1:length(common_factors) 
        f=common_factors[i];
        while f!=0
            for j in 1:length(vars) 
                if leading_monomial(f)==vars[j]
                    x[i,j]=leading_coefficient(f)#R1(leading_coefficient(f));                   
                end
            end   
            f=tail(f); 
        end
    end
    x=X(x);
    C=matrix_space(QQ,length(common_factors),1);
    c=Matrix{typeof(QQ(1))}(undef,length(common_factors),1)
    d=Matrix{typeof(QQ(1))}(undef,length(common_factors),1)
    #P=Int64(101)

    for i in 1:length(common_factors) 
        l1=rand(1:(t*1000));
        l=rand(1:t)//l1;
        x1=R1(l)+O(R1,P^w[i]); 
        p= l-lift(Qx, R1(x1));
        l2=coeff(p,0);
        l3=coeff(lift(Qx, R1(x1)),0);
        while denominator(l2)%P==0 #&& denominator(l3)%P==0 
            l1=rand(1:(t*1000));
            l=rand(1:t)//l1;
            x1=R1(l)+O(R1,P^w[i]); 
            p= l-lift(Qx, R1(x1));
            l2=coeff(p,0);
            l3=coeff(lift(Qx, R1(x1)),0);
        end
    #l=rand(1:t)//l1;
    #x1=R1(l)+O(R1,P^w[i]); 
    #p= l-lift(Qx, R1(x1));
    #c[i,1]=coeff(p,0);
        c[i,1]=l2;
        d[i,1]=coeff(lift(Qx, R1(x1)),0);
    end
    c=C(c)
    d=C(d)
    prob1=Matrix(inv(x)*c)
    prob2=Matrix(inv(x)*d)
    =#
    w_l=[]
    f=prob1[1,1]
    g=prob1[2,1]
    h=prob1[3,1]
    for k in 1:length(CF) 
        push!(w_l,valuation(R1(evaluate(CF[k],[f,g,h]))))
    end

    if w_l==w1
        return prob1,prob2 ;
    else
        return generating_probes(R,P,CF,w1,t);   
    end

#return prob1,prob2;
end

function oneRR(r::Int64,m::Int64)
u=[1,0,m];
v=[0,1,r];
while v[3]>=sqrt(m/2)
    q=floor(u[3]//v[3]);
    r=u-q*v;
    u=v;
    v=r;
end

if abs(v[2])>sqrt(m/2)
    throw(error("given m and r doesn't satisfy the condition in Wang algorithm"));
end
return([v[3],v[2]]);
end

function manyRR(ri::Vector,mi::Vector)
    m=1;
for i in 1:length(mi) 
    m=m*mi[i];
end
c=crt(ri,mi);
return(oneRR(c,m));
end

function pickDen(R::Ring,R1::Field,common_factors::Vector,w::Vector,candidate_den::Vector,nProbes::Int )
    i=0;
    d=0;
    m=0;
    pr=[];
    pr_2=[];
    p=Int(prime(R1));
   while i<nProbes 
       prob1,prob2=generating_probes(R,p,common_factors,w,rand(10000:1000000));
       p1=prob1[1,1]
       q1=prob1[2,1]
       r1=prob1[3,1]
       p2=prob2[1,1]
       q2=prob2[2,1]
       r2=prob2[3,1]
       pw=[];
       for j in 1:length(candidate_den) 
           a=evaluate(candidate_den[j],[p1,q1,r1]) 
           push!(pw,valuation(R1(a)))
       end
       d=candidate_den[findall(x->x==pw[argmax(pw)],pw)[1]];
       if length(findall(x->x==pw[argmax(pw)],pw))<=length(gens(R)) #&& valuation(R1(evaluate(d,[p1,q1,r1])))!=0#<=length(gens(R)) #&& evaluate(f,[p2,q2,r2])!=0
           
           push!(pr,[p1,q1,r1]);
           push!(pr_2,[p2,q2,r2]);
           d=candidate_den[findall(x->x==pw[argmax(pw)],pw)[1]];
           #print(d)
           m=valuation(R1(evaluate(d,[p1,q1,r1])))
               i=i+1;
       end
       
   end
   return([pr,pr_2,d,m]); #this depend on candidate den need to change when  candidate den is changes
end


function evaluateBB(p::Vector,v::Vector)
    a=evaluate(p[1],v); 
    b=evaluate(p[2],v);
    return(a/b);
end

function pickWeights(R::Ring,R1::Field,BB::Vector,common_den::Vector,candidate_w::Vector,exp_vec::Vector,computed_den::Vector,computed_num)
    #extracting terms
    p=Int(prime(R1));
    common_factors=Vector{typeof(R(1))}(undef,0);
    candidate_den=[];
    for i in 1:length(common_den) 
        push!(common_factors,R(common_den[i][1]))
    end
    q8,q9=getIndDen(R,common_factors);
    for i in 1:length(exp_vec) 
        c=1;
        for j in 1:length(common_factors) 
            c=c*common_factors[j]^exp_vec[i][j];
        end 
        push!(candidate_den,c);
    end
    trimmed_w=[];
    d_i=[];
    W_i=[];
    f=BB[2];
    for i in 1:length(candidate_w) 
        w=candidate_w[i];
        W=[];#these are weights of candidate denominators
        prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000));
        p1=prob1[1,1];
        q1=prob1[2,1];
        r1=prob1[3,1];
        for i in 1:length(exp_vec) 
        #toDo write a function to obtain the correspoding denominator 
            e=0;
            for j in 1:length(exp_vec[i]) 
             #   e=e+valuation(R1(evaluate(common_den[j][1],[p1,q1,r1])))*exp_vec[i][j];
                e=e+w[j]*exp_vec[i][j]
            end
        #push!(W,w[1]*exp_vec[i][1]+w[2]*exp_vec[i][2]+w[3]*exp_vec[i][3]);
            push!(W,e);
        end
        
        a=evaluateBB(BB,[p1,q1,r1]); 
        b=evaluate(f,[p1,q1,r1]);
    
        while valuation(R1(b))==0
            prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000));
            p1=prob1[1,1];
            q1=prob1[2,1];
            r1=prob1[3,1];
            a=evaluateBB(BB,[p1,q1,r1]);
            b=evaluate(f,[p1,q1,r1]);
        end
    
        constructed_term=0;
        if length(computed_den)!=0
            for k in 1:length(computed_den) 
                constructed_term=constructed_term+R1(evaluate(computed_num[k],[p1,q1,r1])/evaluate(computed_den[k],[p1,q1,r1]) )
            end    
        end
    
        R_w=valuation(R1(a-constructed_term))
        d_i=[];
        W_i=[];
        for i in 1:length(exp_vec) 
            lhs=0;
            v1=[];
            v2=[];
            for j in 1:length(w) 
                lhs=lhs-exp_vec[i][j]*w[j];
                # lhs=lhs-valuation(R1(evaluate(common_den[j][1],[p1,q1,r1])))*exp_vec[i][j];
                if exp_vec[i][j]!=0
                    push!(v1,common_den[j][1]);
                end
                if exp_vec[i][j]!=0 || w[j]!=0
                #if exp_vec[i][j]!=0 || valuation(R1(evaluate(common_den[j][1],[p1,q1,r1])))!=0
                    push!(v2,common_den[j][1]);
                end
            end    
            I1=ideal(R,v1);
            I2=ideal(R,v2);
            println(lhs);
            if R_w>lhs && I1==I2
                continue;
            else
                push!(d_i,candidate_den[i]);
                push!(W_i,W[i])
            end
    
        end
        println(W_i)
        if length(findall(x->x==W_i[argmax(W_i)],W_i)) <=length(gens(R))
            d=d_i[argmax(W_i)];
            m=W_i[argmax(W_i)];
            push!(trimmed_w,[w,d,m]);
        end        
end #loop only changes when p and computed den/num changes. Does not depend on candidate den

return(trimmed_w);
end


function numeratorAnsatz(R::Ring,m)
tot_var=gens(R);
resid_var=[];
n=[];
deg_vec=power_monomials(length(tot_var),total_degree(m));
resid_var=tot_var
for i in 0:length(deg_vec) 
    temp=R(1);
    if i==0
        continue
        #push!(n,temp);
    else
        for j in 1:length(deg_vec[i]) 
            temp=temp*resid_var[j]^deg_vec[i][j];
        end
        push!(n,temp);
    end
end

return(n);

end

function interpolateAnsatzTemp(R::Ring,p1::Int,n::Vector,BB::Vector,q::Vector,computed_den::Vector,computed_num::Vector)
 try
    return interpolateAnsatz(R::Ring,p1::Int,n::Vector,BB::Vector,q::Vector,computed_den::Vector,computed_num::Vector)
 catch
    return nothing
 end   
end


function interpolateAnsatz(R::Ring,p1::Int,n::Vector,BB::Vector,q::Vector,computed_den::Vector,computed_num::Vector)
    U=q;
    pr=q[1];
    pr_2=q[2];
    m=q[4];
    d=q[3];
    F=GF(p1);
   
    R1,_= QadicField(ZZ(p1), 1,30);
    Qx, x = QQ["x"];

    B=matrix_space(R1,length(n),1);
    A=matrix_space(R1,length(n),length(n));
    a=Matrix{typeof(R1(1))}(undef,length(n),length(n));
    b=Matrix{typeof(R1(1))}(undef,length(n),1);

    ## We have w, d and n. We generate as many x as length(n) for satisfy (4) for a given fixed prime p
    B_F=matrix_space(QQ,length(n),1);
    A_F=matrix_space(QQ,length(n),length(n));
    a_F=Matrix{typeof(QQ(1))}(undef,length(n),length(n));
    b_F=Matrix{typeof(QQ(1))}(undef,length(n),1);

    B_Q=matrix_space(F,length(n),1);
    A_Q=matrix_space(F,length(n),length(n));
    a_Q=Matrix{typeof(F(1))}(undef,length(n),length(n));
    b_Q=Matrix{typeof(F(1))}(undef,length(n),1);
    for i in 1:length(n)
        p=pr[i+8][1];
        q=pr[i+8][2];
        r=pr[i+8][3];

        p2=pr_2[i+8][1];
        q2=pr_2[i+8][2];
        r2=pr_2[i+8][3];
        constructed_term=0;
            if length(computed_den)!=0
                for k in 1:length(computed_den) 
                    constructed_term=constructed_term+(evaluate(computed_num[k],[p,q,r])/evaluate(computed_den[k],[p,q,r]) )
                end    
            end
    
        b[i,1]=R1(evaluateBB(BB,[p,q,r])-constructed_term)+O(R1,p1^(m-1));
        b_F[i,1]=coeff(lift(Qx,R1(evaluateBB(BB,[p,q,r])-constructed_term+O(R1,p1^(m-1)))),0);
        # b_Q[i,1]=F(numerator(coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r])+O(R1,p1^(m-1)))),0))*
        #print(F(denominator(coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r])+O(R1,p1^(m-1)))),0))))
        de=denominator(b_F[i,1]);
        nu=numerator(b_F[i,1]);
        if F(de)==0
            b_Q[i,1]=F(0);
        else
            de=1/F(de);
            b_Q[i,1]=F(de*nu);
        end
        #b_Q[i,1]=1/F(denominator(coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r])+O(R1,p1^(m-1)))),0)))
        for j in 1:length(n) 
            a[i,j]=R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))+O(R1,p1^(m-1));
            a_F[i,j]=coeff(lift(Qx,R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))+O(R1,p1^(m-1))),0);
            # a[i,j]=R1(evaluate(n[j],[p,q,r]))+O(R1,p1^(m-1));
            # a_F[i,j]=coeff(lift(Qx,R1(evaluate(n[j],[p,q,r]))+O(R1,p1^(m-1))),0);
            de=denominator(a_F[i,j]);
            nu=numerator(a_F[i,j]);
            if F(de)==0
                a_Q[i,j]=F(0);
            else
            de=1/F(de);
            a_Q[i,j]=F(nu*de);
        end
        #    a_Q[i,j]=F(numerator(coeff(lift(Qx,R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))+O(R1,p1^(m-1))),0))*1/F(denominator(coeff(lift(Qx,R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))+O(R1,p1^(m-1))),0))));
    end 
    
end
a=A(a);
b=B(b);
a_F=A_F(a_F);
b_F=B_F(b_F);

a_Q=A_Q(a_Q);
b_Q=B_Q(b_Q);
if det(a_F)==0
    error("Try again with different probes det is zero");
    #q=pickDen(R,R1,common_factors,w,candidate_den,20 );
    #return interpolateAnsatz(R,Int(prime(R1)),n,BB,q,computed_den,computed_num);
   \
else
    inc=inv(a_F)*b_F;
    c=[];
    for i in 1:length(n) 
        #a=(inv(a_F)*b_F)[i];
        #de=denominator(inc[i]);
        
        if F(denominator(inc[i]))==0
            error("Try again with different probes");
        #    println([F(a),0])
        else
        #    nu=1/F(numerator(a));
        #    println([F(a),F(de*nu)])
        push!(c,F(inc[i]))
            #---------------------------------
            p=pr[1][1];
            q=pr[1][2];
            r=pr[1][3];
            constructed_term=0;
            if length(computed_den)!=0
                for k in 1:length(computed_den) 
                    constructed_term=constructed_term+(evaluate(computed_num[k],[p,q,r])/evaluate(computed_den[k],[p,q,r]) )
                end    
            end
    
            lhs=R1(evaluateBB(BB,[p,q,r])-constructed_term)+O(R1,p1^(m-1));
            #lhs=R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1))
            lhs_F=coeff(lift(Qx,lhs),0)
            C=matrix_space(R1,1,length(n))
            C_F=matrix_space(QQ,1,length(n))
            c_l=Matrix{typeof(R1(1))}(undef,1,length(n))
            c_F=Matrix{typeof(QQ(1))}(undef,1,length(n))
            for i in 1:length(n) 
                c_l[1,i]=R1(evaluate(n[i],[p,q,r])/evaluate(d,[p,q,r]))
                c_F[1,i]=coeff(lift(Qx,R1(evaluate(n[i],[p,q,r])/evaluate(d,[p,q,r]))),0)
            end 
            c_l=C(c_l)
            c_F=C_F(c_F)
            println([R1((c_F*inv(a_F)*b_F)[1])+O(R1,p1^(m-1)),R1(lhs_F)+O(R1,p1^(m-1))])
            #---------------------------------
        end
        #push!(c,F(inc[i]))
    end    
    return(c);
end
end 

function getWeights(R::Ring,R1::Field,common_factors1::Vector)
P=Int(prime(R1));
F=GF(P);
Qx, x = QQ["x"];
t=rand(10000:100000);
vars=gens(R);
s=[];
candidate_w=[];
common_factors,i_w=getIndDen(R,common_factors1);
#for i in 1:length(i_w) 
for i in 1:length(common_factors1)
push!(s,0:8);
end
#tot_weights=enumerateTerms([0:4,0:4,0:4]);
tot_weights=enumerateTerms(s);
tot_weights=deleteat!(tot_weights,[1])
print(tot_weights)
#---------------------------------------------------------
for o in 1:length(tot_weights) 
    w1=tot_weights[o];
    w=Vector{typeof(1)}(undef,length(vars));
    for i in 1: length(vars)
        w[i]=w1[i];
    end
    for i in 1:length(i_w)
        w[i]= w1[i_w[i]];
    end
    println(w)
    #undef_indices=findall(x->x===nothing || x ===missing,w);
    #for i in 1:length(undef_indices) 
    #    w[undef_indices[i]]=rand(0:4);
    #end
    #println(w)
    X=matrix_space(QQ,length(common_factors),length(vars));
    x=Matrix(zero_matrix(QQ,length(common_factors),length(vars)))
    
    for i in 1:length(common_factors) 
        f=common_factors[i];
        while f!=0
            for j in 1:length(vars) 
                if leading_monomial(f)==vars[j]
                    x[i,j]=leading_coefficient(f)#R1(leading_coefficient(f));                   
                end
            end   
            f=tail(f); 
        end
    end
    x=X(x);

    C=matrix_space(QQ,length(common_factors),1);
    c=Matrix{typeof(QQ(1))}(undef,length(common_factors),1)
    d=Matrix{typeof(QQ(1))}(undef,length(common_factors),1)
    #P=Int64(101)

    for i in 1:length(common_factors) 
        l1=rand(1:(t*1000));
        l=rand(1:t)//l1;
        x1=R1(l)+O(R1,P^w[i]); 
        p= l-lift(Qx, R1(x1));
        l2=coeff(p,0);
        l3=coeff(lift(Qx, R1(x1)),0);
        while denominator(l2)%P==0 #&& denominator(l3)%P==0 
            l1=rand(1:(t*1000));
            l=rand(1:t)//l1;
            x1=R1(l)+O(R1,P^w[i]); 
            p= l-lift(Qx, R1(x1));
            l2=coeff(p,0);
            l3=coeff(lift(Qx, R1(x1)),0);
        end
    #l=rand(1:t)//l1;
    #x1=R1(l)+O(R1,P^w[i]); 
    #p= l-lift(Qx, R1(x1));
    #c[i,1]=coeff(p,0);
        c[i,1]=l2;
        d[i,1]=coeff(lift(Qx, R1(x1)),0);
    end
    c=C(c)
    d=C(d)
    prob1=Matrix(inv(x)*c)
    prob2=Matrix(inv(x)*d)
    e1=prob1[1,1]
    e2=prob1[2,1]
    e3=prob1[3,1]

    if F(coeff(lift(Qx,R1(e1)),0))!=0 || F(coeff(lift(Qx,R1(e2)),0))!=0||F(coeff(lift(Qx,R1(e3)),0))!=0
        w_j=[]
        for j in 1:length(common_factors1) 
            push!(w_j,valuation(R1(evaluate(common_factors1[j],[e1,e2,e3]))));
        end
        if length(findall(x->x==w_j,candidate_w))==0
        push!(candidate_w,w_j)
        end
    end
end
#=---------------------------------------------------------
p=Int(prime(R1));
Qx, x = QQ["x"];
vars=gens(R);
candidate_w=[];
F=GF(p);
#common_factors=getIndDen(R,common_factors);
for i in 2:length(tot_weights) 
    w=collect(tot_weights[i]);
    prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000))
    p1=prob1[1,1]
    q1=prob1[2,1]
    r1=prob1[3,1]
    w_j=[];
    for j in 1:length(common_factors) 
        push!(w_j,valuation(R1(evaluate(common_factors[j],[p1,q1,r1]))));
    end
    if F(coeff(lift(Qx,R1(p1)),0))!=0 && length(findall(x->x==w_j,candidate_w))==0
        push!(candidate_w,w_j)            
    end    
    #println(valuation(R1(evaluate(common_factors[1],[p1,q1,r1]))))
    
end
=#
return(candidate_w);
end

function getCandidateDenominators(R::Ring,common_den::Vector)
s=[];
common_factors=Vector{typeof(R(1))}(undef,0);
candidate_den=[];
for i in 1:length(common_den) 
    push!(s,0:common_den[i][2])
    push!(common_factors,R(common_den[i][1]))
end
exp_vec=enumerateTerms(s);
exp_vec=deleteat!(exp_vec,[1,length(exp_vec)]);
for i in 1:length(exp_vec) 
   c=1;
   for j in 1:length(common_factors) 
    c=c*common_factors[j]^exp_vec[i][j]
   end
    push!(candidate_den,c);#common_den[1][1]^exp_vec[i][1]*common_den[2][1]^exp_vec[i][2]*common_den[3][1]^exp_vec[i][3]);
end
return([candidate_den,exp_vec,common_factors]);
end

function setBB(R::Ring,common_den::Vector,common_num)
f=1;
for i in 1:length(common_den) 
    f=f*common_den[i][1]^common_den[i][2];
end
BB=[common_num,f];
return(BB);
end

function computeCoef(R::Ring,ansatzC::Vector,ansatzP::Vector,prime_i::Vector)
    try
        f=0;
        for i in 1:length(ansatzP[1]) 
            c_i=Vector{Int64}(undef,0);
            for j in 1:length(prime_i) 
                push!(c_i,lift(ZZ,ansatzC[j][i]))
            end 
            c=manyRR(c_i,prime_i);
            f=f+R(c[1]//c[2])*ansatzP[1][i];
        end
            return f
    catch
        return nothing
    end
end

