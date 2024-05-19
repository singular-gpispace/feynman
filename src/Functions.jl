using Oscar
include("DataTypes.jl")

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
    printNet(G.edges);
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

#This is the print function used to print a labeled graph
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

#Make a labeled graph
##ASSUME: G is a graph and ch either zero or prime
##ASSUME: Base field is QQ
##labeled graph with polynomialvariables qi at the bounded edges and functin filed variables pi at the unbounded edges over a prime filed of characteristic ch

function labelGraph(G::simple_graph,ch::Int)
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
    else
        F=GF(ch);
        P,p,q=polynomial_ring(F,"p"=>(1:ct),"q"=>(1:anzq));
    end
    
    
    ##R=P_<p(1),...,p(ct)>
    Q=complement_of_prime_ideal(ideal(P,p));
    R,iso=localization(P,Q);

    ##Making list of labels
    pidx=1;
    qidx=1;
    lab=[];
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

#Make balancing ideal;Ideal of balancing condition of the graph
##Assume G is a labeled graph
##This returns ideal of balancing condition of the graph
function balancingIdeal(G::labeledgraph)
 v=G.vertices;
 e=G.edges;
 lab=G.labels;
 R=G.over;
 edg=R(0);
 rel=R(0);
 for i in 1:length(v)
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
    
 end   
 for j in 1:length(e)
    if length(e[j])==1
        rel=rel+R(lab[j]);
        
    end 
    
 end
 I=ideal(R,[edg,rel]);
 return I;
end

#substitute the label 'a' in the labeling by 'b'
#Assume G is a labeled graph
#Return a labelled graph with labelling where each 'a' is replaced 'b'
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

#Elliminate variables on bounded edges
##Assume G is a labeled graph
##This function returns a labeled graph with variables of bounded edges elliminated according to balancing condition
function eliminateVariables(G::labeledgraph)
  RP=G.overpoly;
  R=G.over;
  I=balancingIdeal(G); #Ideal in G.over
  n=ngens(I);
  #Need ring homomorphism from R to RP to get ideal 
  H=hom(R,RP,gens(RP));
  v=[];
  for i in 1:n
    push!(v,H(I[i])); 
  end
  I=ideal(RP,v);
  I=standard_basis(I,ordering=invlex(RP),complete_reduction=true);
  G1=G;
  eliminatedVariables=[];
  n=length(I);
  for i in 1:n
    ld=leading_term(I[i]);
    ta=ld-I[i];
    ld=RP(ld);
    push!(eliminatedVariables,ld);
    ta=RP(ta); 
    G1=substituteGraph(G,ld,ta);  
  end  
  G1.elimvar=eliminatedVariables;
  return G1;
end

#Assume R is a ring.
#Return a polynomial ring with j-th variable removed
function removeVariable(R::Ring,j::Int)
    nv=0;
    if j<0|| j>ngens(R)
        throw(error("Index out of range"));
    end
    genR=deleteat!(gens(R),j);
    v=String[];
    for i in 1:length(genR)
        push!(v,string(genR[i]))       
    end
    R,c=polynomial_ring(coefficient_ring(R),v);
    return R
end

function removeParameter(P::Ring,j::Int)
    R=base_ring(P);
    genR=gens(R);
    nv=0;
    if j<0|| j>ngens(R)
        throw(error("Index out of range"));
    end
    genR=deleteat!(gens(R),j);
    v=String[];
    for i in 1:length(genR)
        push!(v,string(genR[i]))       
    end
    #delete jth variable from base ring
    S,c=polynomial_ring(coefficient_ring(R),v);
    #-------
    u=[genR[1]];
    k=1;
    for i in 2:length(genR) 
        J=ideal(R,u);
        L=complement_of_prime_ideal(J);
        if L==inverted_set(P)
            break;
        end
        push!(u,genR[i]); 
        k=k+1;  
    end
    #-------
    x=[];
    for i in 1:k-1 
        push!(x,c[i]);
    end
    #-------
    w=deleteat!(x,j);
    J=ideal(S,w);
    U=complement_of_prime_ideal(J);
    T,iso=localization(S,U);
    return T;
end

#Compute ideal containing the propagators in the Feynman integral associated to labeled graph G
##Assume G is a labeled graph
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

#Compute the ideal containing the denominators of the Feynman integral
#Assume G is a labeled Graph
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
  for i in 1:ngens(J)
    push!(v,H(J[i])); 
  end
  J=ideal(RP,v);

#Ideal J as an ideal of RP
  v=[];
  for i in 1:ngens(infedges)
    push!(v,H(infedges[i])); 
  end
  infedges=ideal(RP,v);

#reduce J w.r.t std basis of infedges
SB_infedges=standard_basis(infedges,ordering=lex(RP));

N=Vector{typeof(RP(1))}(undef,0);
for i in 1:ngens(J)
    push!(N,RP(J[i]));     
end

M=Vector{typeof(RP(1))}(undef,0);
for i in 1:length(SB_infedges)
    push!(M,RP(SB_infedges[i]));     
end

N=reduce(N,M,ordering=lex(RP),complete_reduction=true);

#rewite the ideal J as an ideal of S
T=hom(RP,S,gens(S));

v=[];
for i in 1:length(N)
push!(v,T(N[i])); 
end

J=ideal(S,v);
return J;
end

#This compute ideal containing the irrreducible scalar products, that is those scalar product which are not linearly dependent on the propagators.
#Assume G is a labelled graph
function ISP(G::labeledgraph)
    
end