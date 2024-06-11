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
export ISP
export removeElimVars
export computeBaikovMatrix
export makePoly
export removeVariableLocal
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
#o=deglex(para)*deglex(var);
#od=lex(para,var);
od=invlex(RP);
#--------------------------------------
  I=standard_basis(I,ordering=od,complete_reduction=true);
  #G1=G;
  eliminatedVariables=Vector{typeof(RP(1))}(undef,0);
  G.elimvar=Vector{typeof(RP(1))}(undef,0);
  n=length(I);
  for i in 1:n
    ld=leading_term(I[i],ordering=od);
    ta=ld-I[i];
    ld=RP(ld);

    if length(G.elimvar)==0
        push!(G.elimvar,ld);
    else
        k=0;
        for j in 1:length(G.elimvar)
            if G.elimvar[j]==ld || G.elimvar[j]==-ld 
                k=k+1;
                #print(eliminatedVariables);
            end 
            
        end
        if k==0
            push!(G.elimvar,ld);            
        end
        
    end
    ta=RP(ta); 
    G=substituteGraph(G,ld,ta);  
  end  
  #G.elimvar=Vector{typeof(RP(1))}(undef,0);
  #G.elimvar=eliminatedVariables;
  return G;
end

#Assume R is a ring.
#Return a polynomial ring with j-th variable removed
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
        print(genR);

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
    
    J=groebner_basis(J,ordering=lex(RP),complete_reduction=true);

    v=Vector{typeof(RP(1))}(undef,0);
    for i in 1:length(J)
    push!(v,T(J[i])); 
    end
    J=ideal(RP,v);

    w=Vector{typeof(RP(1))}(undef,0);   
    gens_RP=gens(RP);
    for i in 1:length(gens_RP) 
        Q,h=reduce_with_quotients(gens_RP[i]^2,v,ordering=lex(RP),complete_reduction=true);
        if h!=0
            push!(w,h);
        end
        for j in i+1:length(gens_RP) 
        Q,h=reduce_with_quotients(gens_RP[i]*gens_RP[j],v,ordering=lex(RP),complete_reduction=true);
        if h!=0
            push!(w,h);
        end
        end
    end
    
    K=ideal(RP,w);
    K=standard_basis(K,ordering=lex(RP),complete_reduction=true);
    H=hom(RP,S,gens(S));
    u=Vector{typeof(S(1))}(undef,0);
    for i in 1:length(K) 
        push!(u,H(K[i]));
    end
    K=ideal(S,u);
    return K;
     
end
#This function removes the variables from G.elimvars. This key is generated by the procedure eliminateVariables
#Assume G is a labeled graph
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
    for i in 1:length(G.labels) 
      push!(u,T(G.labels[i]));
    end
    G.labels=u;
    G.over=R1;
    G.overpoly=RP1;
    G.elimvar=[];
     return G; 
end

#computeBaikovMatrix(G)
##Assume G is a simple_graph or G is a labeledgraph where redundent variables have been eliminated by the 
##the procedure reliminatedVariables, and deleted from the ring by the procedure removeElimVars
##Output will  a labelled graph G1, computer the Baikov matrix of defined in G1.baikovover and stores it in G1.baikovmatrix
#---ToDo: Multipy zvar and A and rest of the calculation. But since zvar and A are matrices over different polynomial ring, want to find a way to multiply
    #constructing zvar and A as julia matrix type, we can't multiply get the errror:Incompatible polynomial rings in polynomial operation
    #constructing zvar and A as sparse matrices, we get the error: Unable to coerce polynomial
function computeBaikovMatrix(G)
    if typeof(G)=="graph"
        lG=labelGraph(G,0);
        G1=eliminateVariables(lG);
        G2=removeElimVars(G1);
        return computeBaikovMatrix(G2);
        
    end
    
   # if typeof(G)!="labeledgraph"
    #    throw(error("expected a graph or labeledgraph"));
    # end

    R=G.over;
    RP=G.overpoly;
    P=Feynman.propagators(G);

    I=Feynman.ISP(G);
    PI=P+I;
    idx=0;

    #calculate number of parameters
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

    gram1=Vector{typeof(varRP[1])}(undef,0);

    H=hom(R,RP,gens(RP));

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

    for i in 1:startvars-1 
        for j in i+1:startvars-1 
            PI=PI+ideal(R,[varRP[i]*varRP[j]]);
        end
    
    end

#---getting minimal generating set

    a=Vector{typeof(RP(1))}(undef,0);
    for i in 1:ngens(PI) 
        push!(a,H(PI[i]));
    end
    PI=ideal(RP,a);
    v=Vector{typeof(RP(1))}(undef,0);
    for i in 1:ngens(PI) 
        push!(v,PI[i]);
    end
    PI=ideal(RP,v);

#---------------------------------

    v=Vector{typeof(RP(1))}(undef,0);
    for i in 1:ngens(PI) 
        push!(v,H(PI[i]));
    end

    m=length(para);
    m2=Int(m*(m-1)/2);  
    mt=m2-1;
    n=ngens(PI)-m2;

    Z=Feynman.makePoly(mt,n);
    t=gens(Z);
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

    for i in 1:n    
        zvar[1,i]=z[i]; 
    end

    for i in 1:m2
        zvar[1,i+n]=pq[i];     
    end
    
   
    X=ideal(RP,gram1);
    #is_subset(X,PI)
    H,I=groebner_basis_with_transformation_matrix(PI,ordering=lex(RP),complete_reduction=true);
    T,J=groebner_basis_with_transformation_matrix(X,ordering=lex(RP),complete_reduction=true);
    I1=Matrix{typeof(Z(1))}(undef,ngens(PI),length(H));
    J1=Matrix{typeof(Z(1))}(undef,ngens(X),length(T));
    
    for i in 1:ngens(PI) 
        for j in 1:length(H)
            I1[i,j]=Z(constant_coefficient(I[i,j]));
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

    #is_subset(ideal(RP,w),ideal(RP,v))

    A=Matrix{typeof(Z(1))}(undef,length(v),length(gram1));
    D=Matrix{typeof(Z(1))}(undef,length(v),length(gram1));

    for i in 1:length(gram1)
         u=reduce_with_quotients(gram1[i],v,ordering=lex(RP),complete_reduction=true);
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

function makePoly(n::Int,m::Int)
    Z,t,z=polynomial_ring(QQ,"t"=>(1:n),"z"=>(1:m));
    return Z;
end