using Feynman
using Test
using Oscar


include("funct.jl")


#------------------------------------------
#------------------------------------------
G=Feynman.simple_graph([1,2,3,4],[(1,4),(1,2),(2,3),(3,4),1,2,3,4]);
G=simple_graph([1,2,3,4,5,6],[(6,1),(4,6),(1,2),(3,5),(4,3),(2,5),(5,6),1,2,3,4]);
G=simple_graph([1,2,3,4,5,6,7],[(6,1),(6,4),(1,2),(3,7),(4,3),(2,7),(5,6),(7,5),1,2,3,4,5]);

G=Feynman.labelGraph(G,0);
G=Feynman.eliminateVariables(G);
G=Feynman.removeElimVars(G);
G=Feynman.computeBaikovMatrix(G);
computeIBP(G,3,4);





#--------------------------------------------
RZ=G.baikovover;
    R=G.over;
    gens_RZ=gens(RZ);
    B=G.baikovmatrix;
    k=2;
    D=4;
    cutDeg=4;
    #-----------------------------------------------------------------------------------#
    size(G.baikovmatrix)
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
        #Getting t vector and z vector
        var_t=Vector{typeof(RZ(1))}(undef,0);
        var_z=Vector{typeof(RZ(1))}(undef,0);
        
        for i in 1:mt
            push!(var_t,gens_RZ[i]);
        end
        for i in mt+1:length(gens_RZ) 
            push!(var_z,gens_RZ[i]);
        end
    
        #-----------
        vecNu=[];
        if k>length(var_z)
            error("k must be a nonnegative integer less than to number of Baikov variables.");
        else
            for i in 1:k 
                push!(vecNu,RZ(1));
            end
            for i in k+1:length(var_z) 
                push!(vecNu,RZ(0));
            end
    
        end

        #-----------
        #Computing generators of M1
        t=[];
        n=L*E+Int(L*(L+1)/2);
   
    #------------------------------------------------
    #------------------------------------------------
    ##Compute the partial derivative correctly z_alpha/x(i,j)
    W=Matrix{typeof(RZ(1))}(undef,length(var_z),length(var_z));    
    r=0;
    ij=[];
    for i in 1:E+L-1 
        for j in E+1:E+L 
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
    #------------------------------------------------
    #------------------------------------------------
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
                
                if findall(x->x==[k,i],ij)!=[]
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

    
    #-----------------------------------------------------------------------------

    m=length(var_z);
    
    A=Matrix{typeof(RZ(1))}(undef,2*(length(t[1])-1),2*(length(t[1])-1)+length(t));
    
    #filling vertical block 1
    for j in 1:(length(t[1])-1)
        
        for i in 1:(length(t[1])-1) 
            if i==j
                A[i,j]=RZ(1);
                A[i+m,j]=RZ(1);
                
            else
                A[i,j]=RZ(0);
                A[i+m,j]=RZ(0);
            end
        end 
    end

    for i in length(t[1]):2*(length(t[1])-1)
        for j in 1:(length(t[1])-1) 
            if i-length(t[1])+1==j
                A[i,j]=RZ(1);
            else
                A[i,j]=RZ(0);
            end
        end 
        
    end
    
    #filling vertical block 2
    for i in 1:(length(t[1])-1)
        for j in 1:length(t) 
            A[i,j+length(t[1])-1]=t[j][i];
        end
        
    end

    for i in length(t[1]):2*(length(t[1])-1)
        for j in length(t)+1:length(t)+length(t[1])-1 
            A[i,j]=RZ(0);
        end 
        
    end
    
    #filling vertical block 3
    for i in 1:length(t[1])-1
        for j in 1:(length(t[1])-1) 
            A[i,j+length(t)+length(t[1])-1]=RZ(0);
        end
        
    end
    
    for i in 1:(length(t[1])-1)
        for j in 1:(length(t[1])-1) 
            if i==j
                if j<=k#        j<=length(var_z)
                    A[i+length(t[1])-1,j+length(t)+length(t[1])-1]=RZ(var_z[j]);       
                else
                    A[i+length(t[1])-1,j+length(t)+length(t[1])-1]=RZ(1);#last column is problamatic, as g-i has 9 components and f_1 has 10 components
                end
                
            else
                A[i+length(t[1])-1,j+length(t)+length(t[1])-1]=RZ(0);
            end
        end 
    end
    A
    #-----------------------------------------------------------------------------------#
    # --------------compute module intersection and to get Groebner basis---------------#
    #-----------------------------------------------------------------------------------#
    
    #Defining the polynomial ring RZ in Singular
    v1=Vector{String}(undef,0);
    v2=Vector{String}(undef,0);
    
    for i in 1:length(var_t) 
      push!(v1,"t[$i]");  
    end
    
    for i in 1:length(var_z) 
        push!(v2,String("z[$i]"));  
      end
    
    v=convert(AbstractArray,vcat(v1,v2));
    T,var=Singular.polynomial_ring(Singular.QQ,v,ordering=Singular.ordering_ls(length(v1))*Singular.ordering_dp(),degree_bound=4);
    
    #Build the module using columns of matrix A
    gens_A=[];
    for i in 1: size(A)[2]
        v=Vector{typeof(T(1))}(undef,0);
        u=transpose(A[:,i]);
        for j in 1:length(u) 
          push!(v,T(u[j]))  
        end
        push!(gens_A,Singular.vector(T,v...))
        
    end
    M=Singular.Module(T,gens_A...);
    
    #Compute the syzygies of the module
    G=Singular.syz(M);
    
    #Compute Proj_R^m(sys(M))
    mProj=[];
    for i in 1:Singular.number_of_generators(G) 
        push!(mProj,Singular.vector(T,first(Singular.Base.Array(G[i]),length(t[1])-1)...));
    end
    
    #Compute the Groebner basis G of Proj_R^m(sys(M))
    G=Singular.std(Singular.Module(T,mProj...),complete_reduction=true);
    #G=Singular.jet(G,cutDeg);
    
    #------*****Singular computation is over*****-----------
    
    #Convert generators of Groebner basis (in Singular) to generators in Oscar.
    vecG=[];
    for i in 1:Singular.number_of_generators(G) 
        u=Singular.Base.Array(G[i]);
        v=[];
        for j in 1:length(u) 
            push!(v,RZ(u[j]));
        end
        #Computation of b for given pair (a1,...,am)
        g=RZ(0);
        for l in 1:length(u) 
            g=g+v[l]*derivative(f,length(var_t)+l)
        end
        o=lex(var_z)*neglex(var_t);
        b,h,w=reduce_with_quotients_and_unit(g,[RZ(-f)];ordering=negdegrevlex(RZ),complete_reduction=true);
        if w==0
            push!(v,h[1]);    
            push!(vecG,v);
        end
    end

    #----------------------------------------------------
    prod_z=RZ(1);        
        for j in k+1:m
            prod_z=prod_z*var_z[j];
        end
        u=VecG[1]
        typeof(prod_z*u[3])
        u=vecG[1]
        bj, hj = reduce_with_quotients(prod_z*u[k+1], [var_z[1]], ordering = invlex(RZ));
    #----------------------------------------------------

    #Computation of IBP identities
    
    set_IBP=[];
    for i in 1:length(vecG) 
        u=vecG[i];
        v=[];
        #Computation of bj for j<=k using aj=bjz[j]
        proz=RZ(1);        
        for j in k+1:m
            proz=proz*var_z[j];
        end

        for j in 1:length(u) 
            if j<=k
                if u[j]==0
                    bj=0;
                    hj=0;
                else
                    bj, hj = reduce_with_quotients(u[j], [var_z[j]], ordering = invlex(RZ));
    
                end
                push!(v,bj[1]);
            else
                if u[j]==0
                    bj=0;
                    hj=0;
                else
                   e=j-k;
                    g=proz*u[j];
                    bj, hj = reduce_with_quotients(g, [var_z[e]], ordering = invlex(RZ));
                end
                push!(v,bj[1]);
                
            end
        end

        #Compute Baikov identity associated to generator vecG[i]
        ## (\sum_{i=1}^{m}\frac{\partial a_i}{\partial z_i} - \sum_{i=1}^{k}b_i -\sum_{i=k+1}^{m} \frac{\nu_i a_i}{z_i} - \frac{D-L-E-1}{2}b)
        ##---1---compute m1= (\sum_{i=1}^{m}\frac{\partial a_i}{\partial z_i}
        m1=RZ(0);
        for j in 1:m 
          m1=m1+derivative(u[j],length(var_t)+j);
        end
        ##---2---compute m2=\sum_{i=1}^{k}b_i
        m2=RZ(0);
        for j in 1:k 
            m2=m2+v[j];
        end
        ##---3---compute m3=\sum_{i=k+1}^{m} \frac{\nu_i a_i}{z_i}
        m3=RZ(0);
        for j in (k+1):m 
            m3=m3+vecNu[j]*v[j];
        end
        
        M=(m1-m2-RZ((D-L-E-1)//2)*u[m+1])*proz-m3;
        println(M)
        #Since we multiplied the expression by prod_z, we have to add 1 to each ν_{k+1},...,ν_m
        
        single_IBP=[];
        for j in 1:length(M) 
            degT=degrees(monomial(M,j));
            deg_t=first(degT,length(var_t));
            deg_z=last(degT,length(var_z));
            cTerm=RZ(1);
            for l in 1:length(var_t)
               cTerm=cTerm*var_t[l]^deg_t[l];
            end
            cTerm=cTerm*coeff(M,j);
            iTerm=deg_z-vecNu;
            #Since we multiplied the expression by prod_z, we have to add 1 to each ν_{k+1},...,ν_m
            for l in k+1:length(iTerm) 
                iTerm[l]=iTerm[l]-1;
            end
            push!(single_IBP,[cTerm,iTerm]);
    
        end
        if single_IBP!=[]
            push!(set_IBP,single_IBP);   
        end
    
    
    end
    #if(print)
    set_IBP
    println("IBP identities associated to G setting ν_1=...=ν_k=1 and ν_{k+1}=...=ν_m=0 (Total number of generators without trimming (up to degCut=",cutDeg,")=",length(set_IBP),":")
    for i in 1:Int(length(set_IBP))
        printIBP(set_IBP[i]);    
    end
