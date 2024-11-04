using Feynman
using Test
using Oscar

#include tests
include("funct.jl")



                ############################################################
                # Novel approach to reconstruction of rational functions   #
                ############################################################

#---------------------------------------------------------------------------------------------------------------
##Process Summary : R is the rational function that we want to reconstruct in partial practioned form
###Task1: Find commmon denominator-form denominators of R
R,(x,y,z)=polynomial_ring(QQ,["x","y","z"]);
common_den=[];
push!(common_den,[R(x+y),1]);
push!(common_den,[R(y+z),1]);
push!(common_den,[R(x+z),1]);
#common_num=R(x*y*z);
common_num=R(x^2*(y+2*z)+x*(2*y^2+3*y*z+z^2)+y^2*z+2*y*z^2)
f=1;
    for i in 1:length(common_den) 
    f=f*common_den[i][1]^common_den[i][2];
end
BB=[common_num,f];
###Task2: Enumarate complete set of candidates for the denominators {d_i}
exp_max=[];
s=1;
exp_vec=[];
exp_vec=enumerateTerms([0:1,0:1,0:1])
common_factors=Vector{typeof(R(1))}(undef,0);
for i in 1:length(common_den) 
    push!(common_factors,R(common_den[i][1]))
end
###Task3: Choose p-adic evaluation points
##Task3.1:Obtain w and getting all possible denominators
#3.1.1 Get algebraic independent factors of common denominator form:
is_algebraically_independent_with_relations(common_factors)
#3.1.2 Construct w
candidate_den=[];
W=[];
w=[0,0,1];
#3.1.2 Getting all possible denominators
exp_vec=deleteat!(exp_vec,[1,length(exp_vec)]);
exp_vec
for i in 1:length(exp_vec) 
    #toDo write a function to obtain the correspoding denominator 
    push!(candidate_den,common_den[1][1]^exp_vec[i][1]*common_den[2][1]^exp_vec[i][2]*common_den[3][1]^exp_vec[i][3]);
    push!(W,w[1]*exp_vec[i][1]+w[2]*exp_vec[i][2]+w[3]*exp_vec[i][3]);
end
    candidate_den
    W

    #----There we have the for loop


    #-----------------------------------
    #Task3.2: For given prime p, find a valuatio n point x such that for all i in (1:length(w)), p^{-w_i}=|f_i(x)|_p 
p=113;
tot_weights=enumerateTerms([0:4,0:4,0:4]);
R1,_= QadicField(ZZ(p), 1,30);
Qx, x = QQ["x"];
vars=gens(R);
candidate_w=[];
F=GF(p)
for i in 2:length(tot_weights) 
    w=collect(tot_weights[i])
    
        prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000))
        p1=prob1[1,1]
        q1=prob1[2,1]
        r1=prob1[3,1]
      if F(coeff(lift(Qx,R1(p1)),0))!=0 
            push!(candidate_w,w)
            
        end    
    println(valuation(R1(evaluate(common_factors[1],[p1,q1,r1]))))
    
end
#-------------------------------
##Check the chosen point does the job
###Task4: Filter the candidate denominators {d_i} by using probes at the chosen points
##Evaluate R(x) at x computed above
computed_num=[];
computed_den=[];
trimmed_w=[];
trimmed_w=pickWeights(R,R1,p,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
argmin(x->x[1],trimmed_w)[1]
d=argmin(x->x[1],trimmed_w)[2]
m=argmin(x->x[1],trimmed_w)[3]
w=argmin(x->x[1],trimmed_w)[1]
##Find suitable w such that there is only one candidate for denominator
#output:trimmed_w,d_i,W_i
#Input:R,R1,p,BB,common_den,candidate_w,exp_vec
#funct:pickWeights

d_i=[];
W_i=[];
f=1;
    for i in 1:length(common_den) 
    f=f*common_den[i][1]^common_den[i][2];
end
f
for i in 1:length(candidate_w) 
    w=candidate_w[i]
    W=[];
    for i in 2:length(exp_vec)-1 
        #toDo write a function to obtain the correspoding denominator 
        e=0;
        for j in 1:length(exp_vec[i]) 
            e=e+w[j]*exp_vec[i][j]
        end
        #push!(W,w[1]*exp_vec[i][1]+w[2]*exp_vec[i][2]+w[3]*exp_vec[i][3]);
        push!(W,e);
    end
    prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000));
    p1=prob1[1,1]
    q1=prob1[2,1]
    r1=prob1[3,1]
    #a=evaluate(common_num,[p1,q1,r1])
    a=evaluateBB(BB,[p1,q1,r1]); 
    b=evaluate(f,[p1,q1,r1])
    
    
   # R_w=valuation(R1(a/b))
    
    while valuation(R1(b))==0
    prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000));
    p1=prob1[1,1]
    q1=prob1[2,1]
    r1=prob1[3,1]
    #a=evaluate(common_num,[p1,q1,r1]) 
    a=evaluateBB(BB,[p1,q1,r1]);
    b=evaluate(f,[p1,q1,r1])
       
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
        for i in 2:length(exp_vec)-1 
            lhs=0;
            v1=[];
            v2=[];
            for j in 1:length(w) 
                lhs=lhs-exp_vec[i][j]*w[j];
                if exp_vec[i][j]!=0
                    push!(v1,common_den[j][1]);
                end
                if exp_vec[i][j]!=0 || w[j]!=0
                    push!(v2,common_den[j][1]);
                end
            end    
            I1=ideal(R,v1);
            I2=ideal(R,v2);
            println(lhs);
            if R_w>lhs && I1==I2
                continue;
            else
                push!(d_i,candidate_den[i-1]);
                push!(W_i,W[i-1])
            end
    
        end
        if length(findall(x->x==W_i[argmin(W_i)],W_i))==1
            d=d_i[argmin(W_i)]
            m=W_i[argmin(W_i)]
            push!(trimmed_w,[w,d,m]);
        end        
end #loop only changes when p and computed den/num changes. Does not depend on candidate den


d=d_i[argmin(W_i)]
m=W_i[argmin(W_i)]
w=trimmed_w[argmin(W_i)][1]
##Pick the candidate denominator from the list of denominators
##---------------------------------
##Also Find proeb points such that only have one denominator. 
#output: d,m,pr and pr_2
#input : ring R,p,common_factors,w candidate_den,i 
#funct :pickDen
q=pickDen(R,R1,p,common_factors,w,candidate_den,20 );
pr=q[1];
pr_2=q[2];
m=q[4]
d=q[3]
##-------------------------------------
###Task5: Reconstruct the numerator of one candidate by performing additional probes
#construct ansatz for numerator
#intput:R,m,
#output:n
n=numeratorAnsatz(R,d);
n
tot_var=gens(R);
resid_var=[];
## getting residual variables associated to factor d
for i in 1:length(tot_var) 
     if tot_var[i]!=leading_monomial(d)
        push!(resid_var,tot_var[i]);
    end
end

## Use power_monomials function to costruct ansatz for numerator

n=[];

#deg_vec=power_monomials(length(resid_var),total_degree(d))
deg_vec=power_monomials(3,2)
push!(resid_var,tot_var[3])
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

n

##output: interpolateAnsatz(p1,n,)
p1=113;
interpolateAnsatz(R,p1,n,BB,q)
##We have w, d and n. We generate as many x as length(n) for satisfy (4) for a given fixed prime p
B=matrix_space(R1,length(n),1)
A=matrix_space(R1,length(n),length(n))
a=Matrix{typeof(R1(1))}(undef,length(n),length(n))
b=Matrix{typeof(R1(1))}(undef,length(n),1)
common_num
## We have w, d and n. We generate as many x as length(n) for satisfy (4) for a given fixed prime p
p1=113;
F=GF(p1);
B_F=matrix_space(QQ,length(n),1)
A_F=matrix_space(QQ,length(n),length(n))
a_F=Matrix{typeof(QQ(1))}(undef,length(n),length(n))
b_F=Matrix{typeof(QQ(1))}(undef,length(n),1)
m
B_Q=matrix_space(F,length(n),1)
A_Q=matrix_space(F,length(n),length(n))
a_Q=Matrix{typeof(F(1))}(undef,length(n),length(n))
b_Q=Matrix{typeof(F(1))}(undef,length(n),1)
p=pr[9][1]
    q=pr[9][2]
    r=pr[9][3]
    typeof(r)
    (denominator(r)*r)
    evaluate(f,[p,q,r])
    f
    R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1))
R1(evaluate(n[1],[p,q,r])/evaluate(d,[p,q,r]))+O(R1,p1^(m-1))
lift(Qx,R1(evaluate(n[1],[p,q,r]))+O(R1,p1^(m-1)))
p1^(m-1)
p1
R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1))
for i in 1:length(n)
    p=pr[i+8][1];
    q=pr[i+8][2];
    r=pr[i+8][3];

    p2=pr_2[i+8][1];
    q2=pr_2[i+8][2];
    r2=pr_2[i+8][3];
    
    b[i,1]=R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1));
    b_F[i,1]=coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r])+O(R1,p1^(m-1)))),0)
   # b_Q[i,1]=F(numerator(coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r])+O(R1,p1^(m-1)))),0))*
   #print(F(denominator(coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r])+O(R1,p1^(m-1)))),0))))
   de=denominator(b_F[i,1])
   nu=numerator(b_F[i,1])
   if F(de)==0
    b_Q[i,1]=F(0);
   else
    de=1/F(de)
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
        a_Q[i,j]=F(0)
       else
        de=1/F(de)
        a_Q[i,j]=F(nu*de)
       end
    #    a_Q[i,j]=F(numerator(coeff(lift(Qx,R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))+O(R1,p1^(m-1))),0))*1/F(denominator(coeff(lift(Qx,R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))+O(R1,p1^(m-1))),0))));
    end 
    
end
a=A(a)
b=B(b)

a_F=A_F(a_F)
b_F=B_F(b_F)

a_Q=A_Q(a_Q)
b_Q=B_Q(b_Q)

det(a_F)
det(a)
det(a_Q)
m
w
F(coeff(lift(Qx,R1(49*101^0 + O(R1,101^1))),0))
49*101^0 + O(R1,101^1)
F(49*101^0 + O(101^1))
inv(a)*b
inv(a_F)*b_F
for i in 1:6 
a=(inv(a_F)*b_F)[i];
de=denominator(a);

if F(numerator(a))==0
    println([F(a),0])
else
    nu=1/F(numerator(a));
    println([F(a),F(de*nu)])
end
    
end
#------------------------------------------------------------------------
#------------------------------------------------------------------------
p=pr[1][1];
q=pr[1][2];
r=pr[1][3];
lhs=R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1))
lhs_F=coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1))),0)
    C=matrix_space(R1,1,length(n))
    C_F=matrix_space(QQ,1,length(n))
    c=Matrix{typeof(R1(1))}(undef,1,length(n))
    c_F=Matrix{typeof(QQ(1))}(undef,1,length(n))
    c[1,1]=R1(evaluate(n[1],[p,q,r]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#F(coeff(lift(Qx,R1(evaluate(n[1],[p,q,r]))),0));
   c[1,2]=R1(evaluate(n[2],[p,q,r]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#F(coeff(lift(Qx,R1(evaluate(n[2],[p,q,r]))),0));
   c[1,3]=R1(evaluate(n[3],[p,q,r]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m))#F(coeff(lift(Qx,R1(evaluate(n[3],[p,q,r]))),0));
     for i in 1:length(n) 
        c[1,i]=R1(evaluate(n[i],[p,q,r])/evaluate(d,[p,q,r]))
        c_F[1,i]=coeff(lift(Qx,R1(evaluate(n[i],[p,q,r])/evaluate(d,[p,q,r]))),0)
     end 
   c=C(c)
   c_F=C_F(c_F)
   R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))
   c*inv(a)*b
   c_F*inv(a_F)*b_F
   (c*inv(a)*b)[1]+O(R1,p1^(m-1))

   R1((c_F*inv(a_F)*b_F)[1])+O(R1,p1^(m-1))
   R1(lhs_F)+O(R1,p1^(m-1))

   coeff(lift(Qx,R1(lhs_F)+O(R1,p1^(m-1))),0)
   coeff(lift(Qx,R1((c_F*inv(a_F)*b_F)[1])+O(R1,p1^(m-1))),0)


#---------------rational number reconstruction --------------------------
#------------------------------------------------------------------------
r=44;
m=101;
u=[1,0,m];
v=[0,1,r];
while v[3]>=sqrt(m/2)
    q=floor(u[3]//v[3]);
    r=u-q*v;
    u=v;
    v=r;
end

if abs(v[2])>sqrt(m/2)
    throw(error("given m and r doesnt satisfy the condition in Wang algorithm"));
else
    print([v[3],v[2]]);
end
oneRR(44,101)
#--------------compute r for a given m=p1...pn where a/b=r mod m---------
#------------------------------------------------------------------------
#compute a/b=r mod pi
a=10;
b=11;
mi=[101,103,107];
ri=Vector{Int64}(undef,0)
m=1;
for i in 1:length(mi) 
    push!(ri,ZZ(mod(a*(gcdinv(b,mi[i]))[2],mi[i])));
    m=m*mi[i];
end
typeof(mi)
c=crt(ri,mi)
oneRR(c,m)
typeof([56,29,69])
manyRR(ri,mi)
manyRR([55,72],[107,113])
#------------------------------------------------------------------------
#------------------------------------------------------------------------



