using Feynman
using Test
using Oscar

#include tests
include("funct.jl")



#---------------------------------------------------------------------------------------------------------------
#Row reduction for given IBP set (Assigning values in F_p for t_i and D).

G=Feynman.simple_graph([1,2,3,4],[(1,4),(1,2),(2,3),(3,4),1,2,3,4]);
I=computeIBP(G,[1,1,1,0],8,true);
set_IBP=I.setIBP;
RZ=I.baikovover;
paraind=I.paraind;

printIBP(set_IBP,3)

##Get index set

indV=[];
for i in 1:length(set_IBP) 
    for j in 1:length(set_IBP[i]) 
        if isempty(indV)==true || isempty(findall(x->x==set_IBP[i][j][2],indV))==true
            
            push!(indV,set_IBP[i][j][2]);    
        end
    end
end
indV

rnkV=[];
for i in 1:length(indV) 
    Nprop=0;
    Nid=0;
    r=0;
    s=0;
    for j in 1: length(indV[i])
        Nprop=Nprop+Heviside(indV[i][j]-1/2);
        Nid=Nid+Heviside(indV[i][j]-1/2)*2^(j-1);
        r=r+indV[i][j]*Heviside(indV[i][j]-1/2);
        s=s+abs(indV[i][j])*Heviside(-indV[i][j]+1/2);
    end
    push!(rnkV,vcat(indV[i],[Nprop,Nid,r,s]));
end
rnkV

indV=[vec[1:length(indV[1])] for vec in sort(rnkV,by=x->x[length(indV[1])+1:length(rnkV[1])],rev=true)];
indV
F=GF(101)

y=rand(F,paraind)
x=Vector{typeof(RZ(1))}(undef,0)
x=[]
for i in 1:ngens(RZ)
    if i<=paraind
        push!(x,RZ(lift(ZZ,y[i])));   
    else
        push!(x,RZ(0)); 
    end 
    
end
x

f=set_IBP[1][1][1]

x=convert(AbstractArray,x);

F(Int(coeff(evaluate(f,x),1)))

S=matrix_space(F,length(set_IBP),length(indV));
A=Matrix{typeof(F(1))}(undef,length(set_IBP),length(indV))  

for i in 1:length(set_IBP) 
    for j in 1:length(indV)
        indi=0; 
        for k in 1:length(set_IBP[i]) 
            if set_IBP[i][k][2]==indV[j]
                indi=indi+1;
                A[i,j]=F(Int(coeff(evaluate(set_IBP[i][k][1],x),1)));
            end
        end
        if indi==0
            A[i,j]=F(0);
        end
    end
end

#A=sparse_matrix(S(A); keepzrows= true);
#B,C=echelon_with_transform(A);
A=S(A);
r,A=rref(A)
A

findall(x->x==[1,1,1,0],indV)

##order them



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
common_num=R(x*y*z);

###Task2: Enumarate complete set of candidates for the denominators {d_i}
candidate_den=[];
exp_max=[];
s=1;
exp_vec=[];
exp_vec=enumerateTerms([0:1,0:1,0:1])

for i in 2:length(exp_vec) 
    #toDo write a function to obtain the correspoding denominator 
    push!(candidate_den,common_den[1][1]^exp_vec[i][1]*common_den[2][1]^exp_vec[i][2]*common_den[3][1]^exp_vec[i][3]);
end
    candidate_den

###Task3: Choose p-adic evaluation points
##Task3.1:Obtain w
#3.1.1 Get algebraic independent factors of common denominator form:
common_factors=Vector{typeof(R(1))}(undef,0);
for i in 1:length(common_den) 
    push!(common_factors,R(common_den[i][1]))
end
is_algebraically_independent_with_relations(common_factors)

#3.1.2 Construct w

w=[1,1,4];
common_factors

#Task3.2: For given prime p, find a valuatio n point x such that for all i in (1:length(w)), p^{-w_i}=|f_i(x)|_p 
p=101;
#R1, _ = QadicField(ZZ(101),30);

R1,_= QadicField(ZZ(p), 1,30);
Qx, x = QQ["x"];


x1=R1(rand(1:p-1)//rand(1:p-1))+O(R1,p^4);
y1=R1(rand(1:p-1)//rand(1:p-1))+O(R1,p^1);
z1=R1(rand(1:p-1)//rand(1:p-1))+O(R1,p^4);


p1= lift(Qx, R1(x1))
q1= lift(Qx, R1(y1))
r1= lift(Qx, R1(z1))


k=true;
while k==true
    
    x1=R1(rand(1:p-1)//rand(1:p-1))+O(R1,101^4);
    y1=R1(rand(1:p-1)//rand(1:p-1))+O(R1,101^1);
    z1=R1(rand(1:p-1)//rand(1:p-1))+O(R1,101^4);


p1= lift(Qx, R1(x1))
q1= lift(Qx, R1(y1))
r1= lift(Qx, R1(z1))
j=0;
for i in 1:length(common_factors)
    
    if valuation(R1(evaluate(common_factors[i],[p1,q1,r1])))==w[i]  
    j=j+1;
    println(j,valuation(R1(evaluate(common_factors[i],[p1,q1,r1]))))
    end  
end

if j>1
    k=false;
end

end

##Check the chosen point does the job
valuation(R1(evaluate(common_factors[1],[p1,q1,r1])))
valuation(R1(evaluate(common_factors[2],[p1,q1,r1])))
valuation(R1(evaluate(common_factors[3],[p1,q1,r1])))

###Task4: Filter the candidate denominators {d_i} by using probes at the chosen points

##Evaluate R(x) at x computed above
a=evaluate(common_num,[p1,q1,r1])
f=1;
for i in 1:length(common_den) 
    f=f*common_den[i][1]^common_den[i][2];
end
b=evaluate(f,[p1,q1,r1])
R_w=valuation(R1(a/b))

##Choosing the candidate denominators
d_i=[];
for i in 2:length(exp_vec) 
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
end

end

##Pick the candidate denominator from the list of denominators
d=d_i[1];
for i in 2:length(d_i) 
    println(d_i[i],valuation(R1(evaluate(d_i[i],[p1,q1,r1]))))
if valuation(R1(evaluate(d_i[i],[p1,q1,r1])))<valuation(R1(evaluate(d,[p1,q1,r1])))
    d=d_i[i];
end    
end


###Task5: Reconstruct the numerator of one candidate by performing additional probes
#construct ansatz for numerator
d
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

deg_vec=power_monomials(length(resid_var),total_degree(d))

for i in 0:length(deg_vec) 
    temp=R(1);
    if i==0
        push!(n,temp);
    else
        for j in 1:length(deg_vec[i]) 
            temp=temp*resid_var[j]^deg_vec[i][j];
        end
        push!(n,temp);
    end
end

n

## We have w, d and n. We generate as many x as length(n) for satisfy (4) for a given fixed prime p
p1=101;

F=GF(p1);
B=matrix_space(F,length(n),1);
A=matrix_space(F,length(n),length(n));
a=Matrix{typeof(F(1))}(undef,length(n),length(n))
b=Matrix{typeof(F(1))}(undef,length(n),1)

for i in 1:length(n)
    Qm, m = QQ["m"]     

    x1=1//rand(1:p1-1)+O(R1,101^(w[1]));
    y1=1//rand(1:p1-1)+O(R1,101^(w[2]));
    z1=1//rand(1:p1-1)+O(R1,101^(w[3]));
    Qm, m = QQ["m"]
    p= lift(Qm, R1(x1))
    q= lift(Qm, R1(y1))
    r= lift(Qm, R1(z1))
   
    print(F(coeff(lift(Qm,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))),0)))
    #print(coeff(evaluate(f,[p,q,r]),0))
    #b[i,1]=lift(Qm,R1(coeff(evaluate(common_num,[p,q,r]),0)/coeff(evaluate(f,[p,q,r]),0)))
    #b[i,1]=lift(Qm,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r])))
    b[i,1]=F(coeff(lift(Qm,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))),0));
    for j in 1:length(n) 
        a[i,j]=F(coeff(lift(Qm,R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))),0))
    end
end

a
b
typeof(A)
typeof(a)
a=A(a)
b=B(b)
inv(a)*b

##test for computed results

    x1=1//rand(1:p1-1)+O(R1,101^(w[1]));
    y1=1//rand(1:p1-1)+O(R1,101^(w[2]));
    z1=1//rand(1:p1-1)+O(R1,101^(w[3]));
    Qm, m = QQ["m"]
    p= lift(Qx, R1(x1))
    q= lift(Qx, R1(y1))
    r= lift(Qx, R1(z1))
    lhs=F(coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))),0))    
    C=matrix_space(F,1,length(n))
    c=Matrix{typeof(F(1))}(undef,1,length(n))
  
    c[1,1]=F(coeff(lift(Qx,R1(evaluate(n[1],[p,q,r])/evaluate(d,[p,q,r]))),0));
   c[1,2]=F(coeff(lift(Qx,R1(evaluate(n[2],[p,q,r])/evaluate(d,[p,q,r]))),0));
   c[1,3]=F(coeff(lift(Qx,R1(evaluate(n[3],[p,q,r])/evaluate(d,[p,q,r]))),0));
   c=C(c)
   lhs==c*inv(a)*b

###Task6: Repeat steps 4 and 5 to reconstruct the other terms 