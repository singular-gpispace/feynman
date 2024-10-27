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
for i in 2:length(exp_vec) 
    #toDo write a function to obtain the correspoding denominator 
    push!(candidate_den,common_den[1][1]^exp_vec[i][1]*common_den[2][1]^exp_vec[i][2]*common_den[3][1]^exp_vec[i][3]);
    push!(W,w[1]*exp_vec[i][1]+w[2]*exp_vec[i][2]+w[3]*exp_vec[i][3]);
end
    candidate_den
    W
#Task3.2: For given prime p, find a valuatio n point x such that for all i in (1:length(w)), p^{-w_i}=|f_i(x)|_p 
p=101;
tot_weights=enumerateTerms([0:4,0:4,0:4]);
length(tot_weights)
w=collect(tot_weights[1])
prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000))
    p1=prob1[1,1]
    q1=prob1[2,1]
    r1=prob1[3,1]

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

candidate_w[11]

common_factors
prob1,prob2=generating_probes(R,p,common_factors,candidate_w[11],rand(1000:10000));
p1=prob1[1,1]
q1=prob1[2,1]
r1=prob1[3,1]
F(p1)
valuation(R1(evaluate(common_factors[1],[p1,q1,r1])))
valuation(R1(evaluate(common_factors[2],[p1,q1,r1])))
valuation(R1(evaluate(common_factors[3],[p1,q1,r1])))
R1(r1)
F(p1)
q1
r1





#-------------------------------
##Check the chosen point does the job
###Task4: Filter the candidate denominators {d_i} by using probes at the chosen points
##Evaluate R(x) at x computed above

w=trimmed_w[1][1]
prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000))
p1=prob1[1,1]
q1=prob1[2,1]
r1=prob1[3,1]


a=evaluate(common_num,[p1,q1,r1])
f=1;
for i in 1:length(common_den) 
    f=f*common_den[i][1]^common_den[i][2];
end
b=evaluate(f,[p1,q1,r1])
R_w=valuation(R1(a/b))

##Choosing the candidate denominators
d_i=[];
W_i=[];
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
    push!(W_i,W[i-1])
end

end
d_i
W_i
W
##Find suitable w such that there is only one candidate for denominator
trimmed_w=[];
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
    for i in 2:length(exp_vec) 
        #toDo write a function to obtain the correspoding denominator 
        push!(W,w[1]*exp_vec[i][1]+w[2]*exp_vec[i][2]+w[3]*exp_vec[i][3]);
    end
    prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000));
    p1=prob1[1,1]
    q1=prob1[2,1]
    r1=prob1[3,1]
    a=evaluate(common_num,[p1,q1,r1]) 
    b=evaluate(f,[p1,q1,r1])
    R_w=valuation(R1(a/b))
    d_i=[];
    W_i=[];
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
                push!(W_i,W[i-1])
            end
    
        end
        if length(findall(x->x==W_i[argmin(W_i)],W_i))==1
            d=d_i[argmin(W_i)]
            m=W_i[argmin(W_i)]
            push!(trimmed_w,[w,d,m]);
        end        
end

for i in 1:length(trimmed_w) 
    println(trimmed_w[i][3])
end
trimmed_w


candidate_w
##Pick the candidate denominator from the list of denominators

d=d_i[argmin(W_i)]
m=W_i[argmin(W_i)]
typeof(m)
d=trimmed_w[3][2]
m=trimmed_w[3][3]
w=trimmed_w[3][1]
###Task5: Reconstruct the numerator of one candidate by performing additional probes
#construct ansatz for numerator
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
deg_vec=power_monomials(3,1)
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

## We have w, d and n. We generate as many x as length(n) for satisfy (4) for a given fixed prime p

p1=p;
F=GF(p1);


B=matrix_space(R1,length(n),1);
A=matrix_space(R1,length(n),length(n));
a=Matrix{typeof(R1(1))}(undef,length(n),length(n));
b=Matrix{typeof(R1(1))}(undef,length(n),1);
common_num
d
f
m
p1=101
(denominator(12//101))%101
det(a)
#for i in 1:length(trimmed_w) 
    w=trimmed_w[3][1]
 #while true
    for i in 1:length(n)

        prob1,prob2=generating_probes(R,101,common_factors,w,100000)#rand(1000:10000))
        p=prob1[1,1]
        q=prob1[2,1]
        r=prob1[3,1]
        #println(p,",",q,",",r)
        s=1
        while F(evaluate(f,[p,q,r]))==0 && s<1000 #&& R1(evaluate(f,[p,q,r]))==0 && valuation(R1((evaluate(n[1],[p,q,r]))))<=0 && valuation(R1((evaluate(n[2],[p,q,r]))))<=0 && valuation(R1((evaluate(n[3],[p,q,r]))))<=0 #denominator(evaluate(f,[p,q,r]))%p1==0 && F(evaluate(n[1],[p,q,r]))==0 && F(evaluate(n[2],[p,q,r]))==0 F(evaluate(n[3],[p,q,r]))==0#F(evaluate(n[1],[p,q,r])/evaluate(d,[p,q,r]))==0 && F(evaluate(n[2],[p,q,r])/evaluate(d,[p,q,r]))==0 && F(evaluate(n[3],[p,q,r])/evaluate(d,[p,q,r]))==0 #denominator(evaluate(n[1],[p,q,r]))%p1==0 &&denominator(evaluate(n[2],[p,q,r]))%p1==0 &&denominator(evaluate(n[3],[p,q,r]))%p1==0 && R1(evaluate(d,[p,q,r]))==0#||valuation(R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m+1)))>=0#/evaluate(d,[p,q,r]))+O(R1,p1^(-m+1))==0
            s=s+1;
            prob1,prob2=generating_probes(R,101,common_factors,w,100000)#rand(1000:10000))
        p=prob1[1,1]
        q=prob1[2,1]
        r=prob1[3,1]
        println(p,",",q,",",r)
        end

        
        b[i,1]=R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(-m+1))
       #b[i,1]=R1(evaluate(common_num,[p2,q2,r2])/evaluate(f,[p2,q2,r2]))+O(R1,p1^(-m+1))-(R1(evaluate(common_num,[p2,q2,r2])/evaluate(f,[p2,q2,r2]))+O(R1,p1^(-m+2)))
        for j in 1:length(n) 
            println(p,",",q,",",r)
            println(F(evaluate(n[j],[p,q,r])))
          #  if F(evaluate(n[j],[p,q,r]))==0
          #      a[i,j]= F(evaluate(n[j],[p,q,r]))   
          #  else
          #      a[i,j]=F(evaluate(n[j],[p,q,r]))#/evaluate(d,[p,q,r]))    
          #  end
            #a[i,j]=F(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))
          #  a[i,j]=coeff(lift(Qm,R1(evaluate(n[j],[p,q,r]))+O(R1,1)),0)
          
          a[i,j]=R1(evaluate(n[j],[p,q,r]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#+O(R1,1);
          
          #a[i,j]=R1(evaluate(n[j],[p2,q2,r2]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#+O(R1,1);
          
          #println(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]),R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r])))
        end
    
    end
 #   a=A(a);    
  #  if det(a)!=0
   #     println(w)
    #    break
   # end
#end

for i in 1:length(n)
    prob1,prob2=generating_probes(R,101,common_factors,w,rand(1000:10000))
    p=prob1[1,1]
    q=prob1[2,1]
    r=prob1[3,1]
    println(p,",",q,",",r)
    while evaluate(f,[p,q,r])==0 #&& F(evaluate(n[2],[p,q,r]))==0 && R1(evaluate(d,[p,q,r]))==0#||valuation(R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m+1)))>=0#/evaluate(d,[p,q,r]))+O(R1,p1^(-m+1))==0
    prob1,prob2=generating_probes(R,101,common_factors,w)
    p=prob1[1,1]
    q=prob1[2,1]
    r=prob1[3,1]
    #println(p,",",q,",",r)
    end
    p2=prob1[1,1]
    q2=prob1[2,1]
    r2=prob1[3,1]
    
    b[i,1]=R1(evaluate(common_num,[p2,q2,r2])/evaluate(f,[p2,q2,r2]))+O(R1,p1^(-m+1))
   #b[i,1]=R1(evaluate(common_num,[p2,q2,r2])/evaluate(f,[p2,q2,r2]))+O(R1,p1^(-m+1))-(R1(evaluate(common_num,[p2,q2,r2])/evaluate(f,[p2,q2,r2]))+O(R1,p1^(-m+2)))
    for j in 1:length(n) 
        a[i,j]=F(evaluate(n[j],[p2,q2,r2]))
      #  a[i,j]=coeff(lift(Qm,R1(evaluate(n[j],[p,q,r]))+O(R1,1)),0)
      
      #a[i,j]=R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#+O(R1,1);
      
      #a[i,j]=R1(evaluate(n[j],[p2,q2,r2]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#+O(R1,1);
      
      #println(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r]),R1(evaluate(n[j],[p,q,r])/evaluate(d,[p,q,r])))
    end

end

a
b
typeof(A)
typeof(a)
B=matrix_space(F,length(n),1);
b=Matrix{typeof(F(1))}(undef,length(n),1)
b[1,1]=F(61)
b[2,1]=F(87)
b[3,1]=F(52)
a=A(a);
b=B(b);
det(a)

strong_echelon_form(a)
inv(a)*b
is_invertible(a)
#F(coeff(lift(Qx,(inv(a)*b)[1,1]),0))
##test for computed results
R2,_= QadicField(ZZ(7), 1,30);

R2(1//2)+O(R2,7^4)
valuation(R2(-2401//2))
w
    #---------------------------------------------------------------------------------
    prob1,prob2=generating_probes(R,p1,common_factors,w,10000)
    p=prob1[1,1]
    q=prob1[2,1]
    r=prob1[3,1]
    while evaluate(f,[p,q,r])==0 #|| valuation(R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m+1)))>=0
    #s=1
    #    while F(evaluate(f,[p,q,r]))==0 && s<1000 && R1(evaluate(f,[p,q,r]))==0 && denominator(evaluate(f,[p,q,r]))%p1==0 && F(evaluate(n[1],[p,q,r]))<=0 && F(evaluate(n[2],[p,q,r]))<=0 && F(evaluate(n[3],[p,q,r]))<=0#F(evaluate(n[1],[p,q,r])/evaluate(d,[p,q,r]))==0 && F(evaluate(n[2],[p,q,r])/evaluate(d,[p,q,r]))==0 && F(evaluate(n[3],[p,q,r])/evaluate(d,[p,q,r]))==0 #denominator(evaluate(n[1],[p,q,r]))%p1==0 &&denominator(evaluate(n[2],[p,q,r]))%p1==0 &&denominator(evaluate(n[3],[p,q,r]))%p1==0 && R1(evaluate(d,[p,q,r]))==0#||valuation(R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m+1)))>=0#/evaluate(d,[p,q,r]))+O(R1,p1^(-m+1))==0
    #        s=s+1;
        
    prob1,prob2=generating_probes(R,101,common_factors,w,10000)
    p=prob1[1,1]
    q=prob1[2,1]
    r=prob1[3,1]
 end

    #---------------------------------------------------------------------------------
    
    
 (R1(1201))   
    
  #  lhs=F(coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))),0))    
  #  lhs=coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))),0)    
#F(lhs)
lhs=R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(-m+1))
    C=matrix_space(R1,1,length(n))
    c=Matrix{typeof(R1(1))}(undef,1,length(n))
    p2=prob2[1,1]
    q2=prob2[2,1]
    r2=prob2[3,1]
    c[1,1]=R1(evaluate(n[1],[p,q,r]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#F(coeff(lift(Qx,R1(evaluate(n[1],[p,q,r]))),0));
   c[1,2]=R1(evaluate(n[2],[p,q,r]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#F(coeff(lift(Qx,R1(evaluate(n[2],[p,q,r]))),0));
   c[1,3]=R1(evaluate(n[3],[p,q,r]))#/evaluate(d,[p,q,r]))#+O(R1,p1^(m))#F(coeff(lift(Qx,R1(evaluate(n[3],[p,q,r]))),0));
   c=C(c)
   R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))
   valuation((c*inv(a)*b)[1])
   (inv(a)*b)[1,1]
   c[1,3]*(inv(a)*b)[3,1]
   F(coeff(lift(Qx,(c*inv(a)*b)[1]),0))
   F(coeff(lift(Qx,lhs),0))
   coeff(lift(Qx,(c*inv(a)*b)[1]),0)
   coeff(lift(Qx,lhs),0)
   c*inv(a)*b
   coeff(lift(Qx,lhs),0)
   valuation(lhs)
  lhs==(c*inv(a)*b)[1]
   ###Task6: Repeat steps 4 and 5 to reconstruct the other terms 
   R1(evaluate(common_num,[p2,q2,r2])/evaluate(f,[p2,q2,r2]))#+O(R1,p1^(-m+1))
   R1(1//2)==R1(1201)
   1030301//101



   while coeff(lift(Qx,(c*inv(a)*b)[1]),0)!= coeff(lift(Qx,lhs),0)
    prob=generating_probes(R,p1,common_factors,w)
    p=prob[1,1]
    q=prob[2,1]
    r=prob[3,1]
    while evaluate(f,[p,q,r])==0 || valuation(R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m+1)))<0
    prob=generating_probes(R,101,common_factors,w)
    p=prob[1,1]
    q=prob[2,1]
    r=prob[3,1]
 end
    lhs=R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m+1))
    c=Matrix{typeof(R1(1))}(undef,1,length(n))
     c[1,1]=R1(evaluate(n[1],[p,q,r]/evaluate(d,[p,q,r])))#+O(R1,p1^(m+1))#F(coeff(lift(Qx,R1(evaluate(n[1],[p,q,r]))),0)); 
    c[1,2]=R1(evaluate(n[2],[p,q,r])/evaluate(d,[p,q,r]))#+O(R1,p1^(m+1))#F(coeff(lift(Qx,R1(evaluate(n[2],[p,q,r]))),0));
    c=C(c)
   end