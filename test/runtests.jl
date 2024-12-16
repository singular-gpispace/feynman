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
push!(common_den,[R(x-z),1]);
push!(common_den,[R(x+2*y+z),1]);
#----------------------------------------------------------
push!(common_den,[R(x+y),1]);
push!(common_den,[R(y+z),1]);
push!(common_den,[R(x+z),2]);
#------------------------------------------------------------
n=nvars(R);
vars=gens(R);
  common_factors=Vector{typeof(R(1))}(undef,0);
  for i in 1:length(common_den) 
      push!(common_factors,R(common_den[i][1]));
  end
p=101
  R1,_= QadicField(ZZ(p), 1,30);
  candidate_w=getWeights(R,R1,common_factors);

for i in 1:5
    l=getWeights(R,R1,common_factors);
    if length(l)>length(candidate_w)
        candidate_w=l;
    end
    
end

candidate_w
prob1,prob2=generating_probes(R,101,common_factors,[4,3,3,3],rand(10000:1000000));
p1=prob1[1,1]
    q1=prob1[2,1]
    r1=prob1[3,1]
    for i in 1:4 
       print( valuation(R1(evaluate(common_factors[i],[p1,q1,r1]))))
    
    end
    

#------------------------------------------------------------
#------------------------------------------------------------


#common_num=R(x*y*z);
ind_den,P=getIndDen(R,common_factors)
typeof(P[2])

#input:R,common_den
    n=nvars(R)
    vars=gens(R)
    common_factors=Vector{typeof(R(1))}(undef,0);
    for i in 1:length(common_den) 
        push!(common_factors,R(common_den[i][1]))
    end
    common_factors

    prob1,prob2=generating_probes(R,103,common_factors,[8,3,3],rand(10000:1000000));
    prob2
    p1=prob1[1,1]
    q1=prob1[2,1]
    r1=prob1[3,1]
    valuation(R1(evaluate(common_den[3][1],[p1,q1,r1])))


    #----------------------------------------------------------
common_num=R(x^2*(y+2*z)+x*(2*y^2+3*y*z+z^2)+y^2*z+2*y*z^2)
common_num=R(2*x^2*y + x^2*z + x*y^2 - x*y*z + y*z^2)
common_num=R((x^2+y^2)*(x+y)*(x+z)+z*(x+z)^2*(y+z)+y*(x+y)*(x+z)^2)
common_num=R((x+y)*(x+y)*(x+z)+z*(x+z)^2*(y+z)+y*(x+y)*(x+z)^2)
f=1;
s=[];
common_factors=Vector{typeof(R(1))}(undef,0);
    for i in 1:length(common_den) 
    push!(common_factors,R(common_den[i][1]))
end
BB=[common_num,f];
##Input:R,common_num,common_den
##output::BB
##funct::setBB
BB=setBB(R,common_den,common_num)
###Task2: Enumarate complete set of candidates for the denominators {d_i}
p=103
computed_den=[]
computed_num=[]
R1,_= QadicField(ZZ(p), 1,30);
op=getCandidateDenominators(R,common_den);
candidate_den=op[1];
exp_vec=op[2];
exp_vec
common_factors=op[3];
candidate_den
candidate_w=getWeights(R,R1,common_factors);
length(candidate_w)
candidate_w
candidate_w=deleteat!(candidate_w,[1])
trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
trimmed_w


#------------------------------------------------------------------------
#-----------------------------------------------------------------------\
common_factors1=common_factors
P=Int(prime(R1));
F=GF(P);
Qx, x = QQ["x"];
t=rand(10000:100000);
vars=gens(R);
s=[];
candidate_w=[];
common_factors,i_w=getIndDen(R,common_factors1);
for i in 1:length(i_w) 
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
        w[i]=rand(0:8);
    end
    for i in 1:length(i_w)
        w[i]= w1[i];
    end
    #undef_indices=findall(x->x===nothing || x ===missing,w);
    #for i in 1:length(undef_indices) 
    #    w[undef_indices[i]]=rand(0:4);
    #end
    println(w)
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
    c=C(c);
    d=C(d);
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
        println([w_j,candidate_w,"aaa"])
        #if length(findall(x->x==w_j,candidate_w))==0
        push!(candidate_w,w_j)
        #end
    end

end
candidate_w
#------------------------------------------------------------------------\
#-------------------------------------------------------------------------\
op=getCandidateDenominators(R,common_den);
candidate_den=op[1]
exp_vec=op[2];
exp_vec
common_factors=op[3];

w=argmin(x->x[1],trimmed_w)[1]
i=0;
d=0;
m=0;
pr=[];
pr_2=[];
#for k in 1:length(candidate_w)
trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
candidate_den=deleteat!(candidate_den,[1])

trimmed_w=deleteat!(trimmed_w,findall(x->x[2]==candidate_den[1],trimmed_w))    
w=candidate_w[k]; 
p=Int(prime(R1));
nProbes=20
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
#println(pw)
#d=candidate_den[findall(x->x==pw[argmax(pw)],pw)[1]];
   if length(findall(x->x==pw[argmax(pw)],pw))<=length(gens(R)) #&& valuation(R1(evaluate(d,[p1,q1,r1])))!=0#&& evaluate(f,[p2,q2,r2])!=0
       
       push!(pr,[p1,q1,r1]);
       push!(pr_2,[p2,q2,r2]);
       d=candidate_den[findall(x->x==pw[argmax(pw)],pw)[1]];
       print(d)
       m=valuation(R1(evaluate(d,[p1,q1,r1])))
           i=i+1;
           println([m,d])
   end
   
end

#end
m
#--------------------------------------------------------------------------\
#---------------------------------------------------------------------------\
q=[pr,pr_2,d,m]
prime_i=[103]
ansatzP=[];
 ansatzC=[];
 denD=[];
 argmax(x->x[1],trimmed_w)
 trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[3],trimmed_w),trimmed_w))    
 w=argmin(x->x[3],trimmed_w)[1]
computed_den
candidate_den
candidate_den=deleteat!(candidate_den,[3])
 q=pickDen(R,R1,common_factors,w,candidate_den,20 );
 pr=q[1];
 pr_2=q[2];
 m=q[4]
 d=q[3]
 n=numeratorAnsatz(R,d)
 push!(ansatzP,n);
 push!(denD,d);
 modC=interpolateAnsatz(R,103,n,BB,q,computed_den,computed_num);
 modC
 push!(ansatzC,modC);
 f=0;
 for i in 1:length(ansatzP[1]) 
    c_i=Vector{Int64}(undef,0);
    for j in 1:length(prime_i) 
     push!(c_i,lift(ZZ,ansatzC[j][i]))
    end 
    c=manyRR(c_i,prime_i);
     f=f+R(c[1]//c[2])*ansatzP[1][i];
 end
 f
 push!(computed_num,f);
 push!(computed_den,denD[1]);


 #------------------------------------------------------------------------

#-------------------------------------------------------------------

#-------------------------------------------------------------------




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
    #push!(W,w[1]*exp_vec[i][1]+w[2]*exp_vec[i][2]+w[3]*exp_vec[i][3]);
end
    candidate_den
    W

    




#----------------------------------
#output :candidate_den,exp_vec,BB,common_factors
#input  :R,common_den,common_num
#funct  :getCandidateDen(R,common_den) 
op=getCandidateDenominators(R,common_den)
candidate_den=op[1]
exp_vec=op[2]
common_factors=op[3]
candidate_w
trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
   
#----------------------------------

computed_den
exp_vec
candidate_den
exp_vec=deleteat!(exp_vec,[3])
        candidate_den=deleteat!(candidate_den,[3])
#----------------------------------
trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[1],trimmed_w),trimmed_w))    
#----------------------------------
used_w    
w=argmin(x->x[1],trimmed_w)[1]
    computed_num=[];
    computed_den=[];
    candidate_w=getWeights(R,R1,common_factors)
candidate_w=deleteat!(candidate_w,[1])
    used_w=[]
    for i in 1:length(computed_den) 
        exp_vec=deleteat!(exp_vec,findall(x->x==computed_den[i],candidate_den));
        candidate_den=deleteat!(candidate_den,findall(x->x==computed_den[i],candidate_den));
    end
    for i in 1:length(used_w) 
        candidate_w=deleteat!(candidate_w,findall(x->x==used_w[i],candidate_w))
    end
    trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
    candidate_den
    computed_den
    candidate_w
    used_w
    #update candidate denominators
    fil_can=[];
    for i in 1:length(trimmed_w) 
        if length(findall(x->x==trimmed_w[i][2],fil_can))==0
            push!(fil_can,trimmed_w[i][2])
        end
    end
    fil_can
    candidate_den=fil_can
    #----------------------------
    exp_vec
    candidate_w
    prime_i=[103,107]
    ansatzP=[];
    ansatzC=[];
    denD=[];
    i=1
    k=1

    for l in 1:length(used_w) 
        trimmed_w=deleteat!(trimmed_w,findall(x->x[1]==used_w[l],trimmed_w)) 
    end
trimmed_w
   # trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[1],trimmed_w),trimmed_w));
  #  w=argmin(x->x[1],trimmed_w)[1]
  #  q=pickDen(R,R1,common_factors,w,candidate_den,20 );
  #  pr=q[1];
  #  pr_2=q[2];
  #  m=q[4]
  #  d=q[3]
  #  n=numeratorAnsatz(R,d);
  #  x=interpolateAnsatzTemp(R,p,n,BB,q,computed_den,computed_num)
while i<=length(prime_i) 
    p=prime_i[i];
    R1,_= QadicField(ZZ(p), 1,30);
    #trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);

   # for l in 1:length(computed_den) 
   #     trimmed_w=deleteat!(trimmed_w,findall(x->x[2]==computed_den[l],trimmed_w)) 
   # end

   # trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[1],trimmed_w),trimmed_w));
    w=argmin(x->x[3],trimmed_w)[1]
    q=pickDen(R,R1,common_factors,w,candidate_den,20 );
        pr=q[1];
        pr_2=q[2];
        m=q[4];
        d=q[3];
        n=numeratorAnsatz(R,d);
        modC=[];
        while m==0
            trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[3],trimmed_w),trimmed_w));
            w=argmin(x->x[3],trimmed_w)[1]
            q=pickDen(R,R1,common_factors,w,candidate_den,20 );
            pr=q[1];
            pr_2=q[2];
            m=q[4];
            d=q[3];
            n=numeratorAnsatz(R,d);  
        end
        x=interpolateAnsatzTemp(R,p,n,BB,q,computed_den,computed_num); 
        if x===nothing
            if k==50
                break
            else
                continue;
                k=k+1
            end
            
        else
            modC=x; 
            push!(ansatzP,n);
            push!(denD,d);
            push!(ansatzC,modC);
            i=i+1;
        end
        push!(used_w,w);
end
ansatzC
denD
#funct:computeCoef
#input: ansatzC,ansatzP,prime_i
#output: f or nothing
   f=0;
    for i in 1:length(ansatzP[1]) 
       c_i=Vector{Int64}(undef,0);
       for j in 1:length(prime_i) 
        push!(c_i,lift(ZZ,ansatzC[j][i]))
       end 
       c=manyRR(c_i,prime_i);
        f=f+R(c[1]//c[2])*ansatzP[1][i];
    end
    if f!=0
        push!(computed_num,f);
        push!(computed_den,denD[1]);
    end
 computed_num
computed_den
push!(computed_num,R(0));
        push!(computed_den,denD[1]);
#------------------------------------------------------------------------
#alternate
f=computeCoef(ansatzC,ansatzP,prime_i);
if f===nothing
    push!(computed_num,R(0));
    push!(computed_den,denD[1]);
else
    push!(computed_num,f);
        push!(computed_den,denD[1]);
end
#-----------------Test with a probe---------------------------------------
constructed_term=0;
prob1,prob2=generating_probes(R,103,common_factors,[1,0,0],rand(10000:1000000));
prob2
p=prob1[1,1]
q=prob1[2,1]
r=prob1[3,1]
                for k in 1:length(computed_den) 
                    constructed_term=constructed_term+(evaluate(computed_num[k],[p,q,r])/evaluate(computed_den[k],[p,q,r]) )
                end    
 constructed_term- evaluateBB(BB,[p,q,r])           

#--------------------------------------------------------


    i=1

 
    for i in 1:length(prime_i) 
        p=prime_i[i];
        R1,_= QadicField(ZZ(p), 1,30);
        #candidate_w=getWeights(R,R1,common_factors);
        #candidate_w=deleteat!(candidate_w,[1]);
        #    for j in 1:10 
    ##       l= getWeights(R,R1,common_factors);
    #       if length(l)>length(candidate_w)
    #        candidate_w=l;
    #       end
    #    end
        length(candidate_w)
        trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
        w=argmin(x->x[1],trimmed_w)[1]
        q=pickDen(R,R1,common_factors,w,candidate_den,20 );
        pr=q[1];
        pr_2=q[2];
        m=q[4];
        d=q[3];
        n=numeratorAnsatz(R,d);
        modC=[]
       while m==0
            trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[1],trimmed_w),trimmed_w));
            w=argmin(x->x[1],trimmed_w)[1]
            q=pickDen(R,R1,common_factors,w,candidate_den,20 );
            pr=q[1];
            pr_2=q[2];
            m=q[4];
            d=q[3];
            n=numeratorAnsatz(R,d);  
        end
    
        #while length(modC)==0
        x=interpolateAnsatzTemp(R,p,n,BB,q,computed_den,computed_num); 
        while x===nothing
        if x===nothing #&& length(trimmed_w)!=0
          #  modC=interpolateAnsatz(R,p,n,BB,q,computed_den,computed_num); 
            if length(trimmed_w)!=1
                trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[1],trimmed_w),trimmed_w));    
            end
            #trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[1],trimmed_w),trimmed_w));
            w=argmin(x->x[1],trimmed_w)[1]
            q=pickDen(R,R1,common_factors,w,candidate_den,5 );
            pr=q[1];
            pr_2=q[2];
            m=q[4];
            d=q[3];
            n=numeratorAnsatz(R,d);
            
            x=interpolateAnsatzTemp(R,p,n,BB,q,computed_den,computed_num); 
     #=   else# x===nothing && length(trimmed_w)==0
            candidate_w=getWeights(R,R1,common_factors);
            #candidate_w=deleteat!(candidate_w,[1]);
           # length(candidate_w)
            trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
            if length(trimmed_w)!=1
                trimmed_w=deleteat!(trimmed_w,findall(x->x==argmin(x->x[1],trimmed_w),trimmed_w));    
            end
            w=argmin(x->x[1],trimmed_w)[1]
            q=pickDen(R,R1,common_factors,w,candidate_den,20 );
            pr=q[1];
            pr_2=q[2];
            m=q[4];
            d=q[3];
            n=numeratorAnsatz(R,d);
            x=interpolateAnsatzTemp(R,p,n,BB,q,computed_den,computed_num); 
    #            for j in 1:10 
    #                l= getWeights(R,R1,common_factors);
    #                if length(l)>length(candidate_w)
    #                    candidate_w=l;
    #                end
    #            end
            #    trimmed_w=pickWeights(R,R1,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);
        end
        =#
        end    
        end
        modC=x; 
        #modC=interpolateAnsatz(R,p,n,BB,q,computed_den,computed_num);
        #n=numeratorAnsatz(R,d);
        push!(ansatzP,n);
        push!(denD,d);
        #modC=interpolateAnsatz(R,p,n,BB,q,computed_den,computed_num);
        push!(ansatzC,modC);
    end
ansatzC
denD
   f=0;
   oneRR(34,prime_i[1])
    for i in 1:length(ansatzP[1]) 
       c_i=Vector{Int64}(undef,0);
       for j in 1:length(prime_i) 
        push!(c_i,lift(ZZ,ansatzC[j][i]))
       end 
       c=manyRR(c_i,prime_i);
        f=f+R(c[1]//c[2])*ansatzP[1][i];
    end
    if f!=0
        push!(computed_num,f);
        push!(computed_den,denD[1]);
    end
#    push!(computed_num,f);
 #   push!(computed_den,denD[1]);

computed_num
computed_den


constructed_term=0;
prob1,prob2=generating_probes(R,103,common_factors,[1,0,0],rand(10000:1000000));
prob2
p=prob1[1,1]
q=prob1[2,1]
r=prob1[3,1]
                for k in 1:length(computed_den) 
                    constructed_term=constructed_term+(evaluate(computed_num[k],[p,q,r])/evaluate(computed_den[k],[p,q,r]) )
                end    
 constructed_term- evaluateBB(BB,[p,q,r])           

#----------------------------------
    y*(x-z)*(x+y)+(y^2 + y*z)*(x+y)+x*(x-z)*(y+z)+( 1//3*x*y + 1//3*x*z + 1//3*y^2 + 1//3*y*z)*(y+z)  
    computed_num[1]*R(x+y)+computed_num[4]*R(y+z)
    R,(x,y,z)=polynomial_ring(QQ,["x","y","z"]);
   f=R( y*(x - z)*(x + y)+(1//4*x^2 - 1//4*x*y + 3//4*x*z + 1//2*y^2 - 1//2*y*z)*(x + y)+x*(x-z)*(y+z)+(-5//4*x^2 + 5//4*y^2)*(y+z))
    t=[R((x+y)*(x-z)),R(x+y),R(x-z),R((x-z)*(y+z)),R(y+z)]
   k=0
    for i in 1:5
      k=k+computed_num[i]*t[i];  
    end
    k
   #-----------------------------------
    #Task3.2: For given prime p, find a valuatio n point x such that for all i in (1:length(w)), p^{-w_i}=|f_i(x)|_p 
p=103;
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
#funct:getWeights(R,R1,common_factors)
#output: candidate_w
#---------------------------
p=103;
R1,_= QadicField(ZZ(p), 1,30);
candidate_w=getWeights(R,R1,common_factors);
size(candidate_w)
#---------------------------
Int(prime(R1))
prob1,prob2=generating_probes(R,p,common_factors,w,rand(1000:10000))
evaluate(common_factors[1],prob1)
#-------------------------------
##Check the chosen point does the job
###Task4: Filter the candidate denominators {d_i} by using probes at the chosen points
##Evaluate R(x) at x computed above
computed_num=[];
computed_den=[];
trimmed_w=[];
trimmed_w=pickWeights(R,R1,p,BB,common_den,candidate_w,exp_vec,computed_den,computed_num);

##argmin(x->x[1],trimmed_w)[1]
##d=argmin(x->x[1],trimmed_w)[2]
##m=argmin(x->x[1],trimmed_w)[3]
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
w
q=pickDen(R,R1,common_factors,w,candidate_den,20 );
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
deg_vec=power_monomials(3,1)
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
p1=103;
interpolateAnsatz(R,p1,n,BB,q)
##We have w, d and n. We generate as many x as length(n) for satisfy (4) for a given fixed prime p
B=matrix_space(R1,length(n),1)
A=matrix_space(R1,length(n),length(n))
a=Matrix{typeof(R1(1))}(undef,length(n),length(n))
b=Matrix{typeof(R1(1))}(undef,length(n),1)
common_num
## We have w, d and n. We generate as many x as length(n) for satisfy (4) for a given fixed prime p
p1=103;
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
p=pr[1][1];
q=pr[1][2];
r=pr[1][3];
lhs=R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1))
lhs_F=coeff(lift(Qx,R1(evaluate(common_num,[p,q,r])/evaluate(f,[p,q,r]))+O(R1,p1^(m-1))),0)
    C=matrix_space(R1,1,length(n))
    C_F=matrix_space(QQ,1,length(n))
    c=Matrix{typeof(R1(1))}(undef,1,length(n))
    c_F=Matrix{typeof(QQ(1))}(undef,1,length(n))
     for i in 1:length(n) 
        c[1,i]=R1(evaluate(n[i],[p,q,r])/evaluate(d,[p,q,r]))
        c_F[1,i]=coeff(lift(Qx,R1(evaluate(n[i],[p,q,r])/evaluate(d,[p,q,r]))),0)
     end 
   c=C(c)
   c_F=C_F(c_F)
println([R1((c_F*inv(a_F)*b_F)[1])+O(R1,p1^(m-1)),R1(lhs_F)+O(R1,p1^(m-1))])



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



