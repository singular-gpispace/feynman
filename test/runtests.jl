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

