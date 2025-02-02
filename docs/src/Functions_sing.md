```@meta
CurrentModule = Feynman
```


# FUNCTIONS

```@raw html
<details>
<summary>procedure: printMat(matrix M)</summary>
```
**USAGE**   :  printMat(M); M matrix

**ASSUME**  :  M is a matrix.

**THEORY**  :  This is the print function used by Singular to print a matrix.

**KEYWORDS**: matrix

**Example** :

```singular
ring R=0,(x),lp;
matrix M[2][3]=1,243,3,4,522222,6;
printMat(M);
```
```@raw html
</details>
```

```@raw html
<details>
<summary>procedure: printGraph(graph G)</summary>

**USAGE**   :  printGraph(G); G graph

**ASSUME**  :  G is a graph.

**THEORY**  :  This is the print function used by Singular to print a graph.

**KEYWORDS**:  graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(2,4),list(3,4),list(1),list(2),list(3),list(4)));
G;
```
</details>
```

```@raw html
<details>
<summary>procedure: printLabeledGraph(labeledgraph G)</summary>

**USAGE**   :  printLabeledGraph(G); G labeledgraph

**ASSUME**  :  G is a labeled graph.

**THEORY**  :  This is the print function used by Singular to print a labeled graph.

**KEYWORDS**:  Feynman graph

**Example** :
```singular
ring R=(0),q(1..6),dp;
labeledgraph G = makeLabeledGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(1,2),list(2,4),list(3,4),list(3,4)),R, list (q(1),q(2),q(3),q(4),q(5),q(6)),R);
G;
```
</details>
```

```@raw html
<details>
<summary>procedure: printIBP(oneIBP I)</summary>

**USAGE**   :  printIBP(I); I oneIBP

**ASSUME**  :  I is an IBP identity computed using computeIBP.

**THEORY**  :  This is the print function used by Singular to print an IBP relation.

**KEYWORDS**:  Feynman graph

**Example** :
```singular

```
</details>
```

```@raw html
<details>
<summary>procedure: printsetIBP(setIBP I)</summary>

**USAGE**   :  printIBP(I); I setIBP

**ASSUME**  :  I is the set of IBP identities computed using computeIBP.

**THEORY**  :  This is the print function used by Singular to print setIBP.

**KEYWORDS**:  Feynman graph

**Example** :
```singular

```
</details>
```

```@raw html
<details>
<summary>procedure: makeGraph(list v, list e)</summary>

**USAGE**   :  makeGraph(v,e); v list, e list

**ASSUME**  :  v is a list of integers, e is a list of two element lists of v.

**RETURN**  :  graph with vertices v and edges e

**THEORY**  :  Creates a graph from a list of vertices and edges. The vertices can be any   
                type. The data structure respects the ordering of vertices of edges, so can be used for directed graphs,
**KEYWORDS**:   graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(1,2),list(2,4),list(3,4),list(3,4)));
G;
```
</details>
```

```@raw html
<details>
<summary>procedure: makeLabeledGraph(list v, list e, def R, list lab, def Rpoly)</summary>

**USAGE**   :  makeLabeledGraph(v,e,R,l,P); v list, e list, R ring, l list, P ring 

**ASSUME**  :   v is a list of integers, e is a list of two element lists of pairwise 
                different elements of v, R is a ring, l is a list of labels, P is a ring

**RETURN**  :   labeled graph with vertices v and edges e with labels of the edges in R with 
                infinite edges being constants

**THEORY**  :   

**KEYWORDS**:   Feynman graph

**Example** :
```singular
ring R=(0),q(1..6),dp;
labeledgraph G = makeLabeledGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(1,2),list(2,4),list(3,4),list(3,4)),R, list (q(1),q(2),q(3),q(4),q(5),q(6)),R);
G;
```
</details>
```

```@raw html
<details>
<summary>procedure: labelGraph(graph G, int ch)</summary>

**USAGE**   :  labelGraph(G); G graph

**ASSUME**  :  G is a graph and ch is either zero or a prime.

**RETURN**  :  labeled graph with polynomial variables q_i at the bounded edges and function 
                field variables p_i at the unbounded edges over a prime field of 
                characteristic ch

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(2,4),list(3,4),list(1),list(2),list(3),list(4)));
labeledgraph lG = labelGraph(G,0);
lG;
```
</details>
```

```@raw html
<details>
<summary>procedure: balancingIdeal(labeledgraph G)</summary>

**USAGE**   :  balancingIdeal(G); G labeledgraph

**ASSUME**  :   G is a labeled graph

**RETURN**  :  ideal of balancing condition of the graph, basering is assumed to be G.over

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(2,4),list(3,4),list(1),list(2),list(3),list(4)));
labeledgraph lG = labelGraph(G,0);
def R= lG.over;
setring R;
ideal I = balancingIdeal(lG);
I;
```
</details>
```

```@raw html
<details>
<summary>procedure: eliminateVariables(labeledgraph G)</summary>

**USAGE**   :  eliminateVariables(G); G labeledgraph

**ASSUME**  :   G is a labeled graph

**RETURN**  :  labeled graph with variables of the bounded edges eliminated according to 
                balancing condition and using an ordering $q[i]>p[j]$.

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(2,4),list(3,4),list(1),list(2),list(3),list(4)));
labeledgraph lG = labelGraph(G,0);
eliminateVariables(lG);
```
</details>
```

```@raw html
<details>
<summary>procedure: removeVariable(def R, int j)</summary>

**USAGE**   :  removeVariable(R); R ring

**ASSUME**  :   R is a polynomial ring

**RETURN**  :   polynomial ring with j-th variable removed

**THEORY**  :  

**KEYWORDS**:   ring

**Example** :
```singular
ring R=0,(x,y,z),(lp(2),dp(1));
def S= removeVariable(R,2);
S;
```
</details>
```

```@raw html
<details>
<summary>procedure: removeParameter(def R, int j)</summary>

**USAGE**   :  removeParameter(R); R ring

**ASSUME**  :   R is a polynomial ring

**RETURN**  :  polynomial ring with j-th variable removed

**THEORY**  :  

**KEYWORDS**:   ring

**Example** :
```singular
ring R=(0,p(1),p(2),p(3)),(x,y,z),(lp(2),dp(1));
def S= removeParameter(R,2);
S;
```
</details>
```

```@raw html
<details>
<summary>procedure: substituteGraph(labeledgraph G, poly a, poly b)</summary>

**USAGE**   :   substituteGraph(G); G labeledgraph

**ASSUME**  :   G is a labeled graph

**RETURN**  :   substitute the variable a in the labeling by b

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular

```
</details>
```

```@raw html
<details>
<summary>procedure: feynmanDenominators(labeledgraph G)</summary>

**USAGE**   :  feynmanDenominators(G); G labeledgraph

**ASSUME**  :   G is a labeled graph

**RETURN**  :   ideal containing the propagators in the Feynman integral

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(2,4),list(3,4),list(1),list(2),list(3),list(4)));
labeledgraph lG = labelGraph(G,0);
labeledgraph lGelim = eliminateVariables(lG);
def R = lGelim.over;
setring R;
ideal I = feynmanDenominators(lGelim);
I;
```
</details>
```

```@raw html
<details>
<summary>procedure: propagators(labeledgraph G)</summary>

**USAGE**   :   propagators(G); G labeledgraph

**ASSUME**  :   G is a labeled graph

**RETURN**  :   ideal, containing the denominators in the Feynman integral

**THEORY**  :  

**KEYWORDS**:   Feynman graph 

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4),list(list(1,3),list(1,2),list(2,4),list(3,4),list(1),list(2),list(3),list(4)));
labeledgraph lG = labelGraph(G,0);
labeledgraph lGelim = eliminateVariables(lG);
def R = lGelim.over;
setring R;
ideal I = propagators(lGelim);
I;
```
</details>
```

```@raw html
<details>
<summary>procedure: ISP(labeledgraph G)</summary>

**USAGE**   :  ISP(G); G labeledgraph

**ASSUME**  :   G is a labeled graph

**RETURN**  :  ideal, containing the irreducible scalar products, that is, those scalar 
                product which are not linearly dependent on the propagators.

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(1,2),list(3,6),list(4,5),list(1,6),list(2,3),list(5,6),list(3,4),list(1),list(2),list(5),list(4)));
labeledgraph lG = labelGraph(G,0);
labeledgraph G1 = eliminateVariables(lG);
G1;
ring R= G1.over;
setring R;
R;
ISP(G1);
```
</details>
```

```@raw html
<details>
<summary>procedure: removeElimVars(labeledgraph G)</summary>

**USAGE**   :  removeElimVars(G); G labeledgraph

**ASSUME**  :   G is a labeled graph

**RETURN**  :  Removes the variables from G.elimvars. This key is generated by the procedure 
                eliminateVariables.

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(1,2),list(3,6),list(4,5),list(1,6),list(2,3),list(5,6),list(3,4),list(1),list(2),list(5),list(4)));
labeledgraph lG = labelGraph(G,0);
labeledgraph G1 = eliminateVariables(lG);
labeledgraph G2 = removeElimVars(G1);
G2;
ring R= G2.over;
setring R;
R;
G2;
```
</details>
```

```@raw html
<details>
<summary>procedure: computeBaikovMatrix(def G0)</summary>

**USAGE**   :  computeBaikovMatrix(G); G labeledgraph, or G graph

**ASSUME**  :   G is a Graph, or
                G is a labeled graph where redundant variables have been eliminated by 
                the procedure eliminateVariables, and deleted from the ring by the 
                procedure removeElimVars.

**RETURN**  :   a labeled graph G1, computes the Baikov matrix of G defined in G1.baikovover 
                and stores it in G1.baikovmatrix

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :   
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
ring RB= G1.baikovover;
setring RB;
RB;
matrix B = G1.baikovmatrix;
printMat(B);
```
</details>
```

```@raw html
<details>
<summary>procedure: computeM1(def G0)</summary>

**USAGE**   :  computeM1(G0); G labeledgraph, or G graph

**ASSUME**  :   G is a Graph, or
                G is a labeled graph where redundant variables have been eliminated by 
                the procedure eliminateVariables, and deleted from the ring by the 
                procedure removeElimVars.

**RETURN**  :   The module M1 over G1.baikovover that requires to compute IBP identities 

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
ring RB=G1.baikovover;
RB;
module ML=computeM1(G1);
ML;
```
</details>
```

```@raw html
<details>
<summary>procedure: computeM2(def G0,list Nu)</summary>

**USAGE**   :  computeM2(G,Nu); G labeledgraph, or G graph

**ASSUME**  :   G is a Graph, or
                G is a labeled graph where redundant variables have been eliminated by 
                the procedure eliminateVariables, and deleted from the ring by the 
                procedure removeElimVars.
                Nu is the seed.

**RETURN**  :   The module M2 over G1.baikovover that requires to compute IBP identities  

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
ring RB=G1.baikovover;
RB;
module M2=computeM2(G1,list(1,1,1,0,0,1,0,0,0));
M2;
module M2=computeM2(G1, list(1,1,1,1,1,1,1,-5,0));
M2;
```
</details>
```

```@raw html
<details>
<summary>procedure: computeIBP(def G0,list Nu)</summary>

**USAGE**   :  computeIBP(G0,Nu); G labeledgraph, or G graph

**ASSUME**  :   G is a Graph, or
                G is a labeled graph where redundant variables have been eliminated by 
                the procedure eliminateVariables, and deleted from the ring by the 
                procedure removeElimVars.
                Nu is the seed.

**RETURN**  :   The set of IBPS correspond to G0 and given Nu.

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :

```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
setIBP S=computeIBP(G1,list(1,1,0,1,0,1,0,1,0));
ring R=S.over;
setring R;
S;
oneIBP I=S.IBP[1];
I;
```
</details>
```

```@raw html
<details>
<summary>procedure: getSector(list l)</summary>

**USAGE**   :  getSector(l); l list

**ASSUME**  :   l is a list of integer indices of a Feynman integral

**RETURN**  :   list L, L[1]=s The sector (a list of 1s and 0s) that the corresponding      
                integral belongs L[2]=n The sector in that the integral belongs 

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
list l=list(1,2,-3,-4,0,1);
list s=getSector(l);
s;
```
</details>
```

```@raw html
<details>
<summary>procedure: listCombinations(list L,int r)</summary>

**USAGE**   :  listCombintions(L,r); L list, r int

**ASSUME**  :   r is a positive integer such that r < size(L)

**RETURN**  :   list of r-combinations of the elements in the list L

**THEORY**  :  

**KEYWORDS**:   feynman graph

**Example** :
```singular
ring R=0,(x,y,z),dp;
list L=listCombinations(list(1,2,3,4),3);
L[1];
```
</details>
```

```@raw html
<details>
<summary>procedure: generateWebSectors(list seed)</summary>

**USAGE**   :   generateWebSectors(seed);seed list  

**ASSUME**  :   seed is a list of integer values.

**RETURN**  :   Web structure of the sectors L, where L is the list and L[1] is the sector  
                that correspond to the  given seed and L[i] contain the subsectors of the 
                sectors in L[i-1]. Note that sector maps between the sectors have not been 
                setted. 

**THEORY**  :  

**KEYWORDS**:   feynman graph

**Example** :
```singular
ring R=0,(x,y,z),dp;
list l=list(1,-1,0,1,2,-2);
list w=generateWebSectors(l);
```
</details>
```

```@raw html
<details>
<summary>procedure: isSubList(list l1,list l2)</summary>

**USAGE**   :  isSubList(l1,l2); l1 list, l2 list

**ASSUME**  :  l1 and l2 are list of positive integers

**RETURN**  :  1 if elements in l1 contain in l2
                0 if elements in l1 do not contain in l2

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
ring R=0,(x,y,z),dp;
list l1=list(1,2,3,4,5,6,7);
list l2=list(1,4,6);
list l3=list(1,2,8);
list l4=list(1,4,6);
isSubList(l2,l1);
isSubList(l3,l1);
isSubList(l1,l2);
isSubList(l2,l4);
```
</details>
```

```@raw html
<details>
<summary>procedure: getSectorMap(list L1,list L2)</summary>

**USAGE**   :  getSectorMap(L1,L2); L1 list, L2 list, sector

**ASSUME**  :   L1 and L2 are lists of sectors where the lab field of each sector in both   
                lists are filled(i.e. two layers of a sector web)

**RETURN**  :   L1 where sectorMap  of each sector in the list L1 is filled.

**THEORY**  :  

**KEYWORDS**:   sector,graph,feynman,setIBP

**Example** :
```singular
ring R=0,(x,y,z),dp;
list l=list(1,-1,0,1,2,-2);
list w=generateWebSectors(l);
list w1=getSectorMap(w[1],w[2]);
w1[1].sectorMap;
list w2=getSectorMap(w[2],w[3]);
w2[2].sectorMap;
```
</details>
```

```@raw html
<details>
<summary>procedure: setSectorMap(list sectorWeb)</summary>

**USAGE**   :  setSectorMap(sectorWeb); sectorWeb list, sector

**ASSUME**  :   sectorWeb is an output produced by the function @*generateWebSectors

**RETURN**  :  sectorWeb where the field sectorMap field of each sector in sectorWeb is 
                filled.

**THEORY**  :  

**KEYWORDS**:   sector, generateWebSectors, getSectorMap

**Example** :
```singular
ring R=0,(x,y,z),dp;
list l=list(1,-1,0,1,2,-2);
list w=generateWebSectors(l);
list w1=setSectorMap(w);
```
</details>
```

```@raw html
<details>
<summary>procedure: findSector(list sectorWeb, list currentPosition, list L)</summary>

**USAGE**   :  findSector(sectorWeb,currentPosition,L); sectorWeb list,currentPosition list,
                L list,

**ASSUME**  :   sectorWeb is an output produced by the function generateWebSectors, L is 
                an output produced by the function getSector@

**RETURN**  :   position of the sector in the sectorWeb, where the L belongs. 
                -1, if the sector is not found

**THEORY**  :  

**KEYWORDS**:   sector, generateWebSectors, getSectorMap

**Example** :
```singular
ring R=0,(x,y,z),dp;
list l=list(1,-1,0,1,2,-2);
list w=generateWebSectors(l);
list w1=setSectorMap(w);
list oneInt=list(4,-1,-1,0,-1,-2);
list L=getSector(oneInt);
def pos=findSector(w1,list(1,1),L[2]);
isSubList(w[pos[1]][pos[2]].lab,L[2])==1 && isSubList(L[2],w[pos[1]][pos[2]].lab); //this returns 1, since the given integral is in the sector at pos.

//example for a integral that is not in the set
list oneInt=list(4,1,0,-1,-2,3);
list L=getSector(oneInt);
def pos=findSector(w1,list(1,1),L[2]);
pos;
```
</details>
```

```@raw html
<details>
<summary>procedure: updateOneSector(list sectorWeb, list currentPosition,list oneInt)</summary>

**USAGE**   :   updateOneSector(sectorWeb,currentPosition,oneInt); sectorWeb list, sector 

**ASSUME**  :   sectorWeb is an output produced by the function generateWebSectors, oneInt 
                is a list of indeces of the denominators associated to an integral 
                correspond to a given feynman graph. Also assume the sectorweb isalso  
                associated to the same feynman graph.

**RETURN**  :   updated sectorWeb, where the oneInt is assigned to the targetInts field of 
                the seector correspond to provided oneInt

**THEORY**  :  

**KEYWORDS**:   sector, generateWebSectors, getSectorMap,updateWeb,findSector

**Example** :
```singular
ring R=0,(x,y,z),dp;
list l=list(1,-1,0,1,2,-2);
list w=generateWebSectors(l);
list w1=setSectorMap(w);
list oneInt=list(4,-1,-1,0,-1,-2);
list w2=updateOneSector(w1,list(1,1),oneInt);
list L=getSector(oneInt);
L[2];
w2[3][2].lab;
```
</details>
```

```@raw html
<details>
<summary>procedure: updateWeb(list sectorWeb, list currentPosition,list setInt)</summary>

**USAGE**   :   updateWeb(sectorWeb,currentPosition,setInt); sectorWeb list, sector

**ASSUME**  :   sectorWeb is an output produced by the function generateWebSectors, setInt 
                is a list of indeces of the denominators associated to  integrals correspond 
                to a given feynman graph. Also assume the sectorweb is also  associated to 
                the same feynman graph.

**RETURN**  :   list (sectorWeb,MasterInt,notInWeb) where,
                sectorWeb is the updated web by assingning integrals to correspondng sectors,
                masterInt is the list integrals belong to the sector at currentPosition
                notInWeb is the list of integrals that are not belong the integral family 
                associated the SectorWeb.

**THEORY**  :  

**KEYWORDS**:   generateWebSectors, getSector,findSector

**Example** :

**Example 1:**
```singular
ring R=(0,(t,D)),(x,y,z),dp;
list l=list(1,2,1);
list w=generateWebSectors(l);
list w1=setSectorMap(w);
list setInt=list(list(1,2,3),list(-1,1,2),list(1,1,-1),list(-1,0,-2));
list setInt=list(list(1,2,3));
list setInt=list(l);
  
list L1=pickHighestSector(setInt);
list w2=updateWeb(w1,list(1,1),L1[1]);
w2[2]; //master integrals
w2[3];//integrals not in the web
```
**Example 2:**
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
ring RZ= G1.baikovover;
printMat(G1.baikovmatrix);
  
list setInt=list(list(1,1,1,-1,-3,1,-1,-1,-1),list(1,-1,1,-1,-3,-1,-1,-4,-1));
list web=generateWebSectors(setInt[1]);
list w1=setSectorMap(web); 
web=w1;
list L1=pickHighestSector(setInt);  
  
list w2=updateWeb(web,list(1,1),L1[1]); //updateWeb returns a list w3 with w3[1]=sectorWeb,w3[2]=list of master Integrals, w3[3]=list of integrals that not belong to the current web
web=w2[1]; 
setIBP S=computeIBP(G1,L1[1][1]);
ring R=S.over;
setring R;
list L=getRedIBPs(S,101); //L[1]=list of independent IBPs,L[2]=list of master integrals
list independIBPs=L[1];
list masterAndTailIntgrals=L[2];
size(independIBPs) < size(S.IBP); //number of linearly independent set of IBPs are less than the number of orginal IBPs. So this returns true
  
oneIBP I1=independIBPs[18];     //Here is an example for one IBP i
I1;
list w3=updateWeb(web,list(1,1),masterAndTailIntgrals); //updateWeb returns a list w3 with w3[1]=sectorWeb,w3[2]=list of master Integrals, w3[3]=list of integrals that not belong to the current web
web=w3[1];   
size(web[1][1].targetInts);
```
</details>
```

```@raw html
<details>
<summary>procedure: getHighestSectorIndex(list targetInt)</summary>

**USAGE**   :   pickHighestSector(targetInt); G is a list of list of integers of same length 

**ASSUME**  :   targetInt is the list of target integrals

**RETURN**  :   the intgral that belong to the heighest sector, if all integrals belong to 
                the same sector web; otherwise, it returns a list of collection of integrals 
                each need to be handled using different sector webs,

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular

```
</details>
```

```@raw html
<details>
<summary>procedure: pickHighestSector(list targetInt)</summary>

**USAGE**   :   pickHighestSector(targetInt); G is a list of list of integers of same length

**ASSUME**  :   targetInt is the list of target integrals 

**RETURN**  :   the intgral that belong to the heighest sector, if all integrals belong to              
                the same sector web; otherwise, it returns a list of collection of integrals 
                each need to be handled using different sector webs,

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
setInt=list(list(-1,1,2),list(1,1,-1),list(-1,0,-2),list(1,2,3)); //here we can do the reduction using one web
list L=pickHighestSector(setInt);
size(L);

list setInt=list(list(-1,1,2),list(1,1,-1),list(-1,0,-2)); //here we need more than one web
list L=pickHighestSector(setInt);
size(L);
```
</details>
```

```@raw html
<details>
<summary>procedure: getSortMeasures(list l)</summary>

**USAGE**   :  getSortMeasures(l); l list, 

**ASSUME**  :   l is a list of integers (i.e a seed). 

**RETURN**  :   list of sort measures that are used in Laporta Algorithm

**THEORY**  :  

**KEYWORDS**:   Feynman graph

**Example** :
```singular
setInt=list(list(-1,1,2),list(1,1,-1),list(-1,0,-2),list(1,2,3)); 
getSortMeasures(l);
```
</details>
```

```@raw html
<details>
<summary>procedure: extractCoef(oneIBP I,list ind,list l)</summary>

**USAGE**   :  extractCoef(I,ind,l); I oneIBP,ind list,l list,

**ASSUME**  :   ind is the output of getSortedIntegrals, and l is the list of values over 
                the base field I.baikovover and size(l)=npars(I.baikovover)

**RETURN**  :  list of values where, the i-th element is the evaluation of coefficient 
                function  at values in the list l of the IBP relation oneIBP, whose index is 
                i=ind[i][1].

**THEORY**  :  

**KEYWORDS**:   feynman graph,IBPs

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
setIBP S=computeIBP(G1,list(1,1,0,1,0,1,0,1,0));
ring R=S.over;
setring R;
list ind = getSortedIntegrals(S);
oneIBP I=S.IBP[1];
I;
list rowCorrespondToI=extractCoef(I,ind,list(1,2,9)); 
  
// the nonzero coefficient of the IBP relation correspond to integral I(1,1,0,1,0,0,0,1,0).
// when we use lex ordering to order the used integrals in set of IBPs, this integral correspond to the 82th place.
// so we get only a nonzero value at position 82 and the below, the output will be -14.

rowCorrespondToI[82]; //output will be -14
```
</details>
```

```@raw html
<details>
<summary>procedure: setMat(setIBP S,list val, list ind)</summary>

**USAGE**   :  setMat(S,val); S setIBP,val list

**ASSUME**  :  size(val)=npars(I.baikovover) and val list of integers and  ind is the output 
                of getSortedIntegrals(S)

**RETURN**  :  atrix,where i-th row correspond to the evaluation of coefficient functions of 
                i-th IBP in setIBP. Columns of the matrix correspond to the all used indices 
                in the setIBP which are ordered with respect to the output 
                ofgetSortMeasures. 

**THEORY**  :  

**KEYWORDS**:   feynman graph,IBPs

**Example** :

```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
setIBP S=computeIBP(G1,list(1,1,0,1,0,1,0,1,0));
ring R=S.over;
setring R;
list ind = getSortedIntegrals(S);
matrix N=setMat(S,list(1,2,3),ind);
```
</details>
```

```@raw html

<details>
<summary>procedure: getRedIBPs(setIBP S,int p)</summary>

**USAGE**   :   getRedIBPs(S,p); 

**ASSUME**  :   S is setIBP, and p is a prime number. 

**RETURN**  :   list L, L[1]=indIBP, L[2]=seed where,
                indIBP contain the linearly independent IBP relations of setIBP which are 
                obtained by finite field row reduction over the field Fp. 
                seed contain the indeces correspond to the non-free columns in rref.

**THEORY**  :  

**KEYWORDS**:   feynman graph,IBPs

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
setIBP S=computeIBP(G1,list(1,1,0,1,0,1,0,1,0));
ring R=S.over;
setring R;
list L=getRedIBPs(S,101);
size(L[1])<size(S.IBP);
```
</details>
```

```@raw html
<details>
<summary>procedure: getSortedIntegrals(setIBP I)</summary>

**USAGE**   :  getSortedIntegrals(I); I setIBP,

**ASSUME**  :

**RETURN**  :  list ind where each entry is a pair (indv,sortmeasures),
                indv is the list of indices(seed) appered in the setIBP 
                and sortmeasures is the output of getSortMeasures(indv).
                The function getSortedIntegrals extract the seeds appeared in the IBP 
                identities of the setIBP,
                sort them lexicographically based on the values got from getSortMeasures and 
                return the output.

**THEORY**  :  

**KEYWORDS**: 

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);
setIBP S=computeIBP(G1,list(1,1,0,1,0,1,0,1,0));
ring R=S.over;
setring R;
list L=getSortedIntegrals(S); //L list of pair of sorted integrals and the corresponding sorting measures
L[1];
```
</details>
```

```@raw html
<details>
<summary>procedure: computeManyIBP(def G0,list setNu)</summary>

**USAGE**   :  computeManyIBP(G0,setNu); G0 graph,

**ASSUME**  :   setNu is a list of seeds correspond to the graph G0 which are belong to the 
                same sector 

**RETURN**  :   setIBP S, where it contains all the IBP relations obtained by module    
                intersection and seeding 

**THEORY**  :  

**KEYWORDS**: 

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
labeledgraph G1=computeBaikovMatrix(G);

//here we compute set of IBPs correspond to two seeds seperately
setIBP IBP1=computeIBP(G1,list(1,1,0,1,0,1,0,-1,0));
setIBP IBP2=computeIBP(G1,list(1,1,0,1,0,1,0,-3,0));
size(IBP1.IBP);
size(IBP2.IBP);
  
//here we compute set of IBPs correspond both seeds simultaneously 
//We can use this only when both integrals belongs to the same sector

setIBP S=computeManyIBP(G,list(list(1,1,0,1,0,1,0,-1,0),list(1,1,0,1,0,1,0,-3,0)));
size(S.IBP); 
```
</details>
```

```@raw html
<details>
<summary>procedure: getReducedIBPSystem(graph G,list targetInt )</summary>

**USAGE**   :  getReducedIBPSystem(G,targetInt); targetInt list,G graph,

**ASSUME**  :   G is a graph and targetInt is a list of seeds of target integrals.

**RETURN**  :  ist (reducedIBPs,MI) where  reducedIBPs::setIBP, MI::list.
                reducedIBPs contain the reduced IBP system for the target integrals
                MI contain the master integrals

**THEORY**  :  

**KEYWORDS**:   Feynman graph,IBPs

**Example** :
```singular
graph G = makeGraph(list(1,2,3,4,5,6),list(list(6,1),list(4,6),list(1,2),list(3,5),list(4,3),list(2,5),list(5,6),list(1),list(2),list(3),list(4)));
list targetInt=list(list(1,1,1,-1,-3,-1,-1,-1,-1),list(1,-1,1,-1,-3,-1,-1,-4,-1));
list finalset=getReducedIBPSystem(G,targetInt);
setIBP S=finalset[1];
ring R=S.over;
setring R;
oneIBP I=S.IBP[5];
I;
size(finalset[2]);
```
</details>
```