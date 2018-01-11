# truss-decomposition
In-memory truss-decomposition algorithm from the PVLDB'12 paper "Truss Decomposition in Massive Networks".
Please cite if you use the code.

Usage of `imtd.cpp`:

#### Input
1st line:	n m	// #vertices, #edges
(i+1)th line	u v	// ith edge (u,v)

#### Output
m lines containing:
u v c	// (u,v) belongs to c-class,
		meaning it's in c-truss but not (c+1)-truss.
