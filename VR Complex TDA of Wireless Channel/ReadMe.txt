Author:  Usatyuk Vasiliy 2019, L@Lcrypto.com - Moscow, Russia
Modificate Adam Cutbill source codes for vanila residual Wireless channel Topology and correlation matrix of channel TDA

Supplement source code for article 
V. Usatyuk, "Wireless Channels Topology Invariant as Mathematical Foundation of 
Neural Network Channel Estimation Transfer Learning Properties," 2020 43rd 
International Conference on Telecommunications and Signal Processing (TSP), Milan, Italy, 2020, pp. 105-111


1) Run main.m

Author: Adam Cutbill, 2014 (agc5@sfu.ca) - Vancouver, Canada
Description: This code takes a (small) set of data as an input and
computes the Betti Numbers of that data using a witness complex. 

Relevant Papers/Books:
Silva, V. De, & Carlsson, G. (2004). Topological estimation using witness complexes.
Zomorodian, A. (2010). Fast construction of the Vietoris-Rips complex.
Computers & Graphics, 34(3), 263–271. doi:10.1016/j.cag.2010.03.007
Algorithm from "Computational Topology: An Introduction, 2008 (H.Edelsbruneer & J.Harer)"
 
DISCLAIMER:
I have no affiliations with any of the paper authors, and cannot guarantee 
correctness for this code or explanations. This is simply being provided
as a starting point for those interested in the subject as I was. 
I couldn't find a simple MATLAB program that outlined the process,
so I wrote this after reading some papers about the subject. I hope it helps
others learn too. 
 
OTHER NOTES: There is no persistence included here.
For persistent homology, I'd recommend trying out PLEX and/or Perseus.
GPlot is not written by me in any way (it was also on the FileExchange, by Robert Piche, 2005).  

BACKGROUND
SIMPLICIAL COMPLEXES:
A simplicial complex is an n-dimensional triangular structure that, in
TDA, is made up from the original data. After this, the complex can be
used for analysis of homology and other topological features. The complex
is made up of many connected simplices.

k in a k-simplex specifies it's dimensionality.
For example, a 0-Simplex is simply a point. 2 connected points make an
1-simplex (edge). 3 connected points form a 2-Simplex (a triangle). 
4 fully connected points form a 3-Simplex (tetrahedron). 
n Connected points form an n-1 Simplex.
The dimensionality of the simplicial complex can be much higher than the
original data, depending on how connected the data is and the type of complex.
For example, with 4 2-D points you can form a tetrahedron.



Connectivity is determined by the type of complex building method used. 
The connections are put into a matrix which represents the connectivity.
The simplest complex to understand is likely the Rips complex, which 
connects points based on their Euclidean distance to each other. With
edges formed, the higher order complexes can be built up.

For example: If {1,2}, {1,3},{2,3} {3,4} are 1-simplices, and I want to
build a 2-simplices, I can see that {1,2,3} form a triangle, therefore {1,2,3}
is added to my complex as a 2-Simplex. {3,4} doesn't make a 2-Simplex with
anything because 4 is not connected. Algorithms to do this
efficiently are presented in Zomorodian's paper.

In this code, the witness complex is used to build the edges and a rips
complex algorithm builds the higher order complexes.

NOTE: Not all methods use this bottom up approach.
 
FILLING IN simplices of a k-complex:
The order of the complex built can be chosen (parameter k). No complexes
will be built beyond "k". This saves on time and memory as highly connected
data will create many high dimensional simplices.

If k=1: Edges are connected, but 2-simplices are not filled in.
If k=2: Triangles formed by the edges are now "filled in", but the
volume of 3-Complexes (tetrahedrons/ 4 point connections) are still empty.
If k=3: 3-D Tetrahedrons are filled in (4D Tetrahedrons/5 point
connections are not).
If k=n: The n space is "filled in", but n+1 is empty.

It is important to know if the space is considered "filled-in" or not
from a topological perspective. Filled in areas can be "collapsed" in
topology, but holes can't simply be removed.
     
HOW TO READ BETTI NUMBERS: The betti numbers represent the number of n-D
holes in a space. They are an important topological feature of the complex. 
Betti numbers are determined by analyzing the simplicial complex.
Like simplices, they come in different dimensionalities.
 
Betti 0: Number of connected components. If everything is connected by an
edge, betti0 is 1. If there are two groups of connected edges, betti0 is 2.
If n points are not connected (discrete), Betti0 is n.
 
Betti 1: Number of holes in a surface. For example, for a 1-Complex (k=1) a 
triangle has 3 edges and 1 hole. Therefore, each triangle in the complex
creates a hole.
 
However, if the 2-Complex is used, triangles are considered "filled in", 
meaning there will be no more hole (Betti1=0). However there may be voids
in the volume of a tetrahedron (Betti2).

Betti 2: Voids (Empty volumes).

NOTE: It doesn't always make sense to consider the higher order betti
numbers, depending on the data. For example if the data is known to
represent a 2D surface, Betti-2 to Betti-n may not be of value.

Betti Number Example- Consider noisy points sampled from the boundary of an oval shape:
 .   . .

.        .

 .   . .

Connecting nearby points (e.g. using a rips or witness complex), built from
2-simplices, you would get something like the cases below (depending on the 
threshold for connection). Case A will yield a Betti0=1 and a Betti1=1. Notice,
there is one connected component and one 2D-hole (the inside of the circle). In
Case B, the treshhold (Parameter R) is lower, increasing connectivity.
There are now filled in triangles on the right of the oval.  However, the main hole remains.
Therefore, Betti0=1 and Betti=1; those triangles can be collapsed. Essentially,
this creates a smaller hole, but the overall topology hasn't changed. The
fact that this data is an oval will hold over a range of R ( or epsilon
for a rips complex), before it is too connected and seen as a point
(betti0=1 betti1=0).

CASE A: Low connectivity.   CASE B: Higher connectivity.
 .____._.                     b .____._.
/        \                     /     \ |
.         .     <--->         ./      \|.
\        /                    |       /| 
 .____._.                     .\_____._|.
 