# Efficient PageRank

Author: Shane Peelar (peelar@uwindsor.ca)

If you choose to use this tool in your project, please cite me in your report.  You can use this BiBTeX entry to conveniently do so:

~~~
@article{peelar, title={Efficient PageRank}, url={https://github.com/InBetweenNames/PageRank}, author={Peelar, Shane M}} 
~~~

This code is intended for the 60-538 course Information Retrieval at the University of Windsor.  A Visual Studio 2017 solution is provided for the code.

~~~
	Usage: ./PageRank <association list file> (--dense | --sparse) (-t <n>)
~~~

The association list file a file that describes the directed graph composing the citation network.  It has the form:
~~~
<paperid1> <paperid2>
~~~

Where `<paperid1>` cites `<paperid2>` in the citation graph. (That is, the ordered pair `(paperid1, paperid2)` exists in the directed graph)

If `--dense` is provided, then a dense matrix version of PageRank will be called to compute the ranking.  Be warned this does not scale for large graphs.
In particular, memory usage can go through the roof very quickly!

If `--sparse` is provided, an efficient sparse matrix version of PageRank will be used to compute the ranking.  It can handle much bigger sizes.

The `-t` argument is the "teleport" chance that the surfer will hop to a random point in the graph.