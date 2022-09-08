# NESD

## How to run it? 
1. Navigate to the directory Java files exist.
2. Compile all files in the following way: (the required jar files are in the lib directory)

javac -cp "lib/*" -d bin *.java

3. Run NESD.class with the following input parameters (replace : with ; if on Windows instead of Linux):

java -Xms[required memory]g -cp "bin:lib/*" NESD [input graph] [epsilon] [c]

## Preparing Data

There are three files in this format: 

*basename.graph* <br>
*basename.properties* <br>
*basename.offsets*


Let us see for an example dataset, *cnr-2000*, in 
http://law.di.unimi.it/webdata/cnr-2000

There you can see the following files available for download.

*cnr-2000.graph* <br>
*cnr-2000.properties* <br>
*cnr-2000-t.graph* <br>
*cnr-2000-t.properties* <br>
*...* <br>
(you can ignore the rest of the files)

The first two files are for the forward (regular) *cnr-2000* graph. The other two are for the transpose (inverse) graph. If you only need the forward graph, just download: 

*cnr-2000.graph* <br>
*cnr-2000.properties*

What's missing is the "offsets" file. This can be easily created by running:

__java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L cnr-2000__


### Symmetrizing graph
In order to obtain undirected graphs, for each edge we add its inverse. This can be achieved by taking the union of the graph with its *transpose*. Here we show how to do this for cnr-2000.

Download from http://law.di.unimi.it/datasets.php:

*cnr-2000.graph* <br>
*cnr-2000.properties* <br>
*cnr-2000-t.graph* <br>
*cnr-2000-t.properties*

(The last two files are for the transpose graph.)

Build the offsets:

__java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L cnr-2000__ <br>
__java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -o -O -L cnr-2000-t__

Symmetrize by taking union:

__java -cp "lib/*" it.unimi.dsi.webgraph.Transform union cnr-2000 cnr-2000-t cnr-2000-sym__