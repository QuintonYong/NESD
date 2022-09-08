import it.unimi.dsi.webgraph.ImmutableGraph;
//import it.unimi.dsi.fastutil.BigArrays;
import it.unimi.dsi.fastutil.ints.IntBigArrays;

public class BigArrayGraph {
	ImmutableGraph G;
	String basename;
	int n;
    long m;
	long[] start; //size n+1
	boolean[] visited;
    int[][] dest;
	
	public BigArrayGraph(String basename) throws Exception {
		this.basename = basename;
		webgraph2array();
	}

    public BigArrayGraph(int n) {
        start = new long[n+1];
    }

	public void webgraph2array() throws Exception {
		G = ImmutableGraph.load(basename);
		n = G.numNodes();
		start = new long[n+1];
        System.out.println("Edges: " + G.numArcs());
        dest = IntBigArrays.newBigArray(G.numArcs());
		
		m = 0;
		for(int v=0; v<n; v++) {
			if(v % 1_000_000 == 0)
				System.out.println(v);
			
			start[v] = m;
			int v_deg = G.outdegree(v);
			int[] succ = G.successorArray(v);
			for(int i=0; i<v_deg; i++) {
				//dest[m] = succ[i];
                IntBigArrays.set(dest, m, succ[i]);
				m++;
			}
		}
		start[n] = m;
		System.out.println("Number of edges is " + m);
	}
	
	int numNodes() {
		return n;
	}
	
	int outdegree(int v) {
		return (int)(start[v+1] - start[v]); 
	}

    int get(long index) {
        return IntBigArrays.get(dest, index);
    }

	//Just to test how long a scan takes
	public void scanEdgeArray() {
		long cnt = 0;
		for(int i=0; i<m; i++) {
			cnt = cnt+i;
		}
		System.out.println(cnt);
	}
	
	public static void main(String[] args) throws Exception {		
		String basename = args[0];
		
		System.out.println("Starting " + basename);
		long startTime = System.currentTimeMillis();
		BigArrayGraph bigarraygraph = new BigArrayGraph(basename);
		System.out.println("webgraph2array (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);
		
		// System.out.println("Starting " + basename);
		// startTime = System.currentTimeMillis();
		// bigarraygraph.scanEdgeArray();
		// System.out.println("scanEdgeArray (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);
	}
}
