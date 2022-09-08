import java.util.Random;
import java.util.Arrays;
import java.util.ArrayList;
import it.unimi.dsi.fastutil.ints.IntBigArrays;

public class NESD {
	BigArrayGraph G;
	boolean printprogress = false; 
	long E=0;
	int n; 
	int md; //max degree
    double epsilon_1;
    double epsilon;
    double c;
    double last_label;
	
	public NESD(String basename, double epsilon_1, double c) throws Exception {
		G = new BigArrayGraph(basename);		
		n = G.numNodes();
        this.epsilon_1 = epsilon_1;
        this.epsilon = epsilon_1/3;
        this.c = c;
        last_label = 0;
		
		md = 0; 
		for(int v=0; v<n; v++) {
			int v_deg = G.outdegree(v);
			if(md < v_deg) 
				md = v_deg;			
		}
	}
	
    public double[] KCoreCompute () {
        boolean max_only = true;
        boolean cointoss = true;

        long startTime = System.currentTimeMillis();

        long period_start_time = System.currentTimeMillis();
        // Sampling Alg Setup
        int num_v     = n;
        int[] V       = new int[n]; // vertices in V divided by L and H
        int[] pos     = new int[n]; // position in V of each vertex
    	int[] deg     = new int[n]; // degree of each vertex
        int[] k       = new int[n]; 
        int[] t       = new int[n]; // t(v) for each v in V
        int[] deg_bins = new int[md+1];
        double[] label   = new double[n]; // label of each v in V
        long Sample_len = 0;
        int[][] Sample;
        int[] sampled_H = new int[n];
        long[] Sample_start = new long[n];
        int[] sampled_L = new int[n];

        boolean bin_search = false;
        int pow = 1;
        int min_pow = 1;
        int max_pow = 1;
        double start_l = 0;
        double start_p = 0;
        double l = (double) n / (1+epsilon);
        double p = (double) 2 * (1 + c) * (Math.log(n) + Math.log(Math.log(n) / Math.log(1 + epsilon))) * (1 + epsilon) * (1 + epsilon) / (n + epsilon * epsilon);

        float[] rand_u = new float[md];
        Random rand = new Random();

        // Timing 
        double init_time = 0;
        double splitting_time = 0;
        double sampling_H_time = 0;
        double peeling_L_time = 0;
        double post_time = 0;
    	
    	for(int v=0; v<n; v++) { 
            V[v] = v;
            pos[v] = v;
    		deg[v] = G.outdegree(v);
            deg_bins[deg[v]]++;

            label[v] = -1;
            Sample_start[v] = Sample_len;
            Sample_len += deg[v];

    	}
        //Sample = new int[Sample_len];
        Sample = IntBigArrays.newBigArray(Sample_len);

        if (cointoss) {
            // initialize random
            for (int i = 0; i < md; i++) {
                rand_u[i] = rand.nextFloat();
            }
        }

        // determine starting threshhold
        int md_d = 0;
        int deg_greater = 0;
        for (int d = md; d >= 0; d--) {
            deg_greater += deg_bins[d];
            if (deg_greater >= d) {
                md_d = d;
                break;
            }
        }

        int test_count = 0;
        int start_thresh = md_d;
        System.out.println("Max degree: " + md);
        System.out.println("Sqrt m: " + (int)Math.sqrt(G.G.numArcs()));
        System.out.println("d nodes of deg d: " + md_d);
        while (l > start_thresh) {
            l = l / (1 + epsilon);
            p = p * (1 + epsilon);
            test_count++;
        }
        start_l = l;
        start_p = p;
        init_time = (System.currentTimeMillis() - period_start_time)/1000.0;

        test_count = 0;
        while (p < 1) {
            test_count++;
            // divide V into L and H
            int L_index = 0;
            int H_index = num_v - 1;            
            int curr_v = 0;
            int curr_v_deg = 0;

            //System.out.println(p);            
            //System.out.println("Splitting.");
            period_start_time = System.currentTimeMillis();

            int count = 0;
            while(L_index < H_index) {
                count++;
                curr_v = V[L_index];
                curr_v_deg = deg[curr_v];

                if (curr_v_deg >= l) {
                    pos[V[L_index]] = H_index;
                    pos[V[H_index]] = L_index;
                    V[L_index] = V[H_index];
                    V[H_index] = curr_v;
                    H_index--;
                } else {
                    L_index++;
                }
            }
            if (deg[V[H_index]] < l) {
                H_index++;
            } 

            int Q_index = H_index;

            //System.out.println("Done splitting. " + count);
            splitting_time += (System.currentTimeMillis() - period_start_time)/1000.0;
            
            //System.out.println("Sampling H");
            period_start_time = System.currentTimeMillis();

            // for all v in H, run Sample(v) and set t(v) <-- k(v)
            for(int i = H_index; i < num_v; i++) {
                int v = V[i];
                int kv = (int)(p * deg[v]);
                int tv = kv;
                long v_start = G.start[v];

                int sample_count = 0;
                int deg_v = deg[v];
                while(sample_count < kv) {
                    int u = 0;
                    if (cointoss) { u = G.get(v_start + (int)(deg_v * rand_u[sample_count])); }
                    else          { u = G.get(v_start + sample_count); }
                    
                    if (pos[u] < H_index) { // if u is in L
                        tv--;
                    } else {
                        //Sample[Sample_start[u] + sampled_L[u]] = v;
                        IntBigArrays.set(Sample, Sample_start[u] + sampled_L[u], v);
                        sampled_L[u]++;
                    }
                    sample_count++;
                }

                if (tv < l * kv / deg_v) {
                    int pv = pos[v];
                    V[pos[v]] = V[H_index];
                    V[H_index] = v;
                    pos[V[pv]] = pv;
                    pos[v] = H_index;
                    H_index++;
                }

                k[v] = kv;
                t[v] = tv;
            }

            //System.out.println("Done sampling H.");
            sampling_H_time += (System.currentTimeMillis() - period_start_time)/1000.0;

            //System.out.println("Peeling  L");
            period_start_time = System.currentTimeMillis();
            // peeling
            for (int i = Q_index; i < H_index; i++) { // for each u in L
                int u = V[i];

                for (int j = 0; j < sampled_L[u]; j++) { // for each v in Sample(u)
                    //int v = Sample[Sample_start[u] + j];
                    int v = IntBigArrays.get(Sample, Sample_start[u] + j);
                    
                    if(pos[v] >= H_index) { // if v is in H
                        t[v]--;
                        if (t[v] < l * k[v] / deg[v]) {
                            // move v from H to L (put it at the end of L)
                            int pv = pos[v];
                            V[pos[v]] = V[H_index];
                            V[H_index] = v;
                            pos[V[pv]] = pv;
                            pos[v] = H_index;
                            H_index++;
                        }
                    }
                }
            }
            //System.out.println("Done peeling L");
            peeling_L_time += (System.currentTimeMillis() - period_start_time)/1000.0;

            period_start_time = System.currentTimeMillis();

            if (max_only && H_index < num_v) {
                if (!bin_search) { 
                    pow = pow / 2; 
                } // since pow is next power to divide, first time found, divide by 2 to get current power
                bin_search = true;

                max_pow = pow;
                if (max_pow - min_pow > 1) {
                    pow = (max_pow + min_pow) / 2;
                    double mid_mult = Math.pow((1 + epsilon), pow);
                    l = start_l / mid_mult;
                    p = start_p * mid_mult;
                    //reset all samples
                    for (int i = 0; i < num_v; i++) { sampled_L[i] = 0; }
                    continue;
                }

                System.out.println(test_count);
                System.out.println("init time: " + init_time);
                System.out.println("splitting time: " + splitting_time);
                System.out.println("sampling H time: " + sampling_H_time);
                System.out.println("peeling L time: " + peeling_L_time);
                System.out.println("post time: " + post_time);
                System.out.println("Max core approximate core: " + l);
                System.out.println("H nodes: " + (num_v - H_index) + " / " + num_v + " (n = " + n + ")");
                return label;
            }

            // doing binary search and nothing labeled
            if (max_only && bin_search) {
                min_pow = pow;
                if (max_pow - min_pow > 1) {
                    pow = (max_pow + min_pow) / 2;
                    double mid_mult =  Math.pow((1 + epsilon), pow);
                    l = start_l / mid_mult;
                    p = start_p * mid_mult;
                    //reset all samples
                    for (int i = 0; i < num_v; i++) { sampled_L[i] = 0; }
                    continue;
                }

                System.out.println(test_count);
                System.out.println("init time: " + init_time);
                System.out.println("splitting time: " + splitting_time);
                System.out.println("sampling H time: " + sampling_H_time);
                System.out.println("peeling L time: " + peeling_L_time);
                System.out.println("post time: " + post_time);
                System.out.println("Max core approximate core: " + start_l / Math.pow((1 + epsilon), max_pow));
                System.out.println("H nodes: " + (num_v - H_index) + " / " + num_v + " (n = " + n + ")");
                return label;
            }

            // label every node in H by l
            // for (int i = H_index; i < num_v; i++) {
            //     label[V[i]] = l;
            // }

            //reset all samples
            for (int i = 0; i < num_v; i++) {
                sampled_L[i] = 0;
            }

            //num_v = H_index; // V <-- L
            // l = l / (1 + epsilon);
            // p = p * (1 + epsilon);
            double pow_value = Math.pow((1 + epsilon), pow);
            l = start_l / pow_value;
            p = start_p * pow_value;
            min_pow = pow / 2;
            pow *= 2;
            
            post_time += (System.currentTimeMillis() - period_start_time)/1000.0;

        }
        last_label = l;

        System.out.println("Last label: " + l);

        System.out.println("Done rand part " + (n - num_v) + " / " + n + " labels), Time = " + (System.currentTimeMillis() - startTime)/1000.0);
        System.out.println("init time: " + init_time);
        System.out.println("splitting time: " + splitting_time);
        System.out.println("sampling H time: " + sampling_H_time);
        System.out.println("peeling L time: " + peeling_L_time);
        System.out.println("post time: " + post_time);

        System.out.println("Starting peeling algorithm");
        startTime = System.currentTimeMillis();
    	
    	return label;
    }
    
	public static void main(String[] args) throws Exception {
		//args = new String[] {"uk-2005"};
		
		if(args.length != 3) {
			System.err.println("Usage: java NESD basename epsilon c");
			System.exit(1);
		}
		
		String basename = args[0];
        double epsilon = Double.parseDouble(args[1]);
        double c = Double.parseDouble(args[2]); 
		
		System.out.println("Starting " + basename);
		long startTime = System.currentTimeMillis();
		NESD kc = new NESD(basename, epsilon, c);
		System.out.println(args[0] + ": Time elapsed for loading (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);
		
		startTime = System.currentTimeMillis();
        kc.KCoreCompute();
		System.out.println(args[0] + ": Time elapsed for core computing (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);		
	}
}

