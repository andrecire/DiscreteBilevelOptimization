package BilevelBOBILIB_DISCRETE;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Represents one sampled follower response used to generate blocking cuts.
 */
public class Sample {
	int folObj;
	int[] y;					//Follower solution
	int[] yopt;					//Follower optimal solution in response to a given x
	ArrayList<int[]> block; 	// Blocking description: [constraint index, blocking value]
	String key;

	/**
	 * Builds a sample from a follower solution and precomputes the associated
	 * blocking information.
	 */
	public Sample(DataHandler data, int[] yVal, int obj){
		folObj = obj;	
		block = new ArrayList<>();
		y =  new int[data.P.nVarF];
		int [] A = new int[data.P.nConsF];
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = yVal[i];
			for (int j = 0; j < data.P.nConsF; j++) {
				A[j] += yVal[i]*data.P.Afol[j][i];
			}
		}
		//System.out.println("SY: "+Arrays.toString(A));
		//System.out.println("bfol: "+Arrays.toString(data.bfol));
		key = folObj+" ";
		for (int q = 0; q < data.P.nConsF; q++) {
			if(data.cInteraction[q].equals("B")) { //We care about B constraints 
				int[] blockq = new int[2];
				blockq[0] = q;
				blockq[1] = (int) (data.P.bfol[q]-A[q]+1);
				if(blockq[1] <= data.stateBounds[q]) {
					block.add(blockq);
					key+=Arrays.toString(blockq);
				}
			}
		}


	}

	/**
	 * Prints the sample in a debug-friendly format.
	 */
	public void print() {
		System.out.print("Y SAMPLE OBJ "+folObj+" block: [");
		for (int i = 0; i < block.size(); i++) {
			System.out.print(" ("+block.get(i)[0]+" "+block.get(i)[1]+")");
		}
		System.out.println("]");

	}


}
