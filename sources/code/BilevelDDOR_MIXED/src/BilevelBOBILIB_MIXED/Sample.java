package BilevelBOBILIB_MIXED;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

/**
 * Represents one sampled follower response used to generate blocking cuts.
 */
public class Sample {
	double folObj;
	double[] y;						//Follower solution
	double[] yopt;					//Follower optimal solution in response to a given x
	ArrayList<double[]> block; 		// Blocking description: [constraint index, blocking value]
	String key;

	/**
	 * Builds a sample from a follower solution and precomputes the associated
	 * blocking information.
	 */
	public Sample(DataHandler data, double[] yVal, double obj){
		folObj = obj;	
		block = new ArrayList<>();
		y =  new double[data.P.nVarF];
		double [] A = new double[data.P.nConsF];
		for (int i = 0; i < data.P.nVarF; i++) {
			if(!data.isInterdiction) {
				y[i] = yVal[i];
			}
			else if(data.isInterdiction) {
				if(yVal[i]>0.5) {y[i]=1;}
			}
			for (int j = 0; j < data.P.nConsF; j++) {
				A[j] += AlgorithmHandler.round(yVal[i]*data.P.Afol[j][i]);
			}
		}
		//System.out.println("SY: "+Arrays.toString(A));
		//System.out.println("bfol: "+Arrays.toString(data.bfol));
		key = folObj+" ";
		for (int q = 0; q < data.P.nConsF; q++) {
			if(data.cInteraction[q].equals("B")) { //We care about B constraints 
				double[] blockq = new double[2];
				blockq[0] = q;
				blockq[1] = data.P.bfol[q]-A[q]+1e-4;			// Using an epsilon because of the continuous data
				if(data.isInterdiction) {
					blockq[1] = data.P.bfol[q]-A[q]+1;			// Using epsilon=1 because of the integer data in interdiction constraints
				}
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
