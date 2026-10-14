package BilevelBOBILIB_DISCRETE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

/**
 * Node of the relaxed decision diagram used by the discrete solver.
 *
 * <p>The node stores both the diagram state and auxiliary lower-bound
 * information required for pruning and relaxation.</p>
 */
public class Node {

	String key;             //String representation
	int ID;                 //Position in the DD nodes array
	ArrayList<Arc> out;     //Arcs going out
	ArrayList<Arc> in;      //Arcs coming in
	int layer;
	boolean removed;
	boolean cutGenerated;
	int folResponse;

	// State: current LHS
	int[] S;

	// Not part of the state:   
	int[] Lfol; //Lower bound on the leader contribution to the follower constraints
	int[] Llead; //Lower bound on the leader contribution to the leader constraints
	int minCost; //Shortest path to each node			

	// Constructs a fresh node for the root or for a newly created layer state.
	public Node(DataHandler data) {
		key = "";
		removed = false;
		cutGenerated = false;
		ID = -1;
		S = new int[data.P.nConsF];
		Lfol = new int[data.P.nConsF];
		Llead = new int[data.P.nConsL];
		out = new ArrayList<Arc>();
		in = new ArrayList<Arc>();
		layer = 0;
		folResponse = 0;
		minCost = 0;
	}
	// Copy constructor used when branching or relaxing nodes.
	Node(Node node) {
		key = "";
		out = new ArrayList<>();
		in = new ArrayList<>();
		S = new int[node.S.length];
		Lfol = new int[node.Lfol.length];
		Llead = new int[node.Llead.length];

		for (int i = 0; i < S.length; i++) {
			S[i] = node.S[i];
			Lfol[i] = node.Lfol[i];
		}
		for (int i = 0; i < Llead.length; i++) {
			Llead[i] = node.Llead[i];
		}

		minCost = node.minCost;

	}
	/**
	 * Creates the child reached when the current variable is fixed to zero.
	 */
	public Node createZeroNode() {
		Node node = new Node(this); // State after a zero arc stays the same!!!

		return node;
	}

	/**
	 * Creates the child reached when the current leader variable is fixed to one.
	 */
	public Node createOneNode(int i, DataHandler data) {
		Node node = new Node(this);
		// Update the state after one arc for variable i
		for (int q = 0; q < data.P.nConsF; q++) {
			if(data.cInteraction[q].equals("B")) { //We only keep track of the B constraints 
				node.S[q]+=data.P.Bfol[q][i];
			}
			node.Lfol[q]+=data.P.Bfol[q][i]; // We keep lower bounds for all constraints
		}
		// Leader lower bound
		for (int q = 0; q < data.P.nConsL; q++) {
			node.Llead[q]+=data.P.Aleader[q][i]; // We keep lower bounds for all constraints
		}

		node.minCost += data.P.c1leader[i];
		return node;
	}

	/**
	 * Generates the canonical key used to deduplicate nodes within a layer.
	 */
	public void genKey(int _layer) {
		layer = _layer;
		key = Arrays.toString(S);
		key += " Layer: "+layer;
	}
	/**
	 * Encodes the outgoing arcs of the node as a compact string.
	 */
	public String genOut() {
		String outKey="";
		for (int i = 0; i < out.size(); i++) {
			outKey+="arc"+out.get(i).arcVal+"+"+out.get(i).cost+"+"+out.get(i).head;
		}
		//System.out.println("OutKEY: "+outKey);
		return outKey;
	}


	/**
	 * Prints the node key.
	 */
	public void print() {
		System.out.println(key);
		System.out.println();
	}

	/**
	 * Checks whether the node can be feasibly completed using the supplied
	 * completion bounds.
	 */
	public boolean checkBounds(DataHandler data, int[] completionBoundsL, int[] completionBoundsF) {
		// Check each one of the leader constraints
		for (int q = 0; q < data.P.nConsL; q++) {
			if(Llead[q]+completionBoundsL[q] > data.P.b[q]) {
				// Debug: leader completion bound violated.
				return false;
			}
		}
		for (int q = 0; q < data.P.nConsF; q++) {
			if(data.cInteraction[q].equals("B") && Lfol[q]+completionBoundsF[q] > data.stateBounds[q]) {
				// Debug: follower completion bound violated.
				return false;
			}
		}

		return true;
	}
	/**
	 * Updates stored lower bounds and the best path cost seen for this node.
	 */
	public void updateL(DataHandler data, int[] Lfol, int[] Llead, int minC) {
		for (int q = 0; q < data.P.nConsF; q++) {
			if(Lfol[q] < this.Lfol[q]) {
				this.Lfol[q] = Lfol[q];
			}
		}
		for (int q = 0; q < data.P.nConsL; q++) {
			if(Llead[q] < this.Llead[q]) {
				//System.out.println("Updating Llead "+this.Llead[q]+" to "+Llead[q]);
				this.Llead[q] = Llead[q];
			}
		}
		if(minC < this.minCost) {
			this.minCost = minC;
		}
	}

	/**
	 * Merges two nodes into a relaxed representative.
	 */
	public Node relax(Node node) {
		Node relaxed = new Node(this); // Make a copy and the keep all the max/min component wise
		for (int q = 0; q < S.length; q++) {
			if(node.S[q] > relaxed.S[q]) {
				relaxed.S[q] = node.S[q];
			}
			if(node.Lfol[q] < relaxed.Lfol[q]) {
				relaxed.Lfol[q] = node.Lfol[q];
			}
		}
		for (int q = 0; q < Llead.length; q++) {
			if(node.Llead[q] < relaxed.Llead[q]) {
				relaxed.Llead[q] = node.Llead[q];
			}
		}

		if(node.minCost < relaxed.minCost) {
			relaxed.minCost = node.minCost;
		}
		relaxed.genKey(this.layer);
		// Point all the incoming arcs
		for (int i = 0; i < this.in.size(); i++) {
			relaxed.in.add(this.in.get(i));
		}
		for (int i = 0; i < node.in.size(); i++) {
			relaxed.in.add(node.in.get(i));
		}
		//System.out.println("Relaxing node cost "+node.minCost+" and "+this.minCost);
		return relaxed;

	}
	/**
	 * Normalizes outgoing arc costs and pushes the common offset backward through
	 * the incoming arcs.
	 */
	public void propagateCost(DataHandler data) {
		// Find min arc cost out
		int minCost = 999999999;
		for (int i = 0; i < out.size(); i++) {
			if(out.get(i).cost < minCost) {
				minCost = out.get(i).cost;
			}
		}
		//System.out.println("PROPAGATING COST AT NODE "+key); 
		// Subtract this minCost from all outgoing arcs  
		for (int i = 0; i < out.size(); i++) {
			out.get(i).cost = out.get(i).cost - minCost ;
			//System.out.println("------- COST OUT "+out.get(i).arcVal+" : "+out.get(i).cost);
		}
		// Add this minCost to all ingoing arcs 
		for (int i = 0; i < in.size(); i++) {
			in.get(i).cost = minCost;
			//System.out.println("------- COST IN "+in.get(i).cost);
		}
	}
	/**
	 * Initializes the state when follower constraints depend only on leader
	 * variables; those constraints can then be ignored in the DD state.
	 */
	public void init(DataHandler data) {
		for (int q = 0; q < data.P.nConsF; q++) {
			if(data.cInteraction[q].equals("L")) {
				S[q]+=data.P.bfol[q];
			};
		}

	}

	/**
	 * Updates the state componentwise with the maximum of the current node and a
	 * candidate node.
	 */
	public void reduceState(Node node) {
		for (int q = 0; q < S.length; q++) {
			if(node.S[q] > this.S[q]) {
				this.S[q] = node.S[q];
			}
		}

	}
	/**
	 * Returns the state restricted to linking constraints.
	 */
	public String getS(DataHandler data) {
		String out = "[";
		for (int q = 0; q < S.length; q++) {
			if(data.cInteraction[q].equals("B")) {
				out+=(S[q]+",");
			}
		}
		out+="]";
		return out;
	}


}

