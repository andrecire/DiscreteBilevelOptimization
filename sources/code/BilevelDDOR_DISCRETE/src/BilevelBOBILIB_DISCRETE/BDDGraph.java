package BilevelBOBILIB_DISCRETE;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;

import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

/**
 * In-memory representation of the relaxed binary decision diagram used by the
 * discrete solver.
 */
public class BDDGraph {

	ArrayList<Node> nodes;
	ArrayList<Arc> arcs;
	ArrayList<Node>[] layers;
	
	int lastNode; //Points to the index of the last node in the nodes array
	List<Integer> varOrder;
	// Map containing the nodes already created
	Hashtable<String, Integer> GraphMap;
	
	double ub;
		
	public BDDGraph() {
		nodes = new ArrayList<Node>();
		arcs = new ArrayList<Arc>();
		GraphMap = new Hashtable<String, Integer>();
		ub = Double.POSITIVE_INFINITY;
	}
	
	/**
	 * Prints a detailed view of the current diagram, including outgoing arcs for
	 * every node.
	 */
	public void print(DataHandler data) {
		
		System.out.println("***************************DP GRAPH*****************************************************************");
		System.out.println(" NODES: "+nodes.size()+" arcs: "+arcs.size());
		for (int i = 0; i < nodes.size(); i++) {
			System.out.println("Node "+i+": "+nodes.get(i).key);
			for (int j = 0; j < nodes.get(i).out.size(); j++) {
				String head = nodes.get(nodes.get(i).out.get(j).head).key;
				System.out.print("      arc "+nodes.get(i).out.get(j).arcVal+" variable "+nodes.get(i).out.get(j).varIndex+": "+head);
			}
			
			System.out.println();
			
	}
		/*
		for (int k = 0; k < data.K; k++) {
			System.out.println("ARCS MAP BDD "+k);
			for (int i = 0; i < data.numArcs[k]; i++) {
				System.out.println(data.arcsToAT.get(k)[i]);	
			}
		}
		
		for (int i = 0; i < arcs.size(); i++) {
			System.out.println("Arc "+i+": "+arcs.get(i).tail+" -- "+arcs.get(i).head+" cost: "+arcs.get(i).cost);
			for (int j = 0; j < arcs.get(i).subArcs.length; j++) {
				System.out.println("     "+arcs.get(i).subArcs[j]);
			}
		}
		*/
		
	}
	/**
	 * Prints only the number of nodes and arcs currently stored in the diagram.
	 */
	public void print2(DataHandler data) {
		System.out.println(" NODES: "+nodes.size()+" arcs: "+arcs.size());
		/*for (int i = 0; i < data.nVarL+1; i++) {
			for (int j = 0; j < layers[i].size(); j++) {
				System.out.println("LAYER "+i+" : "+layers[i].get(j).key);
			}
		}
		 */
		
	}
	/**
	 * Prints every arc together with its endpoints, branch value, and cost.
	 */
	public void print3(DataHandler data) {
		System.out.println(" NODES: "+nodes.size()+" arcs: "+arcs.size());
	
		for (int i = 0; i < arcs.size(); i++) {
		//	Arc a = arcs.get(i);
		//	a.getKey(this);
		//	a.print();
			//if(arcs.get(i).cost != 0) {
			System.out.println("ARC "+i+" from "+arcs.get(i).tail+" to "+arcs.get(i).head+" arc val: "+arcs.get(i).arcVal+" cost "+arcs.get(i).cost );
			//}
		}

//		print(data);
		
	}


	/**
	 * Clears all node adjacency lists while keeping the node objects themselves.
	 */
	public void clearInOut() {
		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).in.clear();
			nodes.get(i).out.clear();
		}
		
	}

	/**
	 * Resets every arc cost to zero before a fresh backward propagation pass.
	 */
	public void resetArcCosts() {
		for (int i = 0; i < arcs.size() ; i++) {
				arcs.get(i).cost = 0;
		}
		
	}

	/**
	 * Removes all outgoing arcs from a node and disconnects them from their head
	 * nodes.
	 */
	public void removeOutArcs(int nodeIndex) {
		
		for (int i = 0; i < nodes.get(nodeIndex).out.size(); i++) {
			Arc a = nodes.get(nodeIndex).out.get(i);
			nodes.get(a.head).in.remove(a); //Remove from destination
			arcs.remove(a);	//Remove from arc list
		}
		nodes.get(nodeIndex).out.clear(); //Remove from origin
				
	}

	
	
}
