package BilevelBOBILIB_DISCRETE;

import java.util.ArrayList;

/**
 * Arc of the relaxed decision diagram.
 *
 * <p>Each arc stores its endpoints, the binary decision taken on the current
 * leader variable, and the propagated cost contribution used during the
 * backward dynamic-programming phase.</p>
 */
public class Arc {
	int ID; 		//Position in the arc array
	int tail;
	int head;
	int arcVal;
	int varIndex;
	boolean leaderVar ;
	int cost;
	String key;
	
	
	public Arc(int _tail, int _head, int _arcVal, int _varIndex, boolean leader){
		key="";
		tail = _tail;
		head = _head;
		arcVal = _arcVal;
		varIndex = _varIndex;
		leaderVar = leader ;
		cost = 0;
	}

	/**
	 * Builds a human-readable identifier from the tail and head nodes.
	 */
	public String getKey(BDDGraph bdd) {
		key = bdd.nodes.get(tail).key+" -> "+bdd.nodes.get(head).key;
		return key;
	}
	

	/**
	 * Prints the cached arc identifier.
	 */
	public void print() {
		
		System.out.println(key);
				
	}
	
	
}
