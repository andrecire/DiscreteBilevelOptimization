package BilevelBOBILIB_MIXED;

import java.util.ArrayList;

/**
 * Stores a leader solution together with the corresponding DD path.
 */
public class Solution {
	int[] xbar;			// Leader variables
	Node[] pathNodes; 	// DD Nodes visited by the path
	
	/**
	 * Allocates storage for one solution vector and its path through the DD.
	 */
	public Solution(DataHandler data){
		xbar = new int[data.P.nVarL];
		pathNodes = new Node[data.P.nVarL+1];
	}
	
	
}
