/**
 * Core implementation of the discrete DD-based bilevel solver.
 *
 * <p>This class contains the routines used to build the relaxed decision
 * diagram, solve the associated single-level reformulations, and export final
 * solutions and diagnostics.</p>
 */

package BilevelBOBILIB_DISCRETE;
import ilog.concert.*;
import ilog.cplex.*;
import ilog.cplex.IloCplex.UnknownObjectException;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;
import java.util.Random;

import BilevelBOBILIB_DISCRETE.Node;

public class AlgorithmHandler {


	double zstarIP;
	double zstarBDD;
	double ipGap;

	BDDGraph DD;							//DD for the follower 

	int[] xstar; 							// Leader solution
	int[] ystar; 							// Follower solution
	int LB;
	int UB;
	int [] yLow;							//Bounds for var fixing
	int [] yUp;
	int initialUB;							// UB before refinement
	double LR; 								// Linear relaxation
	boolean optimal;						//Flag for relaxed DD approach

	IloCplex cplex;							// Cplex model 1
	IloCplex cplexAP;						// Cplex model approx
	IloRange[] AforAP;						// For model approx
	IloNumVar[] yAP;						// For model approx
	//	int M;									// Big-M value
	int[]  maxInteraction;					// Maximum values for xL*Bfol
	int folUB; 								// UB for the follower obj (from the DD)
	int folObjM;							// UB for the follower obj (from the HPP)
	double[] MfolConstraint;				// LB on the leader contribution to follower constraints
	int numCuts;
	int conflictCuts;
	int width;								// For the relaxation
	boolean isExact;						// If the DD is exact
	boolean preproMaxout;					// If the prepro time is over
	boolean logical;						// Turns on/off the logical cuts

	Random r;								// For diversity cuts
	int numCores;							//For CPLEX

	ArrayList<Sample> sampleMaster;			//Official sample
	ArrayList<Sample> sampleAux;			//Collect all responses here
	Hashtable<String, Integer> sampleMap;

	double timeDD;
	double timeMaxMin;

	/**
	 * Initializes the solver state for one discrete instance.
	 */
	public AlgorithmHandler(DataHandler data, int numcores, Boolean logicalCuts) throws InterruptedException, IloException {
		zstarIP = (int) Double.POSITIVE_INFINITY;
		zstarBDD = (int) Double.POSITIVE_INFINITY;
		LB = -999999999;
		UB = 999999999;
		folUB = -999999999;
		xstar = new int[data.P.nVarL];
		ystar = new int[data.P.nVarF];
		maxInteraction = new int[data.P.nConsF];
		numCuts = 0;
		conflictCuts = 0;
		optimal = false;
		LR = -1;
		sampleMaster = new ArrayList<Sample>();
		sampleAux = new ArrayList<Sample>();
		sampleMap = new Hashtable<String, Integer>();
		isExact = true;
		cplexAP = null;
		preproMaxout = false;
		numCores = numcores;
		xstar = new int[data.P.nVarL];
		ystar = new int[data.P.nVarF];
		logical = logicalCuts;
	}

	/**
	 * Builds the relaxed decision diagram used by the downstream optimization
	 * procedures.
	 */
	public void genRelaxedDD(DataHandler data, int setWidth) throws IloException {
		double timeDDstart = System.currentTimeMillis();
		// Create DD object
		DD = new BDDGraph();
		width = data.calcWidth(setWidth);
		//Run heuristic for variable ordering
		List<Integer> varOrderL = new ArrayList<>();
		orderingHeuristicL(varOrderL,data);
		System.out.println("Variable Order Leader: " + varOrderL);
		DD = buildBDDRelaxed(data, varOrderL, width);
		//DD.print2(data);
		timeDD = (System.currentTimeMillis()-timeDDstart)/1000;
	}

	/**
	 * Computes a heuristic ordering of leader variables for DD construction.
	 */
	private void orderingHeuristicL(List<Integer> varOrder, DataHandler data) {
		List<Item> order=new ArrayList<>();
		if(data.isInterdiction == false) {
			// Compute sum of coefficients
			int[] sumVarCoef = new int[data.P.nVarL];
			for (int i = 0; i < data.P.nVarL; i++) {
				int sum = 0;
				for (int j = 0; j < data.P.nConsF; j++) {
					sum+=Math.abs(data.P.Bfol[j][i]);
				}
				sumVarCoef[i] = sum;
			}
			for (int i = 0; i < data.P.nVarL; i++) {
				order.add(new Item(i, sumVarCoef[i] ));
			}
			//Collections.sort(order);
			Collections.sort(order, Collections.reverseOrder());
			for (int i = 0; i < order.size(); i++) {
				varOrder.add(order.get(i).pos);
				System.out.println("VAR "+order.get(i).pos+" coeff: "+order.get(i).val);
			}
		}
		else {
			// Compute sum of coefficients
			int[] sumVarCoef = new int[data.P.nVarL];
			for (int i = 0; i < data.P.nVarL; i++) {
				int sum = 0;
				for (int j = 0; j < data.P.nConsL; j++) {
					sum+=Math.abs(data.P.Aleader[j][i]);
				}
				sumVarCoef[i] = sum;
			}
			for (int i = 0; i < data.P.nVarF; i++) {
				order.add(new Item(i, (int) data.P.cfol[i] ));
			}
			//Collections.sort(order);
			Collections.sort(order, Collections.reverseOrder());
			for (int i = 0; i < order.size(); i++) {
				varOrder.add(order.get(i).pos);
				// Debug: inspect the ordering score assigned to this variable.
			}

		}

	}

	public BDDGraph buildBDDRelaxed(DataHandler data, List<Integer> varOrderL, int width) throws IloException {
		int nodeCount = 0;
		BDDGraph dd = new BDDGraph();
		dd.varOrder = varOrderL;
		Queue<Node> currentLayerQueue = new ArrayDeque<>(); 		//Nodes being processed at some layer
		Queue<Node> nextLayerQueue = new ArrayDeque<>();    		//Queue of Nodes to be processed at layer+1
		ArrayList<Node> nextLayerNodes = new ArrayList<Node>();		//All nodes generated (before relaxing!)
		Hashtable<String, Integer> nextLayerGraphMap = new Hashtable<String, Integer>();	//Indexes for new nodes in next layer

		// Create initial node
		Node root = new Node(data);
		root.ID = nodeCount;
		root.init(data);			// Initialize the state: 0 for F and B constraints, and RHS for L constraints
		root.genKey(0);
		dd.GraphMap.put(root.key, root.ID);
		dd.nodes.add(root);
		currentLayerQueue.add(root); //Add it to the queue
		nodeCount++;

		// Create final node
		Node sink = new Node(data);
		sink.key = "END";
		sink.ID = nodeCount;
		dd.GraphMap.put(sink.key, sink.ID);
		dd.nodes.add(sink);
		dd.lastNode = nodeCount;
		nodeCount++;

		// Compute completion bounds to prune
		int[] completionBoundsL = new int [data.P.nConsL];
		for (int q = 0; q < data.P.nConsL; q++) {
			for (int i = 0; i < data.P.nVarL; i++) {
				completionBoundsL[q] += Math.min(0, data.P.Aleader[q][i]);
			}
			for (int i = 0; i < data.P.nVarF; i++) {
				completionBoundsL[q] += Math.min(0, data.P.Bleader[q][i]);
			}
		}
		System.out.println("COMPLETION BOUNDS L: "+Arrays.toString(completionBoundsL));
		int[] completionBoundsF = new int [data.P.nConsF];
		for (int q = 0; q < data.P.nConsF; q++) {
			for (int i = 0; i < data.P.nVarL; i++) {
				completionBoundsF[q] += Math.min(0, data.P.Bfol[q][i]);
			}
		}
		System.out.println("COMPLETION BOUNDS F: "+Arrays.toString(completionBoundsF));

		///////////Process nodes in the queue for every position according to the order(layer by layer!) for the LEADER VARIABLES ***********************************************************************************
		for (int layer = 0; layer < data.P.nVarL; layer++) {
			int varIndex = varOrderL.get(layer);
			int layerNodeCount = 0;
			System.out.println("LEADER LAYER "+layer+" VAR "+varIndex+" QUEUE "+currentLayerQueue.size());
			// Update completion bounds
			for (int q = 0; q < data.P.nConsL; q++) {
				completionBoundsL[q] -= Math.min(0, data.P.Aleader[q][varIndex]);
			}
			for (int q = 0; q < data.P.nConsF; q++) {
				completionBoundsF[q] -= Math.min(0, data.P.Bfol[q][varIndex]);
			}
			//System.out.println("COMPLETION BOUNDS: "+Arrays.toString(completionBounds));
			while(!currentLayerQueue.isEmpty())
			{
				Node node=currentLayerQueue.poll();

				// Create nodes from the parent
				Node noHead = node.createZeroNode();  
				Node yesHead = node.createOneNode(varIndex, data);

				////////////////////////////////////////// Add zero-arc
				//Add node to the potential nodes for next layer
				noHead.genKey(layer+1);
				// Try to prune first
				if(noHead.checkBounds(data, completionBoundsL, completionBoundsF)) {
					// Check if the node has been created before
					if(nextLayerGraphMap.containsKey(noHead.key)==false) {
						nextLayerGraphMap.put(noHead.key, layerNodeCount);
						nextLayerNodes.add(noHead);
						layerNodeCount++;

						// Create arc
						Arc noArc = new Arc(dd.GraphMap.get(node.key), -1, 0, varIndex, true ); //ID to the head not yet assigned!
						dd.arcs.add(noArc);           
						node.out.add(noArc);
						noHead.in.add(noArc);
					}
					else {
						// Update L
						nextLayerNodes.get(nextLayerGraphMap.get(noHead.key)).updateL(data, noHead.Lfol, noHead.Llead, noHead.minCost);
						// Create arc
						Arc noArc = new Arc(dd.GraphMap.get(node.key), -1, 0, varIndex, true );
						dd.arcs.add(noArc);           
						node.out.add(noArc);
						nextLayerNodes.get(nextLayerGraphMap.get(noHead.key)).in.add(noArc);
					}
				}
				//////////////////////////////////////////// Add one-arc 
				if(yesHead.checkBounds(data, completionBoundsL, completionBoundsF)) {
					//Add node to the DD and to the queue for next layer
					yesHead.genKey(layer+1);
					// Check if the node has been created before
					if(nextLayerGraphMap.containsKey(yesHead.key)==false) {
						nextLayerGraphMap.put(yesHead.key, layerNodeCount);
						nextLayerNodes.add(yesHead);
						layerNodeCount++;

						// Create arc
						Arc arc = new Arc(dd.GraphMap.get(node.key), -1, 1, varIndex, true);
						dd.arcs.add(arc);                         
						node.out.add(arc);
						yesHead.in.add(arc);
					}
					else {
						// Update L
						nextLayerNodes.get(nextLayerGraphMap.get(yesHead.key)).updateL(data, yesHead.Lfol, yesHead.Llead, yesHead.minCost);
						// Create arc
						Arc arc = new Arc(dd.GraphMap.get(node.key), -1, 1, varIndex, true);
						dd.arcs.add(arc);                         
						node.out.add(arc);
						nextLayerNodes.get(nextLayerGraphMap.get(yesHead.key)).in.add(arc);
					}
				}
			}
			if(layerNodeCount <= width) { // No need to relax!
				// Add nodes to the DD and to the queue
				for (int i = 0; i < nextLayerNodes.size(); i++) {
					Node n = nextLayerNodes.get(i);
					n.ID = nodeCount;
					// Update correct index for incoming arcs
					for (int j = 0; j < n.in.size(); j++) {
						n.in.get(j).head = n.ID;
					}
					// Add Node
					dd.GraphMap.put(n.key, n.ID);
					dd.nodes.add(n);
					nextLayerQueue.add(n);
					nodeCount++;
				}
				nextLayerGraphMap.clear();
				nextLayerNodes.clear();
			}
			else {//Time to relax!!!
				System.out.println("TIME TO RELAX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WE HAVE "+nextLayerNodes.size());
				isExact = false;
				//relaxCost(data, nextLayerNodes, width);
				//relaxRHS(data, nextLayerNodes, width);
				relaxValApprox(data, nextLayerNodes, width);
				//System.out.println("AFTER RELAX!!! WE HAVE "+nextLayerNodes.size());
				// Add nodes to the DD and to the queue
				for (int i = 0; i < nextLayerNodes.size(); i++) {
					Node n = nextLayerNodes.get(i);
					n.ID = nodeCount;
					// Update correct index for incoming arcs
					for (int j = 0; j < n.in.size(); j++) {
						n.in.get(j).head = n.ID;
					}
					// Add Node
					dd.GraphMap.put(n.key, n.ID);
					dd.nodes.add(n);
					nextLayerQueue.add(n);
					nodeCount++;
				}
				nextLayerGraphMap.clear();
				nextLayerNodes.clear();

			}

			currentLayerQueue=nextLayerQueue;
			nextLayerQueue=new ArrayDeque<>();
			System.gc();
		}

		// Add artificial zero arcs to the end node
		while( !currentLayerQueue.isEmpty() ){
			Node node=currentLayerQueue.poll();
			Arc noArc = new Arc(dd.GraphMap.get(node.key), dd.GraphMap.get("END"), 0, -1, false);
			dd.arcs.add(noArc);
			node.out.add(noArc);
		}
		//Fill out layers
		dd.layers = new ArrayList[data.P.nVarL+1];
		for (int i = 0; i < data.P.nVarL+1; i++) {
			dd.layers[i] = new ArrayList<Node>();
		}
		dd.layers[0].add(root);
		for (int i = 2; i < dd.nodes.size() ; i++) { // 0 and 1 are root node and end node
			dd.layers[dd.nodes.get(i).layer].add(dd.nodes.get(i));
		}

		//dd.print3(data);
		if(!isExact) {
			cplexAP.clearModel();
			cplexAP.end();
			cplexAP = null; 
			System.gc();
		}

		return dd;
	}



	/**
	 * Relaxes one DD layer by merging nodes with similar shortest-path costs.
	 */
	private void relaxCost(DataHandler data, ArrayList<Node> layer, int width) {
		ArrayList<Node> relaxedLayer = new ArrayList<Node>();		//Relaxed layer
		Hashtable<String, Integer> GraphMap = new Hashtable<String, Integer>();	//Indexes for new relaxed nodes in relaxed layer
		//for (int i = 0; i < layer.size(); i++) {
		// Debug: inspect the current layer before merging.
		//}
		//Order and merge!
		List<Item> order = new ArrayList<>();
		// Compute sum of coefficients
		for (int i = 0; i < layer.size(); i++) {
			order.add(new Item(i, layer.get(i).minCost ));
		}
		Collections.sort(order);
		int numMerge = 0;
		for (int i = 0; i < order.size()-1; i+=2) {
			System.out.println("About to MERGE "+layer.get(order.get(i).pos).minCost);			
			if(layer.size()-numMerge > width) {// Need to merge!
				Node relaxed = layer.get(order.get(i).pos).relax(layer.get(order.get(i+1).pos));
				layer.get(order.get(i).pos).removed = true;
				layer.get(order.get(i+1).pos).removed = true;

				//System.out.println("MERGING "+layer.get(order.get(i).pos).minCost+" and "+layer.get(order.get(i+1).pos).minCost+" RELAXED: "+relaxed.minCost);
				if(GraphMap.containsKey(relaxed.key) == false) {
					GraphMap.put(relaxed.key, relaxedLayer.size());
					relaxedLayer.add(relaxed);
				}
				else {//Add the incoming arcs to the existing node
					// Debug: repeated relaxed node detected.
					relaxedLayer.get(GraphMap.get(relaxed.key)).in.addAll(relaxed.in);
					//Update L and minCost
					relaxedLayer.get(GraphMap.get(relaxed.key)).updateL(data, relaxed.Lfol, relaxed.Llead, relaxed.minCost);
				}
				numMerge++;
			}
			if(layer.size()-numMerge <= width) {// Stop merging!!!!
				i = order.size();
			}

		}
		// Add remaining not removed nodes and not repeated
		for (int i = 0; i < layer.size(); i++) {
			if(layer.get(i).removed == false) {
				if(GraphMap.containsKey(layer.get(i).key) == false) {
					relaxedLayer.add(layer.get(i));
				}
				else {
					//Reroute the incoming arcs to the existing node
					relaxedLayer.get(GraphMap.get(layer.get(i).key)).in.addAll(layer.get(i).in);
					//Update L and minCost
					relaxedLayer.get(GraphMap.get(layer.get(i).key)).updateL(data, layer.get(i).Lfol, layer.get(i).Llead, layer.get(i).minCost);
				}
			}
			//else {
			//	System.out.println("NODE REMOVED: "+layer.get(i).key);
			//}
		}

		//for (int i = 0; i < relaxedLayer.size(); i++) {
		//	Debug: inspect the relaxed layer after merging.
		//}

		layer.clear();
		layer.addAll(relaxedLayer);
	}

	/**
	 * Relaxes one DD layer by merging nodes with similar accumulated RHS
	 * information.
	 */
	private void relaxRHS(DataHandler data, ArrayList<Node> layer, int width) {
		ArrayList<Node> relaxedLayer = new ArrayList<Node>();		//Relaxed layer
		Hashtable<String, Integer> GraphMap = new Hashtable<String, Integer>();	//Indexes for new relaxed nodes in relaxed layer
		//for (int i = 0; i < layer.size(); i++) {
		// Debug: inspect the current layer before merging.
		//}
		List<Item> order = new ArrayList<>();
		// Compute sum of coefficients
		int[] sumVarCoef = new int[layer.size()];
		for (int i = 0; i < layer.size(); i++) {
			int sum = 0;
			/*for (int j = 0; j < layer.get(i).S.length; j++) {
				if(data.cInteraction[j].equals("B")) { //Only care about "B" constraints
					sum+=Math.abs(layer.get(i).S[j]);
				}
			}*/
			for (int j = 0; j < layer.get(i).Llead.length; j++) {
				sum+=Math.abs(layer.get(i).Llead[j]);
			}
			sumVarCoef[i] = sum;
		}
		for (int i = 0; i < layer.size(); i++) {
			order.add(new Item(i, sumVarCoef[i] ));
		}
		//Collections.sort(order, Collections.reverseOrder());
		Collections.sort(order);
		int numMerge = 0;
		for (int i = 0; i < order.size()-1; i+=2) {
			//System.out.println("About to MERGE "+order.get(i).val);			
			if(layer.size()-numMerge > width) {// Need to merge!
				Node relaxed = layer.get(order.get(i).pos).relax(layer.get(order.get(i+1).pos));
				layer.get(order.get(i).pos).removed = true;
				layer.get(order.get(i+1).pos).removed = true;

				//System.out.println("MERGING "+layer.get(order.get(i).pos).getS(data)+" and "+layer.get(order.get(i+1).pos).getS(data)+" RELAXED: "+relaxed.getS(data));
				if(GraphMap.containsKey(relaxed.key) == false) {
					GraphMap.put(relaxed.key, relaxedLayer.size());
					relaxedLayer.add(relaxed);
				}
				else {//Add the incoming arcs to the existing node
					// Debug: repeated relaxed node detected.
					relaxedLayer.get(GraphMap.get(relaxed.key)).in.addAll(relaxed.in);
					//Update L and minCost
					relaxedLayer.get(GraphMap.get(relaxed.key)).updateL(data, relaxed.Lfol, relaxed.Llead, relaxed.minCost);
				}
				numMerge++;
			}
			if(layer.size()-numMerge <= width) {// Stop merging!!!!
				i = order.size();
			}

		}
		// Add remaining not removed nodes and not repeated
		for (int i = 0; i < layer.size(); i++) {
			if(layer.get(i).removed == false) {
				if(GraphMap.containsKey(layer.get(i).key) == false) {
					relaxedLayer.add(layer.get(i));
				}
				else {
					//Reroute the incoming arcs to the existing node
					relaxedLayer.get(GraphMap.get(layer.get(i).key)).in.addAll(layer.get(i).in);
					//Update L and minCost
					relaxedLayer.get(GraphMap.get(layer.get(i).key)).updateL(data, layer.get(i).Lfol, layer.get(i).Llead, layer.get(i).minCost);
				}
			}
			//else {
			//	System.out.println("NODE REMOVED: "+layer.get(i).key);
			//}
		}

		//for (int i = 0; i < relaxedLayer.size(); i++) {
		//	Debug: inspect the relaxed layer after merging.
		//}

		layer.clear();
		layer.addAll(relaxedLayer);
	}

	// Sort according to an approximation of the value function
	private void relaxValApprox(DataHandler data, ArrayList<Node> layer, int width) throws IloException {
		ArrayList<Node> relaxedLayer = new ArrayList<Node>();		//Relaxed layer
		Hashtable<String, Integer> GraphMap = new Hashtable<String, Integer>();	//Indexes for new relaxed nodes in relaxed layer
		int[] approx = new int[layer.size()]; 
		int min = 999999999;
		int max = -999999999;
		
		// andre change here
//		int numBuckets = width;
		int numBuckets = 5;

		if(cplexAP == null) {
			// Ini model
			cplexAP = new IloCplex();
			cplexAP.setParam(IloCplex.Param.Threads, numCores);
			cplexAP.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0.005);
			cplexAP.setOut(null);

			//Variables 
			yAP = new IloNumVar[data.P.nVarF];
			for (int i = 0; i < data.P.nVarF; i++) {
				yAP[i] = cplexAP.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
			}
			// add constraints follower
			AforAP = new IloRange[data.P.nConsF];
			for(int j=0; j < data.P.nConsF; j++){
				IloLinearNumExpr exprf = cplexAP.linearNumExpr();
				for (int i = 0; i <data.P.nVarF; i++) {
					exprf.addTerm(data.P.Afol[j][i], yAP[i]);
				}
				AforAP[j] = cplexAP.addLe(exprf, data.P.bfol[j]); //Reduce the RHS according to node state
			}

			// ADD OBJECTIVE FUNCTION	
			IloLinearNumExpr obj = cplexAP.linearNumExpr();
			for(int i=0; i<data.P.nVarF; i++){
				obj.addTerm(yAP[i], data.P.cfol[i]);
			}
			cplexAP.addMinimize(obj,"FollObj");
		}

		for (int i = 0; i < layer.size(); i++) {
			Sample s = solveFollowerApprox(data, layer.get(i).S, AforAP); //This one is not solve to optimality
			approx[i] = s.folObj; 
			//System.out.println("Val function done "+i);
			if(approx[i] < min) {min = approx[i];}
			if(approx[i] > max) {max = approx[i];}
			// Add to the Sample
			//s.print();
			if(s.folObj < folObjM && sampleMap.containsKey(s.key) == false) {
				sampleMaster.add(s);
				sampleMap.put(s.key, sampleMaster.size());
				// Debug: sample added to the pool.
			}

		}

		// ==========================================================================================================
		// ALT. A
		
		// // Compute break points for buckets
		// double[] breakpoints = new double[numBuckets+1];
		// double bucketSize = (max-min)/(numBuckets+0.0);
		// breakpoints[0] = min;
		// for (int k = 0; k < numBuckets; k++) {
		// 	//System.out.println("BreakPoint "+k+": "+breakpoints[k]);
		// 	breakpoints[k+1] = breakpoints[k]+bucketSize;
		// }
		// // Distribute the nodes to each bucket
		// ArrayList<Node>[] buckets = new ArrayList[width];
		// for (int k = 0; k < numBuckets; k++) {
		// 	buckets[k] = new ArrayList<Node>();
		// }
		// for (int i = 0; i < layer.size(); i++) {
		// 	for (int k = 0; k < numBuckets; k++) {
		// 		if(approx[i] >= breakpoints[k]-0.01 && approx[i] <= breakpoints[k+1]+0.01) {
		// 			System.out.println("Node "+i+" val "+approx[i]+" goes to bucket "+k);
		// 			buckets[k].add(layer.get(i));
		// 			k = numBuckets;
		// 		}
		// 	}
		// }

		
		// ==========================================================================================================
		// ALT. B	

		// Instead of computing breakpoints by min/max, do equal-frequency bucketing
		int n = layer.size();
		int B = numBuckets;

		// Sort indices by approx value
		Integer[] idx = new Integer[n];
		for (int i = 0; i < n; i++) idx[i] = i;
		Arrays.sort(idx, (i, j) -> Double.compare(approx[i], approx[j]));

		// Prepare buckets
		ArrayList<Node>[] buckets = new ArrayList[B];
		for (int b = 0; b < B; b++) {
			buckets[b] = new ArrayList<>();
		}

		// Distribute approximately n/B nodes per bucket
		int base = n / B, rem = n % B; // some buckets may have one extra
		int pos = 0;
		for (int b = 0; b < B; b++) {
			int size = base + (b < rem ? 1 : 0);
			for (int c = 0; c < size; c++, pos++) {
				int nodeIdx = idx[pos];
				buckets[b].add(layer.get(nodeIdx));
				System.out.println("Node " + nodeIdx + " val " + approx[nodeIdx] + " goes to bucket " + b);
			}
		}

		// ==========================================================================================================
		// ALT. C
		// Inputs assumed available:
		// List<Node> layer;              // nodes
		// double[] approx;               // values for nodes (same order as 'layer')
		// int numBuckets;                // B

//		int n = layer.size();
//		int B = Math.max(1, numBuckets);
//
//		// ----- Pre-allocate outputs so they are always available elsewhere -----
//		@SuppressWarnings("unchecked")
//		ArrayList<Node>[] buckets = new ArrayList[B];
//		for (int b = 0; b < B; b++) buckets[b] = new ArrayList<>();
//		double[] b_min = new double[B];
//		double[] b_max = new double[B];
//
//		// Collect numeric values with indices and sort
//		class Pair { int idx; double v; Pair(int i, double vv){ idx=i; v=vv; } }
//		ArrayList<Pair> pts = new ArrayList<>(n);
//		for (int i = 0; i < n; i++) if (!Double.isNaN(approx[i])) pts.add(new Pair(i, approx[i]));
//		pts.sort(Comparator.comparingDouble(p -> p.v));
//
//		if (pts.isEmpty()) {
//			// No numeric data: keep buckets empty; mark intervals as NaN, epsilon = 0
//			Arrays.fill(b_min, Double.NaN);
//			Arrays.fill(b_max, Double.NaN);
//			double epsilon = 0.0;
//			System.out.println("No numeric values. B=" + B + " epsilon=" + epsilon);
//		} else {
//			double globalMin = pts.get(0).v;
//			double globalMax = pts.get(pts.size()-1).v;
//
//			// Fast path: if B >= number of points, epsilon = 0 works (singletons)
//			double epsilon;
//			if (B >= pts.size()) {
//				epsilon = 0.0;
//				int k = 0;
//				for (; k < pts.size(); k++) {
//					Pair p = pts.get(k);
//					buckets[k].add(layer.get(p.idx));
//					b_min[k] = p.v;
//					b_max[k] = p.v;   // width 0
//				}
//				// Pad remaining empty buckets with zero-width intervals at the last value
//				double anchor = pts.get(pts.size()-1).v;
//				for (; k < B; k++) { b_min[k] = anchor; b_max[k] = anchor; }
//			} else {
//				// ----- Binary search minimal epsilon so that <= B intervals cover all points -----
//				java.util.function.DoublePredicate feasible = (double eps) -> {
//					int used = 0, i = 0, m = pts.size();
//					while (i < m) {
//						used++;
//						double end = pts.get(i).v + eps;
//						int j = i;
//						while (j + 1 < m && pts.get(j + 1).v <= end) j++;
//						i = j + 1;
//						if (used > B) return false;
//					}
//					return true;
//				};
//
//				double lo = 0.0, hi = globalMax - globalMin;
//				for (int it = 0; it < 60; it++) {
//					double mid = (lo + hi) / 2.0;
//					if (feasible.test(mid)) hi = mid; else lo = mid;
//				}
//				epsilon = hi;
//
//				// ----- Reconstruct buckets with width = epsilon -----
//				int i = 0, m = pts.size(), k = 0;
//				while (i < m && k < B) {
//					double start = pts.get(i).v;
//					double end = start + epsilon;
//					b_min[k] = start;
//					b_max[k] = end;
//
//					int j = i;
//					while (j < m && pts.get(j).v <= end) {
//						buckets[k].add(layer.get(pts.get(j).idx));
//						j++;
//					}
//					i = j;
//					k++;
//				}
//
//				// If we used fewer than B buckets, pad remaining with empty intervals of same width
//				// placed contiguously after the last end (or at globalMin if none used).
//				double anchorStart = (k > 0 ? b_max[k - 1] : globalMin);
//				for (; k < B; k++) {
//					b_min[k] = anchorStart;
//					b_max[k] = anchorStart + epsilon;
//					// keep buckets[k] as the pre-allocated empty list
//					anchorStart = b_max[k]; // place next empty interval after this one
//				}
//			}
//
//			 // ---- Optional logging ----
//			 System.out.println("globalMin=" + globalMin + " globalMax=" + globalMax + " epsilon=" + epsilon);
//			 for (int b = 0; b < B; b++) {
//			     System.out.println("Bucket " + b + " [" + b_min[b] + ", " + b_max[b] + "] size=" + buckets[b].size());
//			 }
//		}
							
		// Values
		System.out.println("Values: ");
		for (int i = 0; i < layer.size(); i++) {
			System.out.print(approx[i] + " ");
		}
		System.out.println("");

		// Buckets
		for (int k = 0; k < numBuckets; k++) {
			System.out.println("Bucket: " + k);
			// System.out.println("\tBucket breakpoint: " + breakpoints[k]);			
			System.out.println("\tBucket size: " + buckets[k].size());
		}
		
		// Now merge!!!!!
		int numMerge = 0;

//		for (int k = numBuckets-1; k >= 0; --k) {
		for (int k = 0; k < numBuckets; k++) {
			if(buckets[k].size() >= 2) {//Merge all into one 
				System.out.println("Processing bucket " + k);

				Node relaxed = buckets[k].get(0).relax(buckets[k].get(1));
				buckets[k].get(0).removed = true;
				buckets[k].get(1).removed = true;
				//System.out.println("MERGING "+buckets[k].get(0).getS(data)+" and "+buckets[k].get(1).getS(data)+" RELAXED: "+relaxed.getS(data));
				numMerge++;
				for (int q = 2; q < buckets[k].size(); q++) {
					relaxed = buckets[k].get(q).relax(relaxed);
					buckets[k].get(q).removed = true;
					//System.out.println("---------MERGING WITH"+buckets[k].get(q).getS(data)+" RELAXED: "+relaxed.getS(data));
					numMerge++;
				}

				if(GraphMap.containsKey(relaxed.key) == false) {
					GraphMap.put(relaxed.key, relaxedLayer.size());
					relaxedLayer.add(relaxed);
				}
				else {//Add the incoming arcs to the existing node
					// Debug: repeated relaxed node detected.
					relaxedLayer.get(GraphMap.get(relaxed.key)).in.addAll(relaxed.in);
					//Update L and minCost
					relaxedLayer.get(GraphMap.get(relaxed.key)).updateL(data, relaxed.Lfol, relaxed.Llead, relaxed.minCost);
				}
			}

			System.out.println("\tAfter processing: " + (layer.size()-numMerge));
			if(layer.size()-numMerge <= width) {// Stop merging!!!!
				break;
				//k = numBuckets;
				//k = -1;
			}
		}
		System.out.println("Final layer size after merging: " + (layer.size()-numMerge));

		// Add remaining not removed nodes and not repeated
		for (int i = 0; i < layer.size(); i++) {
			if(layer.get(i).removed == false) {
				if(GraphMap.containsKey(layer.get(i).key) == false) {
					relaxedLayer.add(layer.get(i));
				}
				else {
					//Reroute the incoming arcs to the existing node
					relaxedLayer.get(GraphMap.get(layer.get(i).key)).in.addAll(layer.get(i).in);
					//Update L and minCost
					relaxedLayer.get(GraphMap.get(layer.get(i).key)).updateL(data, layer.get(i).Lfol, layer.get(i).Llead, layer.get(i).minCost);
				}
			}
			//else {
			//	System.out.println("NODE REMOVED: "+layer.get(i).key);
			//}
		}

		//for (int i = 0; i < relaxedLayer.size(); i++) {
		//	Debug: inspect the relaxed layer after merging.
		//}

		layer.clear();
		layer.addAll(relaxedLayer);

	}

	/**
	 * Rounds a value to four decimal places for consistent numerical reporting.
	 */
	public double round(double value) {
		double rounded;

		rounded = Math.round(value*10000)/10000.0;

		return rounded;
	}

	/**
	 * Solves the discrete non-interdiction reformulation with blocking
	 * constraints, without callbacks.
	 */
	public void singleLevelDDRelaxedNC(DataHandler data, int timeLim, double Atime) throws IloException {
		// Solve the follower problem for each node/arc in the last layer (not solved before)
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				if(DD.nodes.get(nodeIndex).cutGenerated == false) {
					//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
					int folObjValR = solveFollower(data, DD.nodes.get(nodeIndex).S); 
					DD.arcs.get(i).cost = folObjValR;
					DD.nodes.get(nodeIndex).cutGenerated = true;
					DD.nodes.get(nodeIndex).folResponse = folObjValR;
					if (folObjValR > folUB) {folUB = folObjValR;}
					//System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
				}
				else {
					DD.arcs.get(i).cost = Math.min(DD.nodes.get(nodeIndex).folResponse, folObjM);
				}
			}
		}
		// Propagate costs backwards
		for (int i = DD.nodes.size()-1; i >= 2 ; i--) { // 0 and 1 are root node and end node
			DD.nodes.get(i).propagateCost(data);
		}
		System.out.println(" FOLL OBJ UB: "+folUB);

		int timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		double lim = Math.max( 10, timeLim-timeUsed); //For safety!
		cplex.setParam(IloCplex.Param.TimeLimit, lim);
		if(isExact==false) {cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0.1);}		//Don't need to close the gap at the beginning, then switch to 0 gap

		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}

		// Flow variables for the BDDs
		IloNumVar[] w = new IloNumVar[DD.arcs.size()]; 
		for (int i = 0; i < DD.arcs.size(); i++) {
//			w[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "w");
			w[i] = cplex.numVar(0, 1, IloNumVarType.Float, "w");
		}
		// Blocking variables
		IloNumVar[][] v = new IloNumVar[data.P.nConsF][sampleMaster.size()];
		for (int i = 0; i < data.P.nConsF; i++) {
			for(int j = 0; j < sampleMaster.size(); j++) {
				v[i][j] = cplex.numVar(0, 1, IloNumVarType.Bool, "v");
			}
		}
		// block variables for the new VIS
		IloNumVar[][] r = new IloNumVar[data.P.nConsF][data.P.nVarF]; 
		for (int i = 0; i < data.P.nVarF; i++) {
			for(int j = 0; j < data.P.nConsF; j++)
				if(data.P.Afol[j][i] <= 0 ) {// Impossible to block!
					r[j][i] = cplex.numVar(0, 0, IloNumVarType.Float, "r");
				}
				else {//Here we can block
					r[j][i] = cplex.numVar(0, 1, IloNumVarType.Bool, "r");
				}
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}
		if(logical) {
			// New VIs: each follower binary variable with negative obj coefficient set to zero only if blocked 
			int[] varRecord = new int[data.P.nVarF];  //Keep track of which variables are already represented in the Sample
			for (int q = 0; q < sampleMaster.size(); q++) {
				for(int i=0; i < data.P.nVarF; i++){
					if(sampleMaster.get(q).y[i] > 0.9) {
						varRecord[i]++;
					}
				}
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				if(data.P.cfol[i] < 0  ) {
					IloLinearNumExpr exprf = cplex.linearNumExpr();
					exprf.addTerm(1, y[i]);
					for(int j=0; j < data.P.nConsF; j++){
						exprf.addTerm(1, r[j][i]);	//Block it by some constraint
					}
					cplex.addGe(exprf, 1);
					// Now define the blocking variables
					for(int j=0; j < data.P.nConsF; j++){
						int block = (int) (data.P.bfol[j]-data.P.Afol[j][i]+1);
						IloLinearNumExpr exprB = cplex.linearNumExpr();
						for (int q = 0; q <data.P.nVarF; q++) {
							if(q!=i) {exprB.addTerm(data.P.Afol[j][q], y[q]);}
						}
						for (int q = 0; q <data.P.nVarL; q++) {
							exprB.addTerm(data.P.Bfol[j][q], x[q]);
						}
						exprB.addTerm(-block-999999999, r[j][i]);
						//System.out.println("CONSTRAINT "+j+" block: "+block);
						cplex.addGe(exprB, -999999999);	
					}
				}
			}
		}
		// FLOW CONSTRAINTS	
		// Add flow constraints (but drop the last one)
		IloLinearNumExpr[] expr = new IloLinearNumExpr[DD.nodes.size()];
		for (int n = 0; n < DD.nodes.size(); n++) {
			expr[n] = cplex.linearNumExpr();
		}	
		for (int i = 0; i < DD.arcs.size(); i++) {
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			expr[tail].addTerm(1, w[i]);
			expr[head].addTerm(-1, w[i]);
			//System.out.println("FLOW FROM "+DD.nodes.get(tail).key+" to "+DD.nodes.get(head).key+"        "+tail+" to "+head+" X"+DD.arcs.get(i).varIndex+" val: "+DD.arcs.get(i).arcVal);
		}
		cplex.addEq(expr[0], 1, "c0");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){cplex.addEq(expr[n], 0, "flow_"+0+"_"+n);}
		}

		// Link leader variables
		IloLinearNumExpr[] exprLink = new IloLinearNumExpr[data.P.nVarL];
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i] = cplex.linearNumExpr();
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {	
				exprLink[DD.arcs.get(i).varIndex].addTerm(1, w[i]);
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i].addTerm(-1, x[i]);
			cplex.addEq(exprLink[i], 0);
		}

		// Select maximum one "blocker" for each constraint
		for (int j = 0; j < data.P.nConsF; j++) {
			IloLinearNumExpr exprZ = cplex.linearNumExpr();
			for (int q = 0; q < sampleMaster.size(); q++) {
				exprZ.addTerm(1, v[j][q]);
			}
			cplex.addLe(exprZ, 1);
		}

		//Add constraints for each element of the sample
		for (int q = 0; q < sampleMaster.size(); q++) {
			// Add constraint on theta
			IloLinearNumExpr exprV = cplex.linearNumExpr();
			for(int i=0; i < data.P.nVarF; i++){
				exprV.addTerm(y[i], data.P.cfol[i]);
			}
			for (int i = 0; i < sampleMaster.get(q).block.size(); i++) {
				exprV.addTerm(v[sampleMaster.get(q).block.get(i)[0]][q], -(folObjM-sampleMaster.get(q).folObj)); 
				for (int k = 0; k < sampleMaster.size(); k++) {
					if(k!=q) {
						for (int j = 0; j < sampleMaster.get(k).block.size(); j++) {
							if(sampleMaster.get(q).block.get(i)[0] == sampleMaster.get(k).block.get(j)[0] && sampleMaster.get(k).block.get(j)[1] >= sampleMaster.get(q).block.get(i)[1]) {
								exprV.addTerm(v[sampleMaster.get(k).block.get(j)[0]][k], -(folObjM-sampleMaster.get(q).folObj));
							}
						}
					}
				}
			}			
			cplex.addLe(exprV, sampleMaster.get(q).folObj);
			// Add block definitions
			for (int i = 0; i < sampleMaster.get(q).block.size(); i++) { 
				IloLinearNumExpr exprB = cplex.linearNumExpr();
				int index = sampleMaster.get(q).block.get(i)[0];
				exprB.addTerm(v[index][q], MfolConstraint[index]-sampleMaster.get(q).block.get(i)[1]);
				for (int k = 0; k <data.P.nVarL; k++) {
					exprB.addTerm(data.P.Bfol[sampleMaster.get(q).block.get(i)[0]][k], x[k]);
				}
				cplex.addGe(exprB, MfolConstraint[index]);
			}

		}

		///////////////////////////////////////////////////// VALID INEQUALITIES ///////////////////////////////////////////
		IloLinearNumExpr exprB = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprB.addTerm(y[q], data.P.cfol[q]);
		}
		cplex.addLe(exprB, folUB);

		/////////////////////////////////////////////////// DD VALUE FUNCTION CONSTRAINT/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			exprV.addTerm(w[i], -DD.arcs.get(i).cost);
		}
		cplex.addLe(exprV, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}
		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		boolean stop = false;
		int iter = 0;
		int blockRecord[] = new int[data.P.nConsF];  //Keep track of which constraints we are using to block
		while(UB-LB >0 && !stop) {
			cplex.solve();			
			if(iter >= 250 || timeLim-timeUsed <= 900  || Math.abs((UB-LB)/(UB+0.0)) <= 0.05 || (cplex.getStatus() == IloCplex.Status.Feasible || cplex.getStatus() == IloCplex.Status.Optimal)  && round(cplex.getObjValue()) == UB) { //Back to default to solve to optimality!!
				cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-6);
				System.out.println("READY TO GO OPTIMAL IN THE MASTER PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			}
			int numSolutions = cplex.getSolnPoolNsolns();
			if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
				System.out.println("		DD-IP-R is infeasible!");
				LB = UB; //If UB= 99999 then problem is infeasible, otherwise cplex just closed the gap!
				stop = true;
			}
			else if(cplex.getStatus() == IloCplex.Status.Optimal) {
				LB = (int) Math.max(Math.round(cplex.getBestObjValue()), LB);
				System.out.println("******************************************************  MIP ITERATION "+iter+" LB: "+ LB+" UB: "+UB);
				
				// Recover solution
				int[][] xVal = new int[numSolutions][data.P.nVarL];
				int[] yVal = new int[data.P.nVarF];
				for (int k = 0; k < numSolutions; k++) {
					for (int i = 0; i < data.P.nVarL; i++) {				
						if(cplex.getValue(x[i],k) > 0.9) {xVal[k][i] = 1;}
					}
				}
				for (int i = 0; i <data.P.nVarF; i++) {
					yVal[i] = (int)Math.round(cplex.getValue(y[i]));
				}
				// Record block use
				for (int i = 0; i < data.P.nConsF; i++) {
					for(int j = 0; j < sampleMaster.size(); j++) {
						if(cplex.getValue(v[i][j]) > 0.9) {blockRecord[i]++;} 
					}
				}
				System.out.println("x:"+Arrays.toString(xVal[0]));
				System.out.println("y:"+Arrays.toString(yVal));
				System.out.println("Block Record: "+Arrays.toString(blockRecord));
				// Try to find feasible completion for the leader variables
				Sample s = solveFollowerSample(data, getS(data, xVal[0]));
				if(getFolObj(data, yVal) <= s.folObj) {//Current solution is bilevel feasible!!!!
					if(cplex.getObjValue() < UB) {
						UB = (int)Math.round(cplex.getObjValue());
						System.out.println("Current solution is bilevel feasible!!! "+ LB+" UB: "+UB+" fol Obj: "+s.folObj);
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[0][i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = yVal[i];
						}
					}
				}
				else {
					int bigM = folUB-s.folObj;
					//System.out.println("BIG M is "+bigM);
					//Add follower optimality cut
					IloLinearNumExpr exprZ = cplex.linearNumExpr();
					for(int i=0; i<data.P.nVarF; i++){
						exprZ.addTerm(y[i], data.P.cfol[i]);
					}
					int oneCount = 0;
					for(int i=0; i<data.P.nVarL; i++){
						if(xVal[0][i] == 1) {
							exprZ.addTerm(bigM, x[i]);
							oneCount++;
						}
						else {
							exprZ.addTerm(-bigM, x[i]);
						}
					}

					cplex.addLe(exprZ, s.folObj+oneCount*bigM);
					numCuts++;

					//Try to find feasible solution
					int leaderObjVal = solveFixedLeader(data, xVal[0], s);
					if(leaderObjVal < Integer.MAX_VALUE && leaderObjVal < UB) {
						System.out.println("UB UPDATED TO "+leaderObjVal);
						UB = leaderObjVal;
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[0][i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = s.yopt[i];
						}
					}
				}

				// Add blocking constraints
				IloNumVar[] vp = new IloNumVar[s.block.size()];
				for (int i = 0; i < s.block.size(); i++) {
					vp[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
				}
				IloLinearNumExpr exprp = cplex.linearNumExpr();
				for(int i=0; i < data.P.nVarF; i++){
					exprp.addTerm(y[i], data.P.cfol[i]);
				}
				for (int i = 0; i < s.block.size(); i++) {
					exprp.addTerm(vp[i], -(folObjM-s.folObj)); //
				}			
				cplex.addLe(exprp, s.folObj);
				// Add block definitions
				for (int i = 0; i < s.block.size(); i++) { 
					IloLinearNumExpr exprQ = cplex.linearNumExpr();
					exprQ.addTerm(vp[i], MfolConstraint[s.block.get(i)[0]]-s.block.get(i)[1]);
					for (int k = 0; k <data.P.nVarL; k++) {
						exprQ.addTerm(data.P.Bfol[s.block.get(i)[0]][k], x[k]);
					}
					cplex.addGe(exprQ, MfolConstraint[s.block.get(i)[0]]);
				}
				// Add blocking constraints for other solutions in the pool
				for (int q = 1; q < numSolutions; q++) {
					// Solve follower
					Sample sprime = solveFollowerSample(data, getS(data, xVal[q]));
					IloNumVar[] vpQ = new IloNumVar[sprime.block.size()];
					for (int i = 0; i < sprime.block.size(); i++) {
						vpQ[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
					}
					IloLinearNumExpr exprpQ = cplex.linearNumExpr();
					for(int i=0; i < data.P.nVarF; i++){
						exprpQ.addTerm(y[i], data.P.cfol[i]);
					}
					for (int i = 0; i < sprime.block.size(); i++) {
						exprpQ.addTerm(vpQ[i], -(folObjM-sprime.folObj)); //
					}			
					cplex.addLe(exprpQ, sprime.folObj);
					// Add block definitions
					for (int i = 0; i < sprime.block.size(); i++) { 
						IloLinearNumExpr exprQ = cplex.linearNumExpr();
						exprQ.addTerm(vpQ[i], MfolConstraint[sprime.block.get(i)[0]]-sprime.block.get(i)[1]);
						for (int k = 0; k <data.P.nVarL; k++) {
							exprQ.addTerm(data.P.Bfol[sprime.block.get(i)[0]][k], x[k]);
						}
						cplex.addGe(exprQ, MfolConstraint[sprime.block.get(i)[0]]);
					}					
					System.out.println("Solution "+q+" val "+sprime.folObj);
					numCuts++;

					//Try to find feasible solution
					int leaderObjVal = solveFixedLeader(data, xVal[q], sprime);
					if(leaderObjVal < Integer.MAX_VALUE && leaderObjVal < UB) {
						System.out.println("UB UPDATED TO "+leaderObjVal);
						UB = leaderObjVal;
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[q][i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = sprime.yopt[i];
						}
					}

				}
			}
			else {
				System.out.println("ABORTED BECAUSE OF TIME LIMIT ");
				LB = (int) Math.max(Math.round(cplex.getBestObjValue()), LB); 
				if (cplex.getStatus() == IloCplex.Status.Feasible) {
					// Recover solution
					int[] xVal = new int[data.P.nVarL];
					int[] yVal = new int[data.P.nVarF];
					for (int i = 0; i < data.P.nVarL; i++) {
						xVal[i] = (int)Math.round(cplex.getValue(x[i]));
					}
					for (int i = 0; i <data.P.nVarF; i++) {
						yVal[i] = (int)Math.round(cplex.getValue(y[i]));
					}

					// Try to find feasible completion for the leader variables
					Sample s = solveFollowerSample(data, getS(data, xVal));
					if(getFolObj(data, yVal) <= s.folObj && cplex.getObjValue() < UB) {
						UB = (int) Math.round(cplex.getObjValue());
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = yVal[i];
						}
					}
					else {
						//Try to find feasible solution
						int leaderObjVal = solveFixedLeader(data, xVal, s);
						if(leaderObjVal < Integer.MAX_VALUE && leaderObjVal < UB) {
							System.out.println("UB UPDATED TO "+leaderObjVal);
							UB = leaderObjVal;
							// Record solution
							xstar = new int[data.P.nVarL];
							ystar = new int[data.P.nVarF];
							for (int i = 0; i < data.P.nVarL; i++) {
								xstar[i] = xVal[i];
							}
							for (int i = 0; i < data.P.nVarF; i++) {
								ystar[i] = s.yopt[i];
							}
						}
					}
				}
				stop = true;
			}
			if(stop == false) {
				timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
				if(timeUsed >= timeLim-10) {stop = true;}
				else {
					cplex.setParam(IloCplex.Param.TimeLimit, timeLim-timeUsed);
					System.out.println("TIME USED: "+timeUsed);
				}
			}
			iter++;
		}


	}


	// Enhanced with the blocking constraints for interdiction, no callbacks
	public void singleLevelDDInterdiction(DataHandler data, int timeLim, double Atime) throws IloException {

		// Solve the follower problem for each node/arc in the last layer (not solved before)
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				if(DD.nodes.get(nodeIndex).cutGenerated == false) {
					//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
					int folObjValR = solveFollower(data, DD.nodes.get(nodeIndex).S); 
					DD.arcs.get(i).cost = folObjValR;
					DD.nodes.get(nodeIndex).cutGenerated = true;
					DD.nodes.get(nodeIndex).folResponse = folObjValR;
					if (folObjValR > folUB) {folUB = folObjValR;}
					//System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
				}
				else {
					DD.arcs.get(i).cost = Math.min(DD.nodes.get(nodeIndex).folResponse, folObjM);
				}
			}
		}
		// Propagate costs backwards
		for (int i = DD.nodes.size()-1; i >= 2 ; i--) { // 0 and 1 are root node and end node
			DD.nodes.get(i).propagateCost(data);
		}
		System.out.println(" FOLL OBJ UB: "+folUB);

		int timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		double lim = Math.max( 10, timeLim-timeUsed); //For safety!
		cplex.setParam(IloCplex.Param.TimeLimit, lim);
		if(isExact==false) {cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0.1);}		//Don't need to close the gap at the beginning, then switch to 0 gap
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		// block variables for the new VIS
		IloNumVar[][] r = new IloNumVar[data.P.nConsF][data.P.nVarF]; 
		for (int i = 0; i < data.P.nVarF; i++) {
			for(int j = 0; j < data.P.nConsF; j++)
				if(data.cInteraction[j].equals("F") == false) {//block directly with x_i!
					r[j][i] = cplex.numVar(0, 0, IloNumVarType.Float, "r");
				}
				else {//Here we can block
					r[j][i] = cplex.numVar(0, 1, IloNumVarType.Float, "r");
				}
		}
		// Flow variables for the BDDs
		IloNumVar[] w = new IloNumVar[DD.arcs.size()]; 
		for (int i = 0; i < DD.arcs.size(); i++) {
//			w[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "w");
			w[i] = cplex.numVar(0, 1, IloNumVarType.Float, "w");
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}
		if(logical) {
			// New VIs: each follower variable with negative obj coefficient set to zero only if blocked 
			int[] varRecord = new int[data.P.nVarF];  //Keep track of which variables are already represented in the Sample
			for (int q = 0; q < sampleMaster.size(); q++) {
				for(int i=0; i < data.P.nVarF; i++){
					if(sampleMaster.get(q).y[i] > 0.9) {
						varRecord[i]++;
					}
				}
			}
			System.out.println("Sample var record: "+Arrays.toString(varRecord));
			for (int i = 0; i <data.P.nVarF; i++) {
				if(data.P.cfol[i] < 0 && varRecord[i] <= 10.5*sampleMaster.size() ) {//Only add for variables likely to be part of opt response
					System.out.println("Adding VIs for var "+i+" record "+varRecord[i]);
					IloLinearNumExpr exprf = cplex.linearNumExpr();
					exprf.addTerm(1, y[i]);
					exprf.addTerm(1, x[i]);	//Block it directly
					for(int j=0; j < data.P.nConsF; j++){//Block it by some follower constraint
						if(data.cInteraction[j].equals("F")) {exprf.addTerm(1, r[j][i]);}	
					}
					cplex.addGe(exprf, 1);
					// Now define the blocking variables
					for(int j=0; j < data.P.nConsF; j++){
						double block = data.P.bfol[j]-data.P.Afol[j][i]+1;
						IloLinearNumExpr exprB = cplex.linearNumExpr();
						for (int q = 0; q <data.P.nVarF; q++) {
							if(q!=i) {exprB.addTerm(data.P.Afol[j][q], y[q]);}
						}
						for (int q = 0; q <data.P.nVarL; q++) {
							exprB.addTerm(data.P.Bfol[j][q], x[q]);
						}
						exprB.addTerm(-block, r[j][i]);
						//System.out.println("CONSTRAINT "+j+" block: "+block);
						cplex.addGe(exprB, 0);	
					}
				}
			}
		}
		// FLOW CONSTRAINTS	
		// Add flow constraints (but drop the last one)
		IloLinearNumExpr[] expr = new IloLinearNumExpr[DD.nodes.size()];
		for (int n = 0; n < DD.nodes.size(); n++) {
			expr[n] = cplex.linearNumExpr();
		}	
		for (int i = 0; i < DD.arcs.size(); i++) {
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			expr[tail].addTerm(1, w[i]);
			expr[head].addTerm(-1, w[i]);
			//System.out.println("FLOW FROM "+DD.nodes.get(tail).key+" to "+DD.nodes.get(head).key+"        "+tail+" to "+head+" X"+DD.arcs.get(i).varIndex+" val: "+DD.arcs.get(i).arcVal);
		}
		cplex.addEq(expr[0], 1, "c0");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){cplex.addEq(expr[n], 0, "flow_"+0+"_"+n);}
		}

		// Link leader variables
		IloLinearNumExpr[] exprLink = new IloLinearNumExpr[data.P.nVarL];
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i] = cplex.linearNumExpr();
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {	
				exprLink[DD.arcs.get(i).varIndex].addTerm(1, w[i]);
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i].addTerm(-1, x[i]);
			cplex.addEq(exprLink[i], 0);
		}

		//Add constraints for each element of the sample
		for (int q = 0; q < sampleMaster.size(); q++) {
			// Add constraint on theta
			IloLinearNumExpr exprV = cplex.linearNumExpr();
			for(int i=0; i < data.P.nVarF; i++){
				exprV.addTerm(y[i], data.P.cfol[i]);
				if(sampleMaster.get(q).y[i] > 0.9) {
					exprV.addTerm(x[i], data.P.cfol[i]);
				}
			}
			cplex.addLe(exprV, sampleMaster.get(q).folObj);
		}

		///////////////////////////////////////////////////// VALID INEQUALITIES ///////////////////////////////////////////
		if(isExact == false) {
			IloLinearNumExpr exprB = cplex.linearNumExpr();
			for(int q=0; q<data.P.nVarF; q++){
				exprB.addTerm(y[q], data.P.cfol[q]);
			}
			cplex.addLe(exprB, folUB);
		}
		/////////////////////////////////////////////////// DD VALUE FUNCTION CONSTRAINT/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			exprV.addTerm(w[i], -DD.arcs.get(i).cost);
		}
		cplex.addLe(exprV, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}
		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		boolean stop = false;
		int[] blockRecord = new int[data.P.nVarL];  //Keep track of which constraints we are using to block
		int iter = 0;
		while(UB-LB >0 && !stop) {
			cplex.solve();
			if(iter >= 250 || timeLim-timeUsed <= 900  || Math.abs((UB-LB)/(UB+0.0)) <= 0.05 || (cplex.getStatus() == IloCplex.Status.Feasible || cplex.getStatus() == IloCplex.Status.Optimal) && round(cplex.getObjValue()) == UB ) { //Back to default to solve to optimality!!				
				cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-6);
				System.out.println("READY TO GO OPTIMAL IN THE MASTER PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			}

			int numSolutions = cplex.getSolnPoolNsolns();
			System.out.println("Number of solutions recovered: "+numSolutions);
			if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
				System.out.println("		DD-IP-R is infeasible!");
				LB = UB; //If UB= 99999 then problem is infeasible, otherwise cplex just closed the gap!
				stop = true;
			}
			else if(cplex.getStatus() == IloCplex.Status.Optimal) {
				LB = (int) Math.max(Math.round(cplex.getBestObjValue()), LB);
				System.out.println("******************************************************  MIP ITERATION "+iter+" LB: "+ LB+" UB: "+UB);

				// Recover solution(s)
				int[][] xVal = new int[numSolutions][data.P.nVarL];
				int[] yVal = new int[data.P.nVarF];
				for (int k = 0; k < numSolutions; k++) {
					for (int i = 0; i < data.P.nVarL; i++) {				
						if(cplex.getValue(x[i],k)>0.9) {xVal[k][i] = 1;}
					}
				}
				for (int i = 0; i <data.P.nVarF; i++) {
					yVal[i] = (int)Math.round(cplex.getValue(y[i]));
				}
				System.out.println("x: "+Arrays.toString(xVal[0]));
				System.out.println("y: "+Arrays.toString(yVal));

				// Try to find feasible completion for the leader variables
				Sample s = solveFollowerSample(data, getS(data, xVal[0]));
				sampleMaster.add(s);
				if(getFolObj(data, yVal) <= s.folObj) {//Current solution is bilevel feasible!!!!
					if(cplex.getObjValue() < UB) {
						UB = (int)Math.round(cplex.getObjValue());
						System.out.println("Current solution is bilevel feasible!!! LB:"+ LB+" UB: "+UB);
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[0][i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = yVal[i];
						}
					}
				}
				else {
					//Try to find feasible solution
					int leaderObjVal = solveFixedLeader(data, xVal[0], s);
					if(leaderObjVal < Integer.MAX_VALUE && leaderObjVal < UB) {
						System.out.println("UB UPDATED TO "+leaderObjVal);
						UB = leaderObjVal;
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[0][i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = s.yopt[i];
						}
					}
				}

				// Add blocking constraints interdiction style
				IloLinearNumExpr exprp = cplex.linearNumExpr();
				for(int i=0; i < data.P.nVarF; i++){
					exprp.addTerm(y[i], data.P.cfol[i]);
					if(s.y[i] > 0.9) {
						exprp.addTerm(x[i], data.P.cfol[i]);
					}
				}
				cplex.addLe(exprp, s.folObj);
				numCuts++;
				// Add blocking constraints interdiction style for other solutions in the pool
				for (int k = 1; k < numSolutions; k++) {
					// Solve follower
					Sample sprime = solveFollowerSample(data, getS(data, xVal[k]));
					IloLinearNumExpr exprE = cplex.linearNumExpr();
					for(int i=0; i < data.P.nVarF; i++){
						exprE.addTerm(y[i], data.P.cfol[i]);
						if(sprime.y[i] > 0.9) {
							exprE.addTerm(x[i], data.P.cfol[i]);
						}
					}
					cplex.addLe(exprE, sprime.folObj);
					System.out.println("Solution "+k+" val "+sprime.folObj);
					numCuts++;

					//Try to find feasible solution
					int leaderObjVal = solveFixedLeader(data, xVal[k], sprime);
					if(leaderObjVal < Integer.MAX_VALUE && leaderObjVal < UB) {
						System.out.println("UB UPDATED TO "+leaderObjVal);
						UB = leaderObjVal;
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[k][i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = sprime.yopt[i];
						}
					}
				}

			}
			else {
				System.out.println("ABORTED BECAUSE OF TIME LIMIT ");
				LB = (int) Math.max(Math.round(cplex.getBestObjValue()), LB); 
				if (cplex.getStatus() == IloCplex.Status.Feasible) {
					// Recover solution
					int[] xVal = new int[data.P.nVarL];
					int[] yVal = new int[data.P.nVarF];
					for (int i = 0; i < data.P.nVarL; i++) {
						xVal[i] = (int)Math.round(cplex.getValue(x[i]));
					}
					for (int i = 0; i <data.P.nVarF; i++) {
						yVal[i] = (int)Math.round(cplex.getValue(y[i]));
					}

					// Try to find feasible completion for the leader variables
					Sample s = solveFollowerSample(data, getS(data, xVal));
					if(getFolObj(data, yVal) <= s.folObj && cplex.getObjValue() < UB) {
						UB = (int) Math.round(cplex.getObjValue());
						// Record solution
						xstar = new int[data.P.nVarL];
						ystar = new int[data.P.nVarF];
						for (int i = 0; i < data.P.nVarL; i++) {
							xstar[i] = xVal[i];
						}
						for (int i = 0; i < data.P.nVarF; i++) {
							ystar[i] = yVal[i];
						}

					}
					else {
						//Try to find feasible solution
						int leaderObjVal = solveFixedLeader(data, xVal, s);
						if(leaderObjVal < Integer.MAX_VALUE && leaderObjVal < UB) {
							System.out.println("UB UPDATED TO "+leaderObjVal);
							UB = leaderObjVal;
							// Record solution
							xstar = new int[data.P.nVarL];
							ystar = new int[data.P.nVarF];
							for (int i = 0; i < data.P.nVarL; i++) {
								xstar[i] = xVal[i];
							}
							for (int i = 0; i < data.P.nVarF; i++) {
								ystar[i] = s.yopt[i];
							}
						}
					}
				}
				stop = true;
			}
			if(stop == false) {
				timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
				if(timeUsed >= timeLim-10) {stop = true;}
				else {
					cplex.setParam(IloCplex.Param.TimeLimit, timeLim-timeUsed);
					System.out.println("TIME USED: "+timeUsed);
				}
			}
			iter++;
		}


	}

	public void DDRelaxation(DataHandler data, int timeLim, double Atime) throws IloException {
		// Don't add any cuts!!!
		sampleMaster.clear();
		sampleMap.clear();

		// Solve the follower problem for each node/arc in the last layer (not solved before)
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				if(DD.nodes.get(nodeIndex).cutGenerated == false) {
					//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
					int folObjValR = solveFollower(data, DD.nodes.get(nodeIndex).S); 
					DD.arcs.get(i).cost = folObjValR;
					DD.nodes.get(nodeIndex).cutGenerated = true;
					DD.nodes.get(nodeIndex).folResponse = folObjValR;
					if (folObjValR > folUB) {folUB = folObjValR;}
					//System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
				}
				else {
					DD.arcs.get(i).cost = Math.min(DD.nodes.get(nodeIndex).folResponse, folObjM);
				}
			}
		}
		// Propagate costs backwards
		for (int i = DD.nodes.size()-1; i >= 2 ; i--) { // 0 and 1 are root node and end node
			DD.nodes.get(i).propagateCost(data);
		}
		System.out.println(" FOLL OBJ UB: "+folUB);

		int timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setParam(IloCplex.Param.TimeLimit, timeLim-timeUsed);
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}

		// Flow variables for the BDDs
		IloNumVar[] w = new IloNumVar[DD.arcs.size()]; 
		for (int i = 0; i < DD.arcs.size(); i++) {
			w[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "w");
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		// FLOW CONSTRAINTS	
		// Add flow constraints (but drop the last one)
		IloLinearNumExpr[] expr = new IloLinearNumExpr[DD.nodes.size()];
		for (int n = 0; n < DD.nodes.size(); n++) {
			expr[n] = cplex.linearNumExpr();
		}	
		for (int i = 0; i < DD.arcs.size(); i++) {
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			expr[tail].addTerm(1, w[i]);
			expr[head].addTerm(-1, w[i]);
			//System.out.println("FLOW FROM "+DD.nodes.get(tail).key+" to "+DD.nodes.get(head).key+"        "+tail+" to "+head+" X"+DD.arcs.get(i).varIndex+" val: "+DD.arcs.get(i).arcVal);
		}
		cplex.addEq(expr[0], 1, "c0");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){cplex.addEq(expr[n], 0, "flow_"+0+"_"+n);}
		}

		// Link leader variables
		IloLinearNumExpr[] exprLink = new IloLinearNumExpr[data.P.nVarL];
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i] = cplex.linearNumExpr();
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {	
				exprLink[DD.arcs.get(i).varIndex].addTerm(1, w[i]);
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i].addTerm(-1, x[i]);
			cplex.addEq(exprLink[i], 0);
		}

		///////////////////////////////////////////////////// VALID INEQUALITIES ///////////////////////////////////////////
		IloLinearNumExpr exprB = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprB.addTerm(y[q], data.P.cfol[q]);
		}
		cplex.addLe(exprB, folUB);

		/////////////////////////////////////////////////// DD VALUE FUNCTION CONSTRAINT/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			exprV.addTerm(w[i], -DD.arcs.get(i).cost);
		}
		cplex.addLe(exprV, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}
		//cplex.setOut(null);
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0);
		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Relaxation is infeasible!");
			LB = 999999999; // problem is infeasible
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			LB = (int) Math.round(cplex.getObjValue());
			System.out.println("******************************************************  DD RELAXED LB: "+ LB);
		}
		else {
			System.out.println("ABORTED BECAUSE OF TIME LIMIT ");
			LB = (int) cplex.getBestObjValue(); 
		}

	}



	// Relaxation for the refine method
	public Solution solveRelaxationDD(DataHandler data) throws IloException {
		// Setup the arc costs
		DD.resetArcCosts();
		// Solve the follower problem for each node/arc in the last layer (not solved before)
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				if(DD.nodes.get(nodeIndex).cutGenerated == false) {
					//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
					int folObjValR = solveFollower(data, DD.nodes.get(nodeIndex).S); 
					DD.arcs.get(i).cost = folObjValR;
					DD.nodes.get(nodeIndex).cutGenerated = true;
					DD.nodes.get(nodeIndex).folResponse = folObjValR;
					if (folObjValR > folUB) {folUB = folObjValR;}
					System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
				}
				else {
					DD.arcs.get(i).cost = DD.nodes.get(nodeIndex).folResponse;
				}
			}
		}
		// Propagate costs backwards
		for (int i = data.P.nVarL; i >= 1; i--) {
			for (int j = 0; j < DD.layers[i].size(); j++) {
				DD.layers[i].get(j).propagateCost(data); 
			}	
		}
		//DD.print3(data);
		System.out.println(" FOLL OBJ LB: "+folUB);

		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setOut(null);

		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}

		// Flow variables for the BDDs
		IloNumVar[] w = new IloNumVar[DD.arcs.size()]; 
		for (int i = 0; i < DD.arcs.size(); i++) {
			w[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}

		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		// FLOW CONSTRAINTS	
		// Add flow constraints (but drop the last one)
		IloLinearNumExpr[] expr = new IloLinearNumExpr[DD.nodes.size()];
		for (int n = 0; n < DD.nodes.size(); n++) {
			expr[n] = cplex.linearNumExpr();
		}	
		for (int i = 0; i < DD.arcs.size(); i++) {
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			expr[tail].addTerm(1, w[i]);
			expr[head].addTerm(-1, w[i]);	
		}
		cplex.addEq(expr[0], 1, "c0");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){cplex.addEq(expr[n], 0, "flow_"+0+"_"+n);}
		}

		// Link leader variables
		IloLinearNumExpr[] exprLink = new IloLinearNumExpr[data.P.nVarL];
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i] = cplex.linearNumExpr();
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {	
				exprLink[DD.arcs.get(i).varIndex].addTerm(1, w[i]);
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i].addTerm(-1, x[i]);
			cplex.addEq(exprLink[i], 0);
		}

		///////////////////////////////////////////////////// VALID INEQUALITIES ///////////////////////////////////////////
		IloLinearNumExpr exprB = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprB.addTerm(y[q], data.P.cfol[q]);
		}
		cplex.addLe(exprB, folUB);

		/////////////////////////////////////////////////// DD VALUE FUNCTION CONSTRAINT/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			exprV.addTerm(w[i], -DD.arcs.get(i).cost);
		}
		cplex.addLe(exprV, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}

		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		cplex.solve();
		Solution sol = new Solution (data);
		int[] ybar = new int[data.P.nVarF];
		if(cplex.getStatus() == IloCplex.Status.Optimal) {
			LB = (int) Math.round(cplex.getObjValue());
			System.out.println("***********************************  Relaxation LB: "+ LB);
			for (int i = 0; i < data.P.nVarL; i++) {
				if(cplex.getValue(x[i])>0.1) {
					sol.xbar[i] = 1;
				}
			}
			for (int i = 0; i < data.P.nVarF; i++) {
				if(cplex.getValue(y[i])>0.1) {
					ybar[i] = 1;
				}
			}

			for (int i = 0; i < DD.arcs.size(); i++) {
				if(cplex.getValue(w[i])>0.1) {
					Arc a = DD.arcs.get(i);
					sol.pathNodes[DD.nodes.get(a.tail).layer] = DD.nodes.get(a.tail);
					//System.out.println("ARC FROM "+DD.nodes.get(a.tail).key+" to "+DD.nodes.get(a.head).key+" val "+a.arcVal+" cost "+a.cost);
				}
			}

			// Check if current solution is bilevel feasible
			int folObjVal = solveFollower(data, getS(data, sol.xbar));
			if(getFolObj(data, ybar) <= folObjVal) {//Current solution is bilevel feasible!!!!
				System.out.println("Current solution is bilevel feasible!");
				UB = LB;
			}

		}
		cplex.clearModel();
		cplex.end();
		cplex = null;
		System.gc();

		return sol;
	}

	// Solves high point relaxation
	public void solveHPoint(DataHandler data, int timeLim) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setParam(IloCplex.Param.TimeLimit, timeLim);

		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}

		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}

		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		HPR is infeasible!");
			LB = 999999999; // problem is infeasible
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			LB = (int) Math.round(cplex.getObjValue());
			System.out.println("***********************************  High point zstar: "+ LB);
		}


	}

	public void solveHPointLR(DataHandler data, int timeLimit) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);

		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Float, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Float, "y["+i+"]");
		}

		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}

		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		HP is infeasible!");
			LB = 999999999; // problem is infeasible
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			LR = cplex.getObjValue();
			System.out.println("***********************************  High point LR: "+ LR);
		}


	}

	private int computeStateBound(DataHandler data, int index) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setOut(null);
		//Variables 
		IloNumVar[] xL = new IloNumVar[data.P.nVarL];
		IloNumVar[] xF = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarF; i++) {
			xF[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			xL[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], xL[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], xF[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], xF[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], xL[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(xL[i], data.P.Bfol[index][i]);
		}
		cplex.addMaximize(obj,"RHS");

		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		State bound compute is infeasible!");	
			return -1;
		}
		else{
			//System.out.println("INDEX: "+index+" State Bound: "+ (int) Math.round(cplex.getObjValue()));
			int[] xbar = new int[data.P.nVarL];
			for (int i = 0; i < data.P.nVarL; i++) {
				if(cplex.getValue(xL[i])>0.1) {
					xbar[i] = 1;
				}
			}
			// Add response to the sample
			Sample s = solveFollowerSample(data, getS(data, xbar));
			sampleAux.add(s);
			return (int) Math.round(cplex.getObjValue());

		}


	}

	public void checkStarSol(DataHandler data) throws IloException {
		if(UB == 999999999) {
			System.out.println("NO FEASIBLE SOLUTION FOUND!!!!");
		}
		else {
			System.out.println("**************************** SOLUTION CHECK ***********************************************");
			int folObj = 0;
			int leaderObj = 0;
			for (int i = 0; i < data.P.nVarL; i++) {
				leaderObj += data.P.c1leader[i]*xstar[i];
			}
			for (int i = 0; i < data.P.nVarF; i++) {
				folObj += data.P.cfol[i]*ystar[i];
				leaderObj += data.P.c2leader[i]*ystar[i];
			}
			System.out.println("LEADER OBJ: "+leaderObj);
			System.out.println("FOLLOWER OBJ: "+folObj);

			// Solve for fixed xstar
			cplex = new IloCplex();
			cplex.setParam(IloCplex.Param.Threads, numCores);
			cplex.setOut(null);
			//Variables 
			IloNumVar[] xL = new IloNumVar[data.P.nVarL];
			IloNumVar[] xF = new IloNumVar[data.P.nVarF];
			for (int i = 0; i < data.P.nVarL; i++) {
				xL[i] = cplex.numVar(xstar[i], xstar[i], IloNumVarType.Bool, "xL["+i+"]");
			}
			for (int i = 0; i < data.P.nVarF; i++) {
				xF[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "xF["+i+"]");
			}

			// add constraints follower
			for(int j=0; j < data.P.nConsF; j++){
				IloLinearNumExpr exprf = cplex.linearNumExpr();
				for (int i = 0; i <data.P.nVarL; i++) {
					exprf.addTerm(data.P.Bfol[j][i], xL[i]);
				}
				for (int i = 0; i <data.P.nVarF; i++) {
					exprf.addTerm(data.P.Afol[j][i], xF[i]);
				}
				cplex.addLe(exprf, data.P.bfol[j]);
			}

			// ADD OBJECTIVE FUNCTION	
			IloLinearNumExpr obj = cplex.linearNumExpr();

			for(int i=0; i<data.P.nVarF; i++){
				obj.addTerm(xF[i], data.P.cfol[i]);
			}
			cplex.addMinimize(obj,"FollObj");

			// Optimize model
			cplex.solve();

			if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
				System.out.println("		Check Solution is infeasible!");	
			}
			else{
				System.out.println("Follower optimal response obj: "+ round(cplex.getObjValue()));
				if(Math.abs(folObj - round(cplex.getObjValue())) > 0.0001)
					System.out.println("WARNING FOLLOWER OBJECTIVES MISMATCH: "+Math.abs(folObj - round(cplex.getObjValue())));
			}
		}

	}

	public int solveFollower(DataHandler data, int[] S) throws IloException {
		IloCplex cplexA;							// Cplex aux model
		cplexA = new IloCplex();
		cplexA.setParam(IloCplex.Param.Threads, numCores);
		cplexA.setOut(null);
		//Variables 
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplexA.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexA.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			cplexA.addLe(exprf, data.P.bfol[j]-S[j]); //Reduce the RHS according to node state
			//System.out.println("REMAINING RHS: "+(data.P.bfol[j]-S[j]));
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplexA.linearNumExpr();
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.cfol[i]);
		}
		cplexA.addMinimize(obj,"FollObj");

		// Optimize model
		cplexA.solve();
		int fObj = folObjM;
		if(cplexA.getStatus() == IloCplex.Status.Infeasible || cplexA.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Follower response is infeasible!");
			return fObj; 
		}
		else{
			fObj = (int) Math.round(cplexA.getObjValue());
			//System.out.println("Follower optimal response obj: "+fObj);
			/*for (int i = 0; i < data.nVar; i++) {
				if(cplexA.getValue(y[i])>0.1) {
					System.out.print(i+" ");
				}
			}
			System.out.println();*/
		}

		cplexA.clearModel();
		cplexA.end();
		cplexA = null;
		System.gc();
		return fObj;
	}

	public Sample solveFollowerApprox(DataHandler data, int[] S, IloRange[] A) throws IloException {
		//Update rhs
		for(int j=0; j < data.P.nConsF; j++){
			A[j].setUB(data.P.bfol[j]-capS(data, S)[j]); //Reduce the RHS according to node state
		}
		// Optimize model
		cplexAP.solve();
		int fObj = folObjM;

		if(cplexAP.getStatus() == IloCplex.Status.Infeasible || cplexAP.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Follower response is infeasible!");
			Sample s = new Sample(data, new int[data.P.nVarF], fObj);
			return s; 
		}
		else{
			fObj = (int) Math.round(cplexAP.getObjValue());
			int[] yval = new int[data.P.nVarF];
			//System.out.println("Follower optimal response obj: "+fObj);
			for (int i = 0; i < data.P.nVarF; i++) {
				if(cplexAP.getValue(yAP[i])>0.9) {
					yval[i] = 1;
				}
			}
			Sample s = new Sample(data, yval, fObj);
			return s;
		}

	}

	public int solveFixedLeader(DataHandler data, int[] xVal, Sample s) throws IloException {
		// Solve for fixed xVal
		IloCplex cplexA;							// Cplex aux model
		cplexA = new IloCplex();
		cplexA.setParam(IloCplex.Param.Threads, numCores);
		cplexA.setOut(null);
		//Variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplexA.numVar(xVal[i], xVal[i], IloNumVarType.Bool, "xL["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplexA.numVar(0, 1, IloNumVarType.Bool, "xF["+i+"]");
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplexA.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplexA.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexA.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplexA.addLe(exprf, data.P.bfol[j]);
		}
		// Add value function constraint
		IloLinearNumExpr expr = cplexA.linearNumExpr();
		for(int i=0; i<data.P.nVarF; i++){
			expr.addTerm(y[i], data.P.cfol[i]);
		}
		cplexA.addLe(expr, s.folObj);

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplexA.linearNumExpr();

		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}

		cplexA.addMinimize(obj,"LeaderObj");

		// Optimize model
		cplexA.solve();
		int ub = Integer.MAX_VALUE;
		if(cplexA.getStatus() == IloCplex.Status.Infeasible || cplexA.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			//System.out.println("		Solution completion is infeasible!");
			return ub;
		}
		else{
			ub = (int) Math.round(cplexA.getObjValue());
			s.yopt = new int[data.P.nVarF]; 
			for (int i = 0; i < data.P.nVarF; i++) {
				s.yopt[i] = (int)Math.round(cplexA.getValue(y[i]));
			}
		}

		cplexA.clearModel();
		cplexA.end();
		cplexA = null;
		System.gc();
		return ub;
	}

	/**
	 * Evaluates the follower objective value for a follower solution.
	 */
	private int getFolObj(DataHandler data, int[] yVal) {
		int obj = 0;
		for (int i = 0; i < data.P.nVarF; i++) {
			obj+=data.P.cfol[i]*yVal[i];
		}
		return obj;
	}

	/**
	 * Computes the linking-state vector induced by a leader solution.
	 */
	private int[] getS(DataHandler data, int[] xVal) {
		int [] S = new int[data.P.nConsF];
		for (int i = 0; i < data.P.nVarL; i++) {
			for (int j = 0; j < data.P.nConsF; j++) {
				S[j] += xVal[i]*data.P.Bfol[j][i];
			}
		}
		//System.out.println("REAL S: "+Arrays.toString(S));
		return S;
	}

	/**
	 * Caps a state vector at the precomputed state bounds.
	 */
	private int[] capS(DataHandler data, int[] s) {
		int [] S = new int[data.P.nConsF];
		for (int j = 0; j < data.P.nConsF; j++) {
			S[j] = Math.min(s[j], data.stateBounds[j]);
		}

		return S;
	}

	/**
	 * Computes a valid upper bound on the follower objective.
	 */
	public void getM(DataHandler data) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		//cplex.setOut(null);
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		// add constraints leader
		for(int j = 0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}

		// add constraints follower
		for(int j = 0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.cfol[i]);
		}

		cplex.addMaximize(obj,"OBJ");

		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		HPP is infeasible!");
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			folObjM = (int) Math.round(cplex.getObjValue());
			System.out.println("***********************************  BigM: "+ folObjM);
		}

		cplex.clearModel();
		cplex.end();
		cplex = null;
		System.gc();

		//Compute bounds to prune DD
		MfolConstraint = new double[data.P.nConsF];
		System.out.print("STATE BOUNDS: [");
		for (int j = 0; j < data.P.nConsF; j++) {
			data.stateBounds[j] = Integer.MAX_VALUE;
			MfolConstraint[j] = -Integer.MAX_VALUE;
		}
		for (int j = 0; j < data.P.nConsF; j++) {
			if(data.cInteraction[j].equals("B")) {
				if(!data.isInterdiction) {
					data.stateBounds[j] = computeStateBound(data, j);
					MfolConstraint[j] = computeFolConstraintM(data, j);
				}
				else { //For interdiction problems any binary is bounded by 0 and 1! 
					data.stateBounds[j] = 1;
					MfolConstraint[j] = 0;
				}
			}
			System.out.print(" "+data.stateBounds[j]);
		}


		System.out.println("]");

	}

	public void reduceDD(DataHandler data, double atime) throws IloException {
		System.out.println("BEFORE REDUCING "+DD.nodes.size()+" "+DD.arcs.size());
		double timeMaxMinStart = System.currentTimeMillis();
		data.iniNodes = DD.nodes.size();
		data.iniArcs = DD.arcs.size();
		//Setup cplex
		IloCplex cplexA;							// Cplex aux model
		cplexA = new IloCplex();
		cplexA.setParam(IloCplex.Param.Threads, numCores);
		cplexA.setOut(null);
		//Variables 
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplexA.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		// add constraints follower
		IloRange[] A = new IloRange[data.P.nConsF];
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexA.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			A[j] = cplexA.addLe(exprf, 0); //Reduce the RHS according to node state
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplexA.linearNumExpr();
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.cfol[i]);
		}
		cplexA.addMinimize(obj,"FollObj");

		// Solve the follower problem for each node/arc in the last layer
		int countFeasible = 0;
		for (int i = 0; i < DD.layers[data.P.nVarL].size() ; i++) {
			Node n = DD.layers[data.P.nVarL].get(i);
			int folObjValR = folObjM;
			System.out.println("STATE: "+Arrays.toString(n.S));
			System.out.println("Lowrb: "+Arrays.toString(n.Llead));
			//Set RHS coefficients
			for(int j=0; j < data.P.nConsF; j++){
				A[j].setUB(data.P.bfol[j]-capS(data, n.S)[j]); //Reduce the RHS according to node state
			}
			// Optimize model
			cplexA.solve();
			if(cplexA.getStatus() == IloCplex.Status.Infeasible || cplexA.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
				if(!preproMaxout) {
					solveFollowerProject(data, capS(data, n.S));
					folObjValR = solveFollowerMaxMinS(data, n, atime);
				}

			}
			else if((int) Math.round(cplexA.getObjValue()) >= folObjM && !preproMaxout) {
				folObjValR = solveFollowerMaxMinS(data, n, atime);
			}
			else{
				folObjValR  = (int) Math.round(cplexA.getObjValue());
				countFeasible++;
			}
			n.cutGenerated = true;
			n.folResponse = folObjValR;
			if (folObjValR > folUB) {folUB= folObjValR;}
			System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
		}
		// Discard CPLEX
		cplexA.clearModel();
		cplexA.end();
		cplexA = null;
		System.gc();
		System.out.println("FEASIBLE FOLLOWER PROBLEMS: "+countFeasible);
		timeMaxMin = (System.currentTimeMillis()-timeMaxMinStart)/1000;
		// Reduce the last layer with special condition
		for (int i = 0; i < DD.layers[data.P.nVarL].size(); i++) {
			if(DD.layers[data.P.nVarL].get(i).removed==false) {
				for (int j = i+1; j < DD.layers[data.P.nVarL].size(); j++) {
					if(i!=j && DD.layers[data.P.nVarL].get(j).removed==false) {
						if(DD.layers[data.P.nVarL].get(i).folResponse == DD.layers[data.P.nVarL].get(j).folResponse) {
							//Remove node j
							DD.layers[data.P.nVarL].get(j).removed = true;
							//Update state 
							DD.layers[data.P.nVarL].get(i).reduceState(DD.layers[data.P.nVarL].get(j));
							//Redirect incoming arcs
							for (int a = 0; a < DD.layers[data.P.nVarL].get(j).in.size(); a++) {
								DD.layers[data.P.nVarL].get(j).in.get(a).head = DD.GraphMap.get(DD.layers[data.P.nVarL].get(i).key);
								DD.layers[data.P.nVarL].get(i).in.add(DD.layers[data.P.nVarL].get(j).in.get(a));

							}
							//Remove outgoing arcs
							for (int a = 0; a < DD.layers[data.P.nVarL].get(j).out.size(); a++) {
								DD.arcs.remove(DD.layers[data.P.nVarL].get(j).out.get(a));
							}

						}

					}

				}
			}

		}
		// Now reduce as normal!
		for (int layer = data.P.nVarL-1; layer >= 0; layer--) {
			//System.out.println("LAYER "+layer);
			for (int i = 0; i < DD.layers[layer].size(); i++) {
				if(DD.layers[layer].get(i).removed==false) {
					String key1 = DD.layers[layer].get(i).genOut();
					for (int j = i+1; j < DD.layers[layer].size(); j++) {
						if(i!=j && DD.layers[layer].get(j).removed==false) {
							String key2 = DD.layers[layer].get(j).genOut();
							if(key1.equals(key2)) {
								// Debug: identical keys found while comparing samples.
								//Remove node j
								DD.layers[layer].get(j).removed = true;
								//Update state 
								DD.layers[layer].get(i).reduceState(DD.layers[layer].get(j));
								//Redirect incoming arcs
								for (int a = 0; a < DD.layers[layer].get(j).in.size(); a++) {
									DD.layers[layer].get(j).in.get(a).head = DD.GraphMap.get(DD.layers[layer].get(i).key);
									DD.layers[layer].get(i).in.add(DD.layers[layer].get(j).in.get(a));
								}
								//Remove outgoing arcs
								for (int a = 0; a < DD.layers[layer].get(j).out.size(); a++) {
									DD.arcs.remove(DD.layers[layer].get(j).out.get(a));
								}
							}
						}
					}
				}
			}
		}

		//Remove nodes and update arc pointers
		Hashtable<Integer, Integer> IDMap = new Hashtable<Integer, Integer>();
		Hashtable<String, Integer> newMap = new Hashtable<String, Integer>();
		ArrayList<Node> newNodes = new ArrayList<Node>();
		for (int i = 0; i < DD.nodes.size(); i++) {
			if(DD.nodes.get(i).removed == false) {
				IDMap.put(DD.GraphMap.get(DD.nodes.get(i).key), newNodes.size());
				newMap.put(DD.nodes.get(i).key, newNodes.size());
				newNodes.add(DD.nodes.get(i));
			}
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			DD.arcs.get(i).tail = IDMap.get(DD.arcs.get(i).tail);
			DD.arcs.get(i).head = IDMap.get(DD.arcs.get(i).head);
		}
		DD.nodes.clear();
		DD.nodes = newNodes;
		DD.GraphMap.clear();
		DD.GraphMap = newMap;

		//Update layers
		DD.layers = new ArrayList[data.P.nVarL+1];
		for (int i = 0; i < data.P.nVarL+1; i++) {
			DD.layers[i] = new ArrayList<Node>();
		}
		DD.layers[0].add(DD.nodes.get(0)); //Add root!
		for (int i = 2; i < DD.nodes.size() ; i++) { // 0 and 1 are root node and end node
			DD.layers[DD.nodes.get(i).layer].add(DD.nodes.get(i));
		}

		//Merge samples
		for (int q= 0; q < sampleAux.size(); q++) {
			if(sampleMap.containsKey(sampleAux.get(q).key) == false) {
				sampleMaster.add(sampleAux.get(q));
				sampleMap.put(sampleAux.get(q).key, sampleMaster.size());
			}
		}
		if(isExact) {sampleMaster.clear();} // No need to keep any samples if exact!
		numCuts = sampleMaster.size();
		System.out.println("AFTER REDUCING "+DD.nodes.size()+" "+DD.arcs.size()+" SAMPLE SIZE: "+sampleMaster.size());

	}


	private void solveFollowerProject(DataHandler data, int[] S) throws IloException {
		System.out.println("SOLVING PROJECTION");
		IloCplex cplexB;							// Cplex aux model
		cplexB = new IloCplex();
		cplexB.setParam(IloCplex.Param.Threads, numCores);
		cplexB.setOut(null);
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		IloNumVar[] pi = new IloNumVar[data.P.nConsF];  // To relax the constraints if needed
		IloNumVar[] z = new IloNumVar[data.P.nConsF];  // To relax the constraints if needed
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplexB.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplexB.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nConsF; i++) {
			pi[i] = cplexB.numVar(0, Double.MAX_VALUE, IloNumVarType.Float);
			z[i] = cplexB.numVar(0, 1, IloNumVarType.Bool);
		}
		// HPP Polytope
		for(int j = 0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplexB.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplexB.addLe(exprf, data.P.b[j]);
		}
		for(int j = 0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexB.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplexB.addLe(exprf, data.P.bfol[j]);
		}
		// add constraints follower for the S given
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexB.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
				exprf.addTerm(-1, pi[j]);
			}
			cplexB.addLe(exprf, data.P.bfol[j]-S[j]); //Reduce the RHS according to node state
			//System.out.println("REMAINING RHS: "+(data.P.bfol[j]-S[j]));
		}

		// add activation
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexB.linearNumExpr();
			exprf.addTerm(1, pi[j]);
			exprf.addTerm(-99999, z[j]);
			cplexB.addLe(exprf, 0); 
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplexB.linearNumExpr();
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.cfol[i]);
		}
		for(int j=0; j < data.P.nConsF; j++){
			obj.addTerm(pi[j], 99);
			obj.addTerm(z[j],  99999);
		}
		cplexB.addMinimize(obj,"FollObj");
		// Optimize model
		cplexB.solve();
		if (cplexB.getStatus() == IloCplex.Status.Feasible || cplexB.getStatus() == IloCplex.Status.Optimal) {
			// Debug: inspect the objective value of the projected follower model.
			int fObj = 0;
			int[] yval = new int[data.P.nVarF];
			for (int i = 0; i < data.P.nVarF; i++) {
				if(cplexB.getValue(y[i])>0.1) {
					yval[i] = 1;
					fObj += data.P.cfol[i];
				}
			}
			Sample s = new Sample(data, yval, fObj);
			sampleAux.add(s);
			/*if(sampleMap.containsKey(s.key) == false) {
				sampleMaster.add(s);
				sampleMap.put(s.key, sampleMaster.size());
			}
			s.print();
			 */
		}
		cplexB.clearModel();
		cplexB.end();
		cplexB = null;
		System.gc();


	}

	private int solveFollowerMaxMinS(DataHandler data, Node node, double atime) throws IloException {
		System.out.println("FOLLOWER MAX MIN SAMPLE SIZE: "+sampleMaster.size());
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setParam(IloCplex.Param.TimeLimit, 300);
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 5e-2);
		cplex.setOut(null);

		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		IloNumVar[][] v = new IloNumVar[data.P.nConsF][sampleMaster.size()];
		IloNumVar theta = cplex.numVar(-Integer.MAX_VALUE, folObjM, IloNumVarType.Int, "theta");
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		for (int i = 0; i < data.P.nConsF; i++) {
			for(int j = 0; j < sampleMaster.size(); j++) {
				v[i][j] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
		}
		// block variables for the new VIS
		IloNumVar[][] r = new IloNumVar[data.P.nConsF][data.P.nVarF]; 
		for (int i = 0; i < data.P.nVarF; i++) {
			for(int j = 0; j < data.P.nConsF; j++)
				if(data.P.Afol[j][i] <= 0 ) {// Impossible to block!
					r[j][i] = cplex.numVar(0, 0, IloNumVarType.Float, "r");
				}
				else {//Here we can block
					r[j][i] = cplex.numVar(0, 1, IloNumVarType.Bool, "r");
				}
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}
		// New VIs: each follower binary variable with negative obj coefficient set to zero only if blocked 
		for (int i = 0; i <data.P.nVarF; i++) {
			if(data.P.cfol[i] < 0 ) {
				IloLinearNumExpr exprf = cplex.linearNumExpr();
				exprf.addTerm(1, y[i]);
				for(int j=0; j < data.P.nConsF; j++){
					exprf.addTerm(1, r[j][i]);	//Block it by some constraint
				}
				cplex.addGe(exprf, 1);
				// Now define the blocking variables
				for(int j=0; j < data.P.nConsF; j++){
					int block = (int) (data.P.bfol[j]-data.P.Afol[j][i]+1);
					IloLinearNumExpr exprB = cplex.linearNumExpr();
					for (int q = 0; q <data.P.nVarF; q++) {
						if(q!=i) {exprB.addTerm(data.P.Afol[j][q], y[q]);}
					}
					for (int q = 0; q <data.P.nVarL; q++) {
						exprB.addTerm(data.P.Bfol[j][q], x[q]);
					}
					exprB.addTerm(-block-999999999, r[j][i]);
					//System.out.println("CONSTRAINT "+j+" block: "+block);
					cplex.addGe(exprB, -999999999);	
				}
			}
		}

		// Link to the DD state
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, node.S[j]);
			cplex.addGe(exprf, node.Lfol[j]);
		}
		//Add constraints for each element of the sample
		for (int q = 0; q < sampleMaster.size(); q++) {
			//sample.get(q).print();
			// Add constraint on theta
			IloLinearNumExpr exprV = cplex.linearNumExpr();
			exprV.addTerm(theta, 1);
			for (int i = 0; i < sampleMaster.get(q).block.size(); i++) {
				exprV.addTerm(v[sampleMaster.get(q).block.get(i)[0]][q], -(folObjM-sampleMaster.get(q).folObj)); //
			}			
			cplex.addLe(exprV, sampleMaster.get(q).folObj);
			// Add block definitions
			for (int i = 0; i < sampleMaster.get(q).block.size(); i++) { 
				IloLinearNumExpr exprB = cplex.linearNumExpr();
				int index = sampleMaster.get(q).block.get(i)[0];
				exprB.addTerm(v[index][q], MfolConstraint[index]-sampleMaster.get(q).block.get(i)[1]);
				for (int k = 0; k <data.P.nVarL; k++) {
					exprB.addTerm(data.P.Bfol[sampleMaster.get(q).block.get(i)[0]][k], x[k]);
				}
				cplex.addGe(exprB, MfolConstraint[index]);
			}

		}

		// add theta constraint
		IloLinearNumExpr objC = cplex.linearNumExpr();
		objC.addTerm(theta, 1);
		for(int i=0; i<data.P.nVarF; i++){
			objC.addTerm(y[i], -data.P.cfol[i]);
		}
		cplex.addLe(objC, 0);

		// add UB constraint
		IloLinearNumExpr objD = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			objD.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			objD.addTerm(y[i], data.P.c2leader[i]);
		}
		cplex.addLe(objD, UB);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(theta, 1);
		cplex.addMaximize(obj,"OBJ");
		//cplex.solve();


		for (int q = 0; q < 10; q++) {
			// Optimize model
			cplex.solve();
			if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
				q=100;
			}
			else {
				System.out.println("ITERATION "+q+" THETA = "+cplex.getObjValue());
				// Recover x solution
				int[] xval = new int[data.P.nVarL];
				for (int i = 0; i < data.P.nVarL; i++) {
					if(cplex.getValue(x[i])>0.1) {
						xval[i] = 1;
					}
				}
				Sample s = solveFollowerSample(data, getS(data, xval));

				if(s!= null) {
					s.print();
					sampleAux.add(s);
					// Add constraint on theta
					IloNumVar[] vp = new IloNumVar[s.block.size()];
					for (int i = 0; i < s.block.size(); i++) {
						vp[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
					}
					IloLinearNumExpr exprp = cplex.linearNumExpr();
					exprp.addTerm(theta, 1);
					for (int i = 0; i < s.block.size(); i++) {
						exprp.addTerm(vp[i], -(folObjM-s.folObj)); //
					}			
					cplex.addLe(exprp, s.folObj);
					// Add block definitions
					for (int i = 0; i < s.block.size(); i++) { 
						IloLinearNumExpr exprB = cplex.linearNumExpr();
						exprB.addTerm(vp[i], MfolConstraint[s.block.get(i)[0]]-s.block.get(i)[1]);
						for (int k = 0; k <data.P.nVarL; k++) {
							exprB.addTerm(data.P.Bfol[s.block.get(i)[0]][k], x[k]);
						}
						cplex.addGe(exprB, MfolConstraint[s.block.get(i)[0]]);
					}
					if(data.isInterdiction) {//Add interdiction version
						IloLinearNumExpr exprI = cplex.linearNumExpr();
						exprI.addTerm(theta, 1);
						for(int i=0; i < data.P.nVarF; i++){
							if(s.y[i] > 0.9) {
								exprI.addTerm(x[i], data.P.cfol[i]);
							}
						}
						cplex.addLe(exprI, s.folObj);
					}
				}
			}
			double timeUsed = (System.currentTimeMillis()-atime)/1000.0;
			if(timeUsed >= 1200) {q=20; preproMaxout = true;}
		}
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		MAXMINDD is infeasible!");

		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("THETA: "+cplex.getValue(theta));
			return (int) Math.ceil(cplex.getBestObjValue());
		}
		else {
			System.out.println("THETA: "+cplex.getBestObjValue());
			return (int) Math.ceil(cplex.getBestObjValue());
		}

		cplex.clearModel();
		cplex.end();
		cplex = null;
		System.gc();

		return -999999999;   // This makes the DD node unreachable !!!		

	}




	public void refineMS(DataHandler data, int iter) throws IloException {
		double Atime = System.currentTimeMillis();
		//Initialize sample
		if(sampleMaster.size()== 0) {solveFollowerProject(data, data.stateBounds);}

		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setParam(IloCplex.Param.TimeLimit, 600);
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 5e-2);
		//cplex.setOut(null);
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		IloNumVar[][] w = new IloNumVar[data.P.nConsF][sampleMaster.size()];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		for (int i = 0; i < data.P.nConsF; i++) {
			for(int j = 0; j < sampleMaster.size(); j++) {
				w[i][j] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
		}
		IloNumVar theta = cplex.numVar(-Integer.MAX_VALUE, folObjM, IloNumVarType.Int, "theta"); 
		// block variables for the new VIS
		IloNumVar[][] r = new IloNumVar[data.P.nConsF][data.P.nVarF]; 
		for (int i = 0; i < data.P.nVarF; i++) {
			for(int j = 0; j < data.P.nConsF; j++)
				if(data.P.Afol[j][i] <= 0 ) {// Impossible to block!
					r[j][i] = cplex.numVar(0, 0, IloNumVarType.Float, "r");
				}
				else {//Here we can block
					r[j][i] = cplex.numVar(0, 1, IloNumVarType.Bool, "r");
				}
		}
		// add constraints leader
		for(int j = 0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j = 0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}
		// New VIs: each follower binary variable with negative obj coefficient set to zero only if blocked 
		/*		for (int i = 0; i <data.P.nVarF; i++) {
			if(data.P.cfol[i] < 0 ) {
				IloLinearNumExpr exprf = cplex.linearNumExpr();
				exprf.addTerm(1, y[i]);
				for(int j=0; j < data.P.nConsF; j++){
					exprf.addTerm(1, r[j][i]);	//Block it by some constraint
				}
				cplex.addGe(exprf, 1);
				// Now define the blocking variables
				for(int j=0; j < data.P.nConsF; j++){
					int block = (int) (data.P.bfol[j]-data.P.Afol[j][i]+1);
					IloLinearNumExpr exprB = cplex.linearNumExpr();
					for (int q = 0; q <data.P.nVarF; q++) {
						if(q!=i) {exprB.addTerm(data.P.Afol[j][q], y[q]);}
					}
					for (int q = 0; q <data.P.nVarL; q++) {
						exprB.addTerm(data.P.Bfol[j][q], x[q]);
					}
					exprB.addTerm(-block-999999999, r[j][i]);
					//System.out.println("CONSTRAINT "+j+" block: "+block);
					cplex.addGe(exprB, -999999999);	
				}
			}
		}*/
		//Add constraints for each element of the sample
		for (int q = 0; q < sampleMaster.size(); q++) {
			//sample.get(q).print();
			// Add constraint on theta
			IloLinearNumExpr expr = cplex.linearNumExpr();
			expr.addTerm(theta, 1);
			for (int i = 0; i < sampleMaster.get(q).block.size(); i++) {
				expr.addTerm(w[sampleMaster.get(q).block.get(i)[0]][q], -(folObjM-sampleMaster.get(q).folObj)); //
			}			
			cplex.addLe(expr, sampleMaster.get(q).folObj);
			// Add block definitions
			for (int i = 0; i < sampleMaster.get(q).block.size(); i++) { 
				IloLinearNumExpr exprB = cplex.linearNumExpr();
				int index = sampleMaster.get(q).block.get(i)[0];
				exprB.addTerm(w[index][q], MfolConstraint[index]-sampleMaster.get(q).block.get(i)[1]);
				//System.out.println("w"+s.block.get(i)[0]+"*"+s.block.get(i)[1]);
				for (int k = 0; k <data.P.nVarL; k++) {
					exprB.addTerm(data.P.Bfol[index][k], x[k]);
					//System.out.println("Constraint "+s.block.get(i)[0]+" Bij: "+data.P.Bfol[s.block.get(i)[0]][k]);
				}
				cplex.addGe(exprB, MfolConstraint[index]);
			}
			if(data.isInterdiction) {//Add interdiction version
				IloLinearNumExpr exprV = cplex.linearNumExpr();
				exprV.addTerm(theta, 1);
				for(int i=0; i < data.P.nVarF; i++){
					if(sampleMaster.get(q).y[i] > 0.9) {
						exprV.addTerm(x[i], data.P.cfol[i]);
					}
				}
				cplex.addLe(exprV, sampleMaster.get(q).folObj);
			}

		}

		// add theta constraint
		IloLinearNumExpr objC = cplex.linearNumExpr();
		objC.addTerm(theta, 1);
		for(int i=0; i<data.P.nVarF; i++){
			objC.addTerm(y[i], -data.P.cfol[i]);
		}
		cplex.addLe(objC, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		obj.addTerm(theta, 1);
		cplex.addMaximize(obj,"OBJ");

		for (int q = 0; q < iter; q++) {
			// Optimize model
			cplex.solve();
			System.out.println("ITERATION "+q+" THETA = "+cplex.getBestObjValue());
			// Recover x solution
			int[] xval = new int[data.P.nVarL];
			for (int i = 0; i < data.P.nVarL; i++) {
				if(cplex.getValue(x[i])>0.1) {
					xval[i] = 1;
				}
			}
			Sample s = solveFollowerSample(data, getS(data, xval));
			if(sampleMap.containsKey(s.key) == false) {
				sampleMaster.add(s);
				sampleMap.put(s.key, sampleMaster.size());
			}
			if(s!= null) {
				s.print();

				// Add constraint on theta
				IloNumVar[] v = new IloNumVar[s.block.size()];
				for (int i = 0; i < s.block.size(); i++) {
					v[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
				}
				IloLinearNumExpr expr = cplex.linearNumExpr();
				expr.addTerm(theta, 1);
				for (int i = 0; i < s.block.size(); i++) {
					expr.addTerm(v[i], -(folObjM-s.folObj)); //
				}			
				cplex.addLe(expr, s.folObj);
				// Add block definitions
				for (int i = 0; i < s.block.size(); i++) { 
					IloLinearNumExpr exprB = cplex.linearNumExpr();
					exprB.addTerm(v[i], MfolConstraint[s.block.get(i)[0]]-s.block.get(i)[1]);
					for (int k = 0; k <data.P.nVarL; k++) {
						exprB.addTerm(data.P.Bfol[s.block.get(i)[0]][k], x[k]);
					}
					cplex.addGe(exprB, MfolConstraint[s.block.get(i)[0]]);
				}

				if(data.isInterdiction) {//Add interdiction version
					IloLinearNumExpr exprV = cplex.linearNumExpr();
					exprV.addTerm(theta, 1);
					for(int i=0; i < data.P.nVarF; i++){
						if(s.y[i] > 0.9) {
							exprV.addTerm(x[i], data.P.cfol[i]);
						}
					}
					cplex.addLe(exprV, s.folObj);
				}
			}
			double timeUsed = (System.currentTimeMillis()-Atime)/1000.0;
			if(timeUsed >= 600) {q=iter; preproMaxout = true;}
		}
		cplex.solve();

		for (int q = 0; q < sampleMaster.size(); q++) {
			sampleMaster.get(q).print();
		}

		folObjM = (int) Math.ceil(cplex.getBestObjValue()); //Because we are not solving each problem to optimality

		System.out.println("UPDATED BIG M: "+folObjM);
		cplex.clearModel();
		cplex.end();
		cplex = null;
		System.gc();

	}

	public Sample solveFollowerSample(DataHandler data, int[] S) throws IloException {
		IloCplex cplexA;							// Cplex aux model
		cplexA = new IloCplex();
		cplexA.setParam(IloCplex.Param.Threads, numCores);
		cplexA.setOut(null);
		//Variables 
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplexA.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexA.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			cplexA.addLe(exprf, data.P.bfol[j]-S[j]); //Reduce the RHS according to node state
			//System.out.println("REMAINING RHS: "+(data.P.bfol[j]-S[j]));
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplexA.linearNumExpr();
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.cfol[i]);
		}
		cplexA.addMinimize(obj,"FollObj");
		// Optimize model
		cplexA.solve();
		if (cplexA.getStatus() == IloCplex.Status.Feasible || cplexA.getStatus() == IloCplex.Status.Optimal) {
			int fObj = 0;
			int[] yval = new int[data.P.nVarF];
			for (int i = 0; i < data.P.nVarF; i++) {
				if(cplexA.getValue(y[i])>0.1) {
					yval[i] = 1;
					fObj += data.P.cfol[i];
				}
			}
			Sample s = new Sample(data, yval, fObj);
			cplexA.clearModel();
			cplexA.end();
			cplexA = null;
			System.gc();
			return s;
		}
		cplexA.clearModel();
		cplexA.end();
		cplexA = null;
		System.gc();
		System.out.println("THIS SHOULD NEVER EVER HAPPEN!!!!");
		return null;
	}
	public Sample solveFollowerDiversityI(DataHandler data, int[] blockRecord) throws IloException {
		IloCplex cplexA;							// Cplex aux model
		cplexA = new IloCplex();
		cplexA.setParam(IloCplex.Param.Threads, numCores);
		cplexA.setOut(null);
		//Variables 
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplexA.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
			if( blockRecord[i] >= 1) {y[i] = cplexA.numVar(0, 0, IloNumVarType.Bool, "y["+i+"]");}
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplexA.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			cplexA.addLe(exprf, data.P.bfol[j]);
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplexA.linearNumExpr();
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.cfol[i]);
		}
		cplexA.addMinimize(obj,"FollObj");
		// Optimize model
		cplexA.solve();
		if (cplexA.getStatus() == IloCplex.Status.Feasible || cplexA.getStatus() == IloCplex.Status.Optimal) {
			int fObj = 0;
			int[] yval = new int[data.P.nVarF];
			for (int i = 0; i < data.P.nVarF; i++) {
				if(cplexA.getValue(y[i])>0.1) {
					yval[i] = 1;
					fObj += data.P.cfol[i];
				}
			}
			Sample s = new Sample(data, yval, fObj);
			return s;
		}
		cplexA.clearModel();
		cplexA.end();
		cplexA = null;
		System.gc();
		return null;
	}

	// Linear relaxation of the DD model
	public void singleLevelDDLR(DataHandler data, int timeLim, double Atime) throws IloException {
		// Solve the follower problem for each node/arc in the last layer (not solved before)
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				if(DD.nodes.get(nodeIndex).cutGenerated == false) {
					//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
					int folObjValR = solveFollower(data, capS(data,DD.nodes.get(nodeIndex).S)); 
					DD.arcs.get(i).cost = folObjValR;
					DD.nodes.get(nodeIndex).cutGenerated = true;
					DD.nodes.get(nodeIndex).folResponse = folObjValR;
					if (folObjValR > folUB) {folUB = folObjValR;}
					System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
				}
				else {
					DD.arcs.get(i).cost = Math.min(DD.nodes.get(nodeIndex).folResponse, folObjM);
				}
			}
		}
		// Propagate costs backwards
		//for (int i = DD.nodes.size()-1; i >= 2 ; i--) { // 0 and 1 are root node and end node
		//	DD.nodes.get(i).propagateCost(data);
		//}
		System.out.println(" FOLL OBJ UB: "+folUB);

		int timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setParam(IloCplex.Param.TimeLimit, timeLim-timeUsed);
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Float, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Float, "y["+i+"]");
		}

		// Flow variables for the BDDs
		IloNumVar[] w = new IloNumVar[DD.arcs.size()]; 
		for (int i = 0; i < DD.arcs.size(); i++) {
			w[i] = cplex.numVar(0, 1, IloNumVarType.Float, "w");
		}

		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		// FLOW CONSTRAINTS	
		// Add flow constraints (but drop the last one)
		IloLinearNumExpr[] expr = new IloLinearNumExpr[DD.nodes.size()];
		for (int n = 0; n < DD.nodes.size(); n++) {
			expr[n] = cplex.linearNumExpr();
		}	
		for (int i = 0; i < DD.arcs.size(); i++) {
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			expr[tail].addTerm(1, w[i]);
			expr[head].addTerm(-1, w[i]);
			//System.out.println("FLOW FROM "+DD.nodes.get(tail).key+" to "+DD.nodes.get(head).key+"        "+tail+" to "+head+" X"+DD.arcs.get(i).varIndex+" val: "+DD.arcs.get(i).arcVal);
		}
		cplex.addEq(expr[0], 1, "c0");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){cplex.addEq(expr[n], 0, "flow_"+0+"_"+n);}
		}

		// Link leader variables
		IloLinearNumExpr[] exprLink = new IloLinearNumExpr[data.P.nVarL];
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i] = cplex.linearNumExpr();
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {	
				exprLink[DD.arcs.get(i).varIndex].addTerm(1, w[i]);
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i].addTerm(-1, x[i]);
			cplex.addEq(exprLink[i], 0);
		}

		/////////////////////////////////////////////////// DD VALUE FUNCTION CONSTRAINT/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			exprV.addTerm(w[i], -DD.arcs.get(i).cost);
			// Debug: inspect the propagated arc cost.
		}
		cplex.addLe(exprV, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}
		//cplex.setOut(null);
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0);
		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		boolean stop = false;
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		DD-IP-R is infeasible!");
			LB = UB; //If UB= 99999 then problem is infeasible, otherwise cplex just closed the gap!
			stop = true;
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("****************************************************** LR VAL: "+ cplex.getObjValue());
			LR = cplex.getObjValue();

			for (int i = 0; i < data.P.nVarL; i++) {
				if(cplex.getValue(x[i])>0.001) {
					System.out.println("x"+i+" = "+cplex.getValue(x[i]));
				}
			}

			for (int i = 0; i < data.P.nVarF; i++) {
				if(cplex.getValue(y[i])>0.001) {
					System.out.println("y"+i+" = "+cplex.getValue(y[i]));
				}
			}

			for (int i = 0; i < DD.arcs.size(); i++) {
				if(cplex.getValue(w[i])>0.001) {
					Arc a = DD.arcs.get(i);
					System.out.println("ARC FROM "+DD.nodes.get(a.tail).key+" to "+DD.nodes.get(a.head).key+" flow: "+ cplex.getValue(w[i]));
				}
			}


		}	



	}

	// Linear relaxation of the binary model (no DDs!!!)
	public void singleLevelBinaryLR(DataHandler data, int timeLim, double Atime) throws IloException {
		//Capture all the states and their value function
		ArrayList<int[]> states = new ArrayList<int[]>();
		ArrayList<Integer> valFun = new ArrayList<Integer>();

		// Solve the follower problem for state in the last layer
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
				int folObjValR = solveFollower(data, DD.nodes.get(nodeIndex).S);
				states.add(DD.nodes.get(nodeIndex).S);
				valFun.add(folObjValR);
				System.out.println("State "+Arrays.toString(DD.nodes.get(nodeIndex).S)+" function "+folObjValR);
			}
		}

		int timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setParam(IloCplex.Param.TimeLimit, timeLim-timeUsed);
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Float, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Float, "y["+i+"]");
		}
		// State selector variables
		IloNumVar[] w = new IloNumVar[states.size()]; 
		for (int i = 0; i < states.size(); i++) {
			w[i] = cplex.numVar(0, 1, IloNumVarType.Float, "w");
		}

		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		/////////////////////////////////////////////////// VALUE FUNCTION CONSTRAINTS/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < states.size(); i++) {
			exprV.addTerm(w[i], -valFun.get(i));
		}
		cplex.addLe(exprV, 0);

		//Select only one state
		IloLinearNumExpr exprQ = cplex.linearNumExpr();
		for (int i = 0; i < states.size(); i++) {
			exprQ.addTerm(w[i] , 1);
		}
		cplex.addEq(exprQ, 1);

		//Link state
		for (int j = 0; j < data.P.nConsF; j++) {
			IloLinearNumExpr exprZ = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprZ.addTerm(data.P.Bfol[j][i], x[i]);
			}
			for (int i = 0; i < states.size(); i++) {
				exprZ.addTerm(w[i], -states.get(i)[j]);
			}
			cplex.addEq(exprZ, 0);
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}
		//cplex.setOut(null);
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0);
		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		boolean stop = false;
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Binary Model is infeasible!");
			LB = UB; //If UB= 99999 then problem is infeasible, otherwise cplex just closed the gap!
			stop = true;
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("****************************************************** LR VAL: "+ cplex.getObjValue());
			LR = cplex.getObjValue();
			for (int i = 0; i < states.size(); i++) {
				if(cplex.getValue(w[i])>0.01) System.out.println("w_"+i+" = "+cplex.getValue(w[i])+" state: "+Arrays.toString(states.get(i)));
			}
			for (int i = 0; i < data.P.nVarL; i++) {
				if(cplex.getValue(x[i])>0.01) System.out.println("x_"+(i+1)+" = "+cplex.getValue(x[i]));
			}
			for (int i = 0; i < data.P.nVarF; i++) {
				if(cplex.getValue(y[i])>0.01) System.out.println("y_"+(i+1)+" = "+cplex.getValue(y[i]));
			}


		}	



	}

	private double computeFolConstraintM(DataHandler data, int index) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setOut(null);
		//Variables 
		IloNumVar[] xL = new IloNumVar[data.P.nVarL];
		IloNumVar[] xF = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarF; i++) {
			xF[i] = cplex.numVar(0, 1, IloNumVarType.Float, "y["+i+"]");
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			xL[i] = cplex.numVar(0, 1, IloNumVarType.Float, "x["+i+"]");
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], xL[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], xF[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], xF[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], xL[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(xL[i], data.P.Bfol[index][i]);
		}
		cplex.addMinimize(obj,"RHS");

		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		FOL constraint M compute is infeasible!");	
			return -1;
		}
		else{
			//System.out.println("INDEX: "+index+" State Bound: "+ (int) Math.round(cplex.getObjValue()));
			return cplex.getObjValue();

		}


	}

	/**
	 * Exports the incumbent solution in the JSON format used by the replication
	 * package.
	 */
	public void ExportSol(DataHandler data, String name, PrintStream ps) {
		if(UB == 999999999) {
			System.out.println("NO FEASIBLE SOLUTION FOUND!!!!");
		}
		else {
			ps.println("{");
			ps.println("\t\"name\": \""+name+"\",");
			ps.println("\t\"bilevel_type\": \"optimistic\",");
			ps.println("\t\"status\": \"\",");
			ps.println("\t\"difficulty\": \"\",");
			ps.println("\t\"objective_value\": "+UB+",");
			ps.println("\t\"upper_level_decisions\": {");
			for (int i = 0; i < data.P.nVarL-1; i++) {
				ps.println("\t\t\""+data.P.leaderVarNames[i]+"\": "+xstar[i]+",");
			}
			ps.println("\t\t\""+data.P.leaderVarNames[data.P.nVarL-1]+"\": "+xstar[data.P.nVarL-1]);
			ps.println("\t},");
			ps.println("\t\"lower_level_decisions\": {");
			for (int i = 0; i < data.P.nVarF-1; i++) {
				ps.println("\t\t\""+data.P.followerVarNames[i]+"\": "+ystar[i]+",");
			}
			ps.println("\t\t\""+data.P.followerVarNames[data.P.nVarF-1]+"\": "+ystar[data.P.nVarF-1]);
			ps.println("\t}");
			ps.println("}");
		}
	}

	/**
	 * Solves the interdiction root-node reformulation with blocking constraints,
	 * without callbacks.
	 */
	public void singleLevelDDInterdictionRoot(DataHandler data, int timeLim, double Atime, String file) throws IloException {
		// Don't add any cuts!!!
		sampleMaster.clear();
		sampleMap.clear();

		// Solve the follower problem for each node/arc in the last layer (not solved before)
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				if(DD.nodes.get(nodeIndex).cutGenerated == false) {
					//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
					int folObjValR = solveFollower(data, DD.nodes.get(nodeIndex).S); 
					DD.arcs.get(i).cost = folObjValR;
					DD.nodes.get(nodeIndex).cutGenerated = true;
					DD.nodes.get(nodeIndex).folResponse = folObjValR;
					if (folObjValR > folUB) {folUB = folObjValR;}
					//System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
				}
				else {
					DD.arcs.get(i).cost = Math.min(DD.nodes.get(nodeIndex).folResponse, folObjM);
				}
			}
		}
		// Propagate costs backwards
		for (int i = DD.nodes.size()-1; i >= 2 ; i--) { // 0 and 1 are root node and end node
			DD.nodes.get(i).propagateCost(data);
		}
		System.out.println(" FOLL OBJ UB: "+folUB);

		int timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		double lim = Math.max( 10, timeLim-timeUsed); //For safety!
		cplex.setParam(IloCplex.Param.TimeLimit, lim);
		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}

		// Flow variables for the BDDs
		IloNumVar[] w = new IloNumVar[DD.arcs.size()]; 
		for (int i = 0; i < DD.arcs.size(); i++) {
			w[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "w");
		}
		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}

		// FLOW CONSTRAINTS	
		// Add flow constraints (but drop the last one)
		IloLinearNumExpr[] expr = new IloLinearNumExpr[DD.nodes.size()];
		for (int n = 0; n < DD.nodes.size(); n++) {
			expr[n] = cplex.linearNumExpr();
		}	
		for (int i = 0; i < DD.arcs.size(); i++) {
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			expr[tail].addTerm(1, w[i]);
			expr[head].addTerm(-1, w[i]);
			//System.out.println("FLOW FROM "+DD.nodes.get(tail).key+" to "+DD.nodes.get(head).key+"        "+tail+" to "+head+" X"+DD.arcs.get(i).varIndex+" val: "+DD.arcs.get(i).arcVal);
		}
		cplex.addEq(expr[0], 1, "c0");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){cplex.addEq(expr[n], 0, "flow_"+0+"_"+n);}
		}

		// Link leader variables
		IloLinearNumExpr[] exprLink = new IloLinearNumExpr[data.P.nVarL];
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i] = cplex.linearNumExpr();
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {	
				exprLink[DD.arcs.get(i).varIndex].addTerm(1, w[i]);
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i].addTerm(-1, x[i]);
			cplex.addEq(exprLink[i], 0);
		}

		///////////////////////////////////////////////////// VALID INEQUALITIES ///////////////////////////////////////////
		IloLinearNumExpr exprQ = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprQ.addTerm(y[q], data.P.cfol[q]);
		}
		cplex.addLe(exprQ, folUB);

		/////////////////////////////////////////////////// DD VALUE FUNCTION CONSTRAINT/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			exprV.addTerm(w[i], -DD.arcs.get(i).cost);
		}
		cplex.addLe(exprV, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}
		cplex.addMinimize(obj,"OBJ");
		cplex.solve();
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		DD-IP-R is infeasible!");
			LB = UB; //If UB= 99999 then problem is infeasible, otherwise cplex just closed the gap!
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			int LB1 = (int) Math.round(cplex.getBestObjValue());
			double time1 = (System.currentTimeMillis()-Atime)/1000; 
			System.out.println("******************************************************  MIP ROOT LB: "+ LB+" UB: "+UB);

			//// NOW ADD THE VIs and run again

			// block variables for the new VIS
			IloNumVar[][] r = new IloNumVar[data.P.nConsF][data.P.nVarF]; 
			for (int i = 0; i < data.P.nVarF; i++) {
				for(int j = 0; j < data.P.nConsF; j++)
					if(data.cInteraction[j].equals("F") == false) {//block directly with x_i!
						r[j][i] = cplex.numVar(0, 0, IloNumVarType.Float, "r");
					}
					else {//Here we can block
						r[j][i] = cplex.numVar(0, 1, IloNumVarType.Float, "r");
					}
			}
			// New VIs: each follower variable with negative obj coefficient set to zero only if blocked 
			for (int i = 0; i <data.P.nVarF; i++) {
				if(data.P.cfol[i] < 0 ) {
					IloLinearNumExpr exprf = cplex.linearNumExpr();
					exprf.addTerm(1, y[i]);
					exprf.addTerm(1, x[i]);	//Block it directly
					for(int j=0; j < data.P.nConsF; j++){//Block it by some follower constraint
						if(data.cInteraction[j].equals("F")) {exprf.addTerm(1, r[j][i]);}	
					}
					cplex.addGe(exprf, 1);
					// Now define the blocking variables
					for(int j=0; j < data.P.nConsF; j++){
						double block = data.P.bfol[j]-data.P.Afol[j][i]+1;
						IloLinearNumExpr exprB = cplex.linearNumExpr();
						for (int q = 0; q <data.P.nVarF; q++) {
							if(q!=i) {exprB.addTerm(data.P.Afol[j][q], y[q]);}
						}
						for (int q = 0; q <data.P.nVarL; q++) {
							exprB.addTerm(data.P.Bfol[j][q], x[q]);
						}
						exprB.addTerm(-block, r[j][i]);
						//System.out.println("CONSTRAINT "+j+" block: "+block);
						cplex.addGe(exprB, 0);	
					}
				}
			}
			cplex.solve();
			int LB2 = (int) Math.round(cplex.getBestObjValue());
			System.out.println("RESULT;"+file+";"+LB1+";"+time1+";"+LB2+";"+(System.currentTimeMillis()-Atime)/1000);

		}	
	}

	// Enhanced with the blocking constraints, no callbacks
	public void singleLevelDDRelaxedNCRoot(DataHandler data, int timeLim, double Atime, String file) throws IloException {
		// Don't add any cuts!!!
		sampleMaster.clear();
		sampleMap.clear();

		// Solve the follower problem for each node/arc in the last layer (not solved before)
		for (int i = DD.arcs.size()-1; i >=0 ; i--) {
			if(DD.arcs.get(i).varIndex==-1) { //Arcs in the last layer
				int nodeIndex = DD.arcs.get(i).tail;
				if(DD.nodes.get(nodeIndex).cutGenerated == false) {
					//System.out.println("NODE OF INTEREST "+DD.nodes.get(nodeIndex).key);
					int folObjValR = solveFollower(data, DD.nodes.get(nodeIndex).S); 
					DD.arcs.get(i).cost = folObjValR;
					DD.nodes.get(nodeIndex).cutGenerated = true;
					DD.nodes.get(nodeIndex).folResponse = folObjValR;
					if (folObjValR > folUB) {folUB = folObjValR;}
					//System.out.println("FOLLOWER OBJ RELAXED STATE: "+folObjValR);
				}
				else {
					DD.arcs.get(i).cost = Math.min(DD.nodes.get(nodeIndex).folResponse, folObjM);
				}
			}
		}
		// Propagate costs backwards
		for (int i = DD.nodes.size()-1; i >= 2 ; i--) { // 0 and 1 are root node and end node
			DD.nodes.get(i).propagateCost(data);
		}
		System.out.println(" FOLL OBJ UB: "+folUB);

		int timeUsed = (int) ((System.currentTimeMillis()-Atime)/1000);
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		double lim = Math.max( 10, timeLim-timeUsed); //For safety!
		cplex.setParam(IloCplex.Param.TimeLimit, lim);


		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}

		// Flow variables for the BDDs
		IloNumVar[] w = new IloNumVar[DD.arcs.size()]; 
		for (int i = 0; i < DD.arcs.size(); i++) {
			w[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "w");
		}

		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}
		// New VIs: each follower binary variable with negative obj coefficient set to zero only if blocked 
		int[] varRecord = new int[data.P.nVarF];  //Keep track of which variables are already represented in the Sample
		for (int q = 0; q < sampleMaster.size(); q++) {
			for(int i=0; i < data.P.nVarF; i++){
				if(sampleMaster.get(q).y[i] > 0.9) {
					varRecord[i]++;
				}
			}
		}

		// FLOW CONSTRAINTS	
		// Add flow constraints (but drop the last one)
		IloLinearNumExpr[] expr = new IloLinearNumExpr[DD.nodes.size()];
		for (int n = 0; n < DD.nodes.size(); n++) {
			expr[n] = cplex.linearNumExpr();
		}	
		for (int i = 0; i < DD.arcs.size(); i++) {
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			expr[tail].addTerm(1, w[i]);
			expr[head].addTerm(-1, w[i]);
			//System.out.println("FLOW FROM "+DD.nodes.get(tail).key+" to "+DD.nodes.get(head).key+"        "+tail+" to "+head+" X"+DD.arcs.get(i).varIndex+" val: "+DD.arcs.get(i).arcVal);
		}
		cplex.addEq(expr[0], 1, "c0");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){cplex.addEq(expr[n], 0, "flow_"+0+"_"+n);}
		}

		// Link leader variables
		IloLinearNumExpr[] exprLink = new IloLinearNumExpr[data.P.nVarL];
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i] = cplex.linearNumExpr();
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {	
				exprLink[DD.arcs.get(i).varIndex].addTerm(1, w[i]);
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			exprLink[i].addTerm(-1, x[i]);
			cplex.addEq(exprLink[i], 0);
		}

		///////////////////////////////////////////////////// VALID INEQUALITIES ///////////////////////////////////////////
		IloLinearNumExpr exprQ = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprQ.addTerm(y[q], data.P.cfol[q]);
		}
		cplex.addLe(exprQ, folUB);

		/////////////////////////////////////////////////// DD VALUE FUNCTION CONSTRAINT/////////////////////////////////////////////////////////
		IloLinearNumExpr exprV = cplex.linearNumExpr();
		for(int q=0; q<data.P.nVarF; q++){
			exprV.addTerm(y[q], data.P.cfol[q]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			exprV.addTerm(w[i], -DD.arcs.get(i).cost);
		}
		cplex.addLe(exprV, 0);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}
		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		cplex.solve();
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		DD-IP-R is infeasible!");
			LB = UB; //If UB= 99999 then problem is infeasible, otherwise cplex just closed the gap!
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			int LB1 = (int)Math.round(cplex.getBestObjValue());
			double time1 = (System.currentTimeMillis()-Atime)/1000;
			System.out.println("******************************************************  MIP ROOT LB: "+ LB+" UB: "+UB);

			// Add the VIs and run again

			// block variables for the new VIS
			IloNumVar[][] r = new IloNumVar[data.P.nConsF][data.P.nVarF]; 
			for (int i = 0; i < data.P.nVarF; i++) {
				for(int j = 0; j < data.P.nConsF; j++)
					if(data.P.Afol[j][i] <= 0 ) {// Impossible to block!
						r[j][i] = cplex.numVar(0, 0, IloNumVarType.Float, "r");
					}
					else {//Here we can block
						r[j][i] = cplex.numVar(0, 1, IloNumVarType.Bool, "r");
					}
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				if(data.P.cfol[i] < 0) {//Only add for variables likely to be part of opt response
					IloLinearNumExpr exprf = cplex.linearNumExpr();
					exprf.addTerm(1, y[i]);
					for(int j=0; j < data.P.nConsF; j++){
						exprf.addTerm(1, r[j][i]);	//Block it by some constraint
					}
					cplex.addGe(exprf, 1);
					// Now define the blocking variables
					for(int j=0; j < data.P.nConsF; j++){
						int block = (int) (data.P.bfol[j]-data.P.Afol[j][i]+1);
						IloLinearNumExpr exprB = cplex.linearNumExpr();
						for (int q = 0; q <data.P.nVarF; q++) {
							if(q!=i) {exprB.addTerm(data.P.Afol[j][q], y[q]);}
						}
						for (int q = 0; q <data.P.nVarL; q++) {
							exprB.addTerm(data.P.Bfol[j][q], x[q]);
						}
						exprB.addTerm(-block-999999999, r[j][i]);
						//System.out.println("CONSTRAINT "+j+" block: "+block);
						cplex.addGe(exprB, -999999999);	
					}
				}
			}
			cplex.solve();
			int LB2 = (int)Math.round(cplex.getBestObjValue());
			System.out.println("RESULT;"+file+";"+LB1+";"+time1+";"+LB2+";"+(System.currentTimeMillis()-Atime)/1000);

		}	
	}

	/**
	 * Computes summary statistics of the linking constraints for reporting.
	 */
	public void computeExtras(DataHandler data) {
		double[] ratios = new double[data.P.nConsF];
		double min = 999999999;
		double max = -999999999;
		double total = 0;
		int numLinking = 0;

		for (int i = 0; i < data.P.nConsF; i++) {
			if(data.cInteraction[i].equals("B")) {
				double sum = 0;
				for (int j = 0; j < data.P.nVarF; j++) {
					sum+=data.P.Afol[i][j];
				}
				ratios[i] = data.P.bfol[i]/(sum+0.0);
				total+=ratios[i];
				if(ratios[i] < min) {
					min = ratios[i];
				}
				if(ratios[i] > max) {
					max = ratios[i];
				}
				numLinking++;
			}
		}
		System.out.println(min+";"+(total/(numLinking+0.0))+";"+max+";"+numLinking);

	}

	// Solves high point relaxation
	public void solveHPointCuts(DataHandler data, int timeLim) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.Threads, numCores);
		cplex.setParam(IloCplex.Param.TimeLimit, timeLim);

		//Leader and follower variables 
		IloNumVar[] x = new IloNumVar[data.P.nVarL];
		IloNumVar[] y = new IloNumVar[data.P.nVarF];
		for (int i = 0; i < data.P.nVarL; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x["+i+"]");
		}
		for (int i = 0; i < data.P.nVarF; i++) {
			y[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "y["+i+"]");
		}
		// block variables for the new VIS
		IloNumVar[][] r = new IloNumVar[data.P.nConsF][data.P.nVarF]; 
		for (int i = 0; i < data.P.nVarF; i++) {
			for(int j = 0; j < data.P.nConsF; j++)
				if(data.P.Afol[j][i] <= 0 ) {// Impossible to block!
					r[j][i] = cplex.numVar(0, 0, IloNumVarType.Float, "r");
				}
				else {//Here we can block
					r[j][i] = cplex.numVar(0, 1, IloNumVarType.Bool, "r");
				}
		}
		for (int i = 0; i <data.P.nVarF; i++) {
			if(data.P.cfol[i] < 0) {
				IloLinearNumExpr exprf = cplex.linearNumExpr();
				exprf.addTerm(1, y[i]);
				for(int j=0; j < data.P.nConsF; j++){
					exprf.addTerm(1, r[j][i]);	//Block it by some constraint
				}
				cplex.addGe(exprf, 1);
				// Now define the blocking variables
				for(int j=0; j < data.P.nConsF; j++){
					int block = (int) (data.P.bfol[j]-data.P.Afol[j][i]+1);
					IloLinearNumExpr exprB = cplex.linearNumExpr();
					for (int q = 0; q <data.P.nVarF; q++) {
						if(q!=i) {exprB.addTerm(data.P.Afol[j][q], y[q]);}
					}
					for (int q = 0; q <data.P.nVarL; q++) {
						exprB.addTerm(data.P.Bfol[j][q], x[q]);
					}
					exprB.addTerm(-block-999999999, r[j][i]);
					//System.out.println("CONSTRAINT "+j+" block: "+block);
					cplex.addGe(exprB, -999999999);	
				}
			}
		}

		// add constraints leader
		for(int j=0; j < data.P.nConsL; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Aleader[j][i], x[i]);
			}
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Bleader[j][i], y[i]);
			}
			cplex.addLe(exprf, data.P.b[j]);
		}
		// add constraints follower
		for(int j=0; j < data.P.nConsF; j++){
			IloLinearNumExpr exprf = cplex.linearNumExpr();
			for (int i = 0; i <data.P.nVarF; i++) {
				exprf.addTerm(data.P.Afol[j][i], y[i]);
			}
			for (int i = 0; i <data.P.nVarL; i++) {
				exprf.addTerm(data.P.Bfol[j][i], x[i]);
			}
			cplex.addLe(exprf, data.P.bfol[j]);
		}


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for(int i=0; i<data.P.nVarL; i++){
			obj.addTerm(x[i], data.P.c1leader[i]);
		}
		for(int i=0; i<data.P.nVarF; i++){
			obj.addTerm(y[i], data.P.c2leader[i]);
		}

		cplex.addMinimize(obj,"OBJ");

		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		HPR is infeasible!");
			LB = 999999999; // problem is infeasible
		}
		else if(cplex.getStatus() == IloCplex.Status.Optimal) {
			LB = (int) Math.round(cplex.getObjValue());
			System.out.println("***********************************  High point zstar: "+ LB);
		}
	}

	
	void exportInstance(DataHandler data, String name) throws FileNotFoundException {
//		//// CREATE MPS FILE
//		String name=file.substring(0, file.length()-3);
//		name+="plusDD";
		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream(name+".mps", false));
		ps.println("NAME"+" "+name);
		
		
		/////////////////////////////////////////// CONSTRAINTS ///////////////////////////////////////////////////////////////////
		ps.println("ROWS");
		ps.println(" "+"N"+" "+"COST");
		for (int i = 0; i < data.P.nConsL; i++) {
			ps.println(" "+"L"+" "+"LK"+i);
		}
		for (int i = 0; i < data.P.nConsF; i++) {
			ps.println(" "+"L"+" "+"FK"+i);
		}
		for (int n = 0; n < DD.nodes.size(); n++) {
			ps.println(" "+"E"+" "+"FLOW"+n);
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			ps.println(" "+"E"+" "+"LINK"+i);
		}
		ps.println(" "+"L"+" "+"VALFUNCTION");
		ps.println(" "+"L"+" "+"VI");

		
		/////////////////////////////////////////// VARIABLES ///////////////////////////////////////////////////////////////////
		ps.println("COLUMNS");
		// Leader and follower
		ps.println(" MARK0000 'MARKER' 'INTORG'");
		for (int i = 0; i <data.P.nVarL; i++) {
			ps.println(" "+data.P.leaderVarNames[i]+" "+"COST"+" "+(data.P.c1leader[i]));
			for (int j = 0; j < data.P.nConsL; j++) {
				ps.println(" "+data.P.leaderVarNames[i]+" "+"LK"+j+" "+data.P.Aleader[j][i]);
			}
			for (int j = 0; j < data.P.nConsF; j++) {
				ps.println(" "+data.P.leaderVarNames[i]+" "+"FK"+j+" "+data.P.Bfol[j][i]);
			}
			ps.println(" "+data.P.leaderVarNames[i]+" "+"LINK"+i+" "+"-1");
		}
		for (int i = 0; i <data.P.nVarF; i++) {
			ps.println(" "+data.P.followerVarNames[i]+" "+"COST"+" "+(data.P.c2leader[i]));
			for (int j = 0; j < data.P.nConsL; j++) {
				ps.println(" "+data.P.followerVarNames[i]+" "+"LK"+j+" "+data.P.Bleader[j][i]);
			}
			for (int j = 0; j < data.P.nConsF; j++) {
				ps.println(" "+data.P.followerVarNames[i]+" "+"FK"+j+" "+data.P.Afol[j][i]);
			}
			ps.println(" "+data.P.followerVarNames[i]+" "+"VALFUNCTION"+" "+data.P.cfol[i]);
			ps.println(" "+data.P.followerVarNames[i]+" "+"VI"+" "+data.P.cfol[i]);
		}
		ps.println(" MARK0001 'MARKER' 'INTEND'");
		// Flow variables
		for (int i = 0; i < DD.arcs.size(); i++) {
			ps.println(" "+"W"+i+" "+"COST"+" "+"0");
			int tail = DD.arcs.get(i).tail;
			int head = DD.arcs.get(i).head;
			ps.println(" "+"W"+i+" "+"FLOW"+tail+" "+"1");
			ps.println(" "+"W"+i+" "+"FLOW"+head+" "+"-1");
			/// VALUE FUNCTION
			double coeff = -(DD.arcs.get(i).cost+0.0000001); //Tolerance for numeric problems
			ps.println(" "+"W"+i+" "+"VALFUNCTION"+" "+String.format("%.7f", coeff)+" ");
			// LINKING CONSTRAINTS
			if(DD.arcs.get(i).leaderVar==true && DD.arcs.get(i).arcVal==1) {
				ps.println(" "+"W"+i+" "+"LINK"+DD.arcs.get(i).varIndex+" "+"1");
			}
		}

		
		/////////////////////////////////////////// RHS ///////////////////////////////////////////////////////////////////
		ps.println("RHS");
		for (int j = 0; j < data.P.nConsL; j++) {
			ps.println(" "+"rhs"+" "+"LK"+j+" "+data.P.b[j]);
		}
		for (int j = 0; j < data.P.nConsF; j++) {
			ps.println(" "+"rhs"+" "+"FK"+j+" "+data.P.bfol[j]);
		}
		ps.println(" "+"rhs"+" "+"FLOW"+"0"+" "+"1");
		ps.println(" "+"rhs"+" "+"FLOW"+DD.lastNode+" "+"-1");
		for (int n = 1; n < DD.nodes.size(); n++) {
			if(n!=DD.lastNode){
				ps.println(" "+"rhs"+" "+"FLOW"+n+" "+"0");	
			}
		}
		for (int i = 0; i < data.P.nVarL; i++) {
			ps.println(" "+"rhs"+" "+"LINK"+i+" "+"0");
		}
		ps.println(" "+"rhs"+" "+"VALFUNCTION"+" "+"0");
		ps.println(" "+"rhs"+" "+"VI"+" "+folUB);

		
		/////////////////////////////////////////// BOUNDS ///////////////////////////////////////////////////////////////////
		ps.println("BOUNDS");
		for (int i = 0; i <data.P.nVarL; i++) {
			ps.println(" "+"BV bnd"+"	"+data.P.leaderVarNames[i]);
		}
		for (int i = 0; i <data.P.nVarF; i++) {
			ps.println(" "+"BV bnd"+"	"+data.P.followerVarNames[i]);
		}
		for (int i = 0; i < DD.arcs.size(); i++) {
			ps.println(" "+"LO BND1"+"	"+"W"+i + " 0");
			ps.println(" "+"UP BND1"+"	"+"W"+i + " 1");
		}
		ps.println("ENDATA");

		//NOW THE AUX FILE ///////////////////////////////////////////////////////////////////////
		java.io.PrintStream aux = new java.io.PrintStream( new java.io.FileOutputStream(name+".aux", false));
		aux.println("@NUMVARS");
		aux.println(data.P.nVarF);
		aux.println("@NUMCONSTRS");
		aux.println(data.P.nConsF);
		aux.println("@VARSBEGIN");
		for (int i = 0; i <data.P.nVarF; i++) {
			aux.println(data.P.followerVarNames[i]+" "+String.format("%.0f",data.P.cfol[i])+".");
		}
		aux.println("@VARSEND");
		aux.println("@CONSTRSBEGIN");
		for (int i = 0; i < data.P.nConsF; i++) {
			aux.println("FK"+i);
		}
		aux.println("@CONSTRSEND");
		aux.println("@NAME");
		aux.println(name);
		aux.println("@MPS");
		aux.println(name+".mps");
		System.out.println("DONE! "+name);
	}

}
