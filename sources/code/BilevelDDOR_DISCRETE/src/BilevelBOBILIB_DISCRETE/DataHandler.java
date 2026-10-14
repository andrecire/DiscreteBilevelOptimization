/**
 * Stores a parsed discrete bilevel instance together with solver-specific
 * metadata derived from that instance.
 *
 * <p>In addition to the raw model coefficients, this class identifies the
 * follower-constraint interaction pattern and detects whether the instance fits
 * the specialized interdiction structure exploited by the solver.</p>
 */
package BilevelBOBILIB_DISCRETE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;
import java.util.StringTokenizer;

import ilog.concert.IloException;
import ilog.concert.IloLPMatrix;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloLinearNumExprIterator;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class DataHandler {

	// Instance data
	Instance P;
	// Bounds to prune the dd
	int[] stateBounds;
	boolean isInterdiction;	 // True if all the interaction constraints are of the form x_i+y_i <= 1	
	String [] cInteraction;  // F: only follower variables, L: only leader variables, B: both follower and leader variables

	// Some performance measures
	int iniNodes;
	int iniArcs; 
	double density;

	// Empty constructor; the instance is loaded later through readInstance.
	public DataHandler() {

	}

	/**
	 * Computes the decision-diagram width to use for the current run.
	 *
	 * @param setWidth user-provided width; positive values override the heuristic
	 * @return diagram width used by the solver
	 */
	public int calcWidth(int setWidth) {
		if(setWidth > 0) {
			return setWidth;
		}
		else {
			if(P.nVarL < 100) {
				return 100;
			}
			else if(P.nVarL <= 500) {
				return 50;
			}
			else if (P.nVarL <= 1000) {
				return 25;
			}

			return 10;
		}
	}


	/**
	 * Reads a discrete .dd instance and derives the auxiliary structures used by
	 * the solver.
	 *
	 * @param fileName path to the instance file
	 * @throws IOException if the file cannot be read
	 */
	public void readInstance(String fileName) throws IOException {

		P = Instance.readFromFile(fileName);
		stateBounds = new int[P.nConsF];
		P.print();

		// Fill out interactions
		cInteraction = new String[P.nConsF];
		for (int j = 0; j < P.nConsF; j++) {
			int folNZ = 0;			//Fol nonzero coeff
			int leaderNZ = 0;		// Leader nonzero coeff	
			for (int i = 0; i < P.nVarL; i++) {
				if(P.Bfol[j][i]!=0) {leaderNZ++;}
			}
			for (int i = 0; i < P.nVarF; i++) {
				if(P.Afol[j][i]!=0) {folNZ++;}
			}
			if(folNZ==0) {cInteraction[j] = "L";} //Leader only
			else if(leaderNZ==0) {cInteraction[j] ="F";} //Follower only
			else {cInteraction[j] = "B";} //Both
		}

		// Check if problem is "interdiction"
		isInterdiction = true;
		//Check that all variables are binary
		for (int i = 0; i < P.nVarL; i++) {
			if(!P.typeLeaderStr[i].equals("bin")) {isInterdiction = false;}
		}
		//Check that all follower variables are binary or nonnegative integer 
		for (int i = 0; i < P.nVarF; i++) {
			if(P.typeFollowerStr[i].equals("cont")) {isInterdiction = false;}
			if(P.lbFollower[i] < 0) {isInterdiction = false;}
		}
		//Check that for all follower constraints the coefficients are nonzero for follower variables
		for (int i = 0; i < P.nVarF; i++) {
			for (int j = 0; j < P.nConsF; j++) { 
				if(P.Afol[j][i] < 0) {
					isInterdiction = false;
				}
			}
		}
		// Check that all linking follower constraints are of the form x_i+y_i <= 1.
		if(P.nVarL != P.nVarF) {
			isInterdiction = false;
		}
		else {//Same number of variables!!
			for (int i = 0; i < P.nVarL; i++) {
				int numCleader = 0;
				int numCfol = 0;
				int indexLeader = -1;
				int indexFollower = -1;
				for (int j = 0; j < P.nConsF; j++) {
					if(cInteraction[j].equals("B") && P.Bfol[j][i]!=0 ) {
						if(P.Bfol[j][i]!=1) {
							isInterdiction = false;
						}
						else {
							numCleader++;
							indexLeader = j;
						}
					}
					if(cInteraction[j].equals("B") && P.Afol[j][i]!=0 ) {
						if(P.Afol[j][i]!=1) {
							isInterdiction = false;
						}
						else {
							numCfol++;
							indexFollower = j;
						}
					}
				}
				// Check number of coef equal to 1
				if(numCleader != 1 || numCfol!=1) {
					isInterdiction = false;
				}
				//Check the variable x_i is in the same constraint as y_i
				if(indexLeader != indexFollower) {
					isInterdiction = false;
				}

			}
		}

		System.out.println("DONE READING THE INSTANCE. IS IT AN INTERDICTION INSTANCE? "+isInterdiction);


	}



}

