/**
 * Entry point for the mixed DD-based solver.
 *
 * <p>Key assumption for this variant: all leader variables are binary.</p>
 */
package BilevelBOBILIB_MIXED;

import ilog.concert.IloException;

import java.awt.DefaultFocusTraversalPolicy;
import java.io.File;
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Random;


public class Main {

	// public static void main(String[] args) throws IOException, InterruptedException, IloException {

	// 	java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("Stats.txt", true));

	// 	// Read a file with the input data
	// 	String fileName = "dataBOBILIB/ran14x18-disj-8_0_100.dd";
	// 	int timeLimit = 3600;
	// 	int numCores = 1;				
	// 	String outputFolder = ".";
	// 	int algSwitch = 0;				//0: DD approach, 1: DD only root IP no cuts, 2: Only HPR IP, 3: HPR+Logical Cuts
	// 	int setWidth = 0;				//0: calc width default. Any positive value: use as width
	// 	Boolean logicalCuts = true;		//Turn on and off the logical cuts. Only works with algSwith=0. 
		
	// 	File file = new File(fileName);
	// 	String instanceName = file.getName();
	// 	System.out.println(fileName); 
	// 	java.io.PrintStream psSol = new java.io.PrintStream( new java.io.FileOutputStream(outputFolder+"/dd_"+instanceName+"-res.json"));	

	// 	DataHandler data = new DataHandler();
	// 	data.readInstance(fileName);
	// 	//data.outputInstance();

	// 	// Begin the time count						
	// 	double Atime = System.currentTimeMillis();

	// 	runSingleLevelDD(data, Atime, timeLimit, numCores, ps, psSol, instanceName, algSwitch,setWidth, logicalCuts);

	// 	data = null;
	// 	System.gc();
	// }
	/**
	 * Parses command-line arguments, reads the input instance, and launches the
	 * selected algorithmic workflow.
	 */
	public static void main(String[] args) {
        try {
            // Check that the required arguments are provided
            if (args.length < 7) {
                System.err.println("Usage: java dd_solver_mixed.jar <fileName> <timeLimit> <numCores> <outputFolder> <algSwitch> <setWidth> <logicalCuts>");
                System.exit(1);
            }

            // Parse command-line arguments
            String fileName     = args[0];
            int timeLimit       = Integer.parseInt(args[1]);
            int numCores        = Integer.parseInt(args[2]);
            String outputFolder = args[3];
            int algSwitch       = Integer.parseInt(args[4]);
            int setWidth        = Integer.parseInt(args[5]);
            boolean logicalCuts = Boolean.parseBoolean(args[6]);

            // Friendly printing of arguments
            System.out.println("Command-line arguments received:");
            System.out.println("\tfileName    = " + fileName);
            System.out.println("\ttimeLimit   = " + timeLimit);
            System.out.println("\tnumCores    = " + numCores);
            System.out.println("\toutputFolder= " + outputFolder);
            System.out.println("\talgSwitch   = " + algSwitch);
            System.out.println("\tsetWidth    = " + setWidth);
            System.out.println("\tlogicalCuts = " + logicalCuts);

            // Prepare output
            File file = new File(fileName);
            String instanceName = file.getName();
            String solutionFilename = outputFolder + "/dd_" + instanceName + "_" + setWidth + "-res.json";

            // Read input instance
            DataHandler data = new DataHandler();
            data.readInstance(fileName);

            // Begin the time count
            double Atime = System.currentTimeMillis();

            // Run algorithm
            runSingleLevelDD(data, Atime, timeLimit, numCores, instanceName, algSwitch, setWidth, logicalCuts, solutionFilename);

            // Cleanup
            data = null;
            System.gc();

            System.exit(0); // success

        } catch (Exception e) {
            System.err.println("An error occurred while running the program:");
            e.printStackTrace();
            System.exit(1); // failure
        }
    }


	/**
	 * Runs the requested mixed DD workflow and exports the solution when the
	 * exact method is selected.
	 */
	private static void runSingleLevelDD(DataHandler data, double Atime, int timeLimit, int numCores, String file, int algSwitch, int setWidth, Boolean logicalCuts, String solutionFilename) throws InterruptedException, IloException {
		// Create an AlgorithmHandler
		AlgorithmHandler alg = new AlgorithmHandler(data, numCores, logicalCuts);

		// Run with Relaxed DD
		boolean solved = alg.getM(data, timeLimit);
		if(!solved) {System.out.println("COULD NOT EVEN SOLVE THE HPP WITHIN THE TIME LIMIT!!!!!!!!!!");}
		else {
			if(algSwitch <= 1) {
				alg.refineMS(data, 30);
				alg.genRelaxedDD(data, setWidth);
				alg.reduceDD(data, Atime);
				if(data.isInterdiction) {
					if(algSwitch == 0) {alg.singleLevelDDInterdiction(data, timeLimit, Atime);}
					if(algSwitch == 1) {alg.singleLevelDDInterdictionRoot(data, timeLimit, Atime, file);}
				}
				else {
					if(algSwitch == 0) {alg.singleLevelDDRelaxedNC(data, timeLimit, Atime);}
					if(algSwitch == 1) {alg.singleLevelDDRelaxedNCRoot(data, timeLimit, Atime, file);}					
				}
				System.out.println("DDMAXMIN;"+file+";"+data.P.nVarL+";"+data.P.nVarF+";"+data.P.nConsL+";"+data.P.nConsF+";"+(System.currentTimeMillis()-Atime)/1000+";"+alg.timeDD+";"+alg.timeMaxMin+";"+alg.LB+";"+alg.UB+";"+alg.DD.nodes.size()+";"+alg.DD.arcs.size()+";"+alg.numCuts+";"+alg.isExact);
				if(algSwitch == 0) {
					System.out.println("Best x: "+Arrays.toString(alg.xstar));
					System.out.println("Best y: "+Arrays.toString(alg.ystar));
					alg.checkStarSol(data);
					System.out.println("RESULT;"+file+";"+alg.UB+";"+(System.currentTimeMillis()-Atime)/1000+";"+(alg.UB-alg.LB)/(alg.UB+1e-12));

					// Create the JSON solution file.
					System.out.println("Exporting solution...");
					try {
 			            PrintStream psSol = new PrintStream(
			                	new FileOutputStream(solutionFilename)
			            );	
						alg.ExportSol(data, file, psSol);
					} catch (Exception e) {
			            System.err.println("Error exporting solution");
			            e.printStackTrace();
					}					
				}
			}
			else if(algSwitch == 2) {//Run HPR
				alg.solveHPoint(data, timeLimit);
				System.out.print("RESULT;"+file+";"+alg.LB+";"+(System.currentTimeMillis()-Atime)/1000+";EXTRAS;");
				alg.computeExtras(data);
			}
			else if(algSwitch == 3) {//Run HPR with logical cuts
				alg.solveHPointCuts(data, timeLimit);
				System.out.print("RESULT;"+file+";"+alg.LB+";"+(System.currentTimeMillis()-Atime)/1000+";EXTRAS;");
				alg.computeExtras(data);
			}
		}

		

	}
}
