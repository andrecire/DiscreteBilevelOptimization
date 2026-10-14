package BilevelBOBILIB_DISCRETE;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Represents one discrete bilevel optimization instance stored in the custom
 * plain-text `.dd` format.
 * 
 * Structure:
 * - Leader objective: min c1^T x + c2^T y
 * - Follower objective: min cfol^T y
 * - Leader constraints: Aleader x + Bleader y <= b
 * - Follower constraints: Afol y + Bfol x <= bfol
 * 
 * Notes:
 * - Problems are always minimization problems
 * - Inequalities are always in the form <=
 * - Variable bounds are included as inequalities in the constraints
 */

public class Instance {

	// Problem dimensions
    public int nVarL, nVarF;   // Number of leader and follower variables
    public int nConsL, nConsF; // Number of leader and follower constraints

    // Constraint matrices
    public double[][] Aleader; // Leader matrix on leader vars (A)
    public double[][] Bleader; // Leader matrix on follower vars (B)

    public double[][] Afol;    // Follower matrix on follower vars (A)
    public double[][] Bfol;    // Follower matrix on leader vars (B)

    // Objective vectors
    public double[] c1leader;  // Coefficients on x in leader objective
    public double[] c2leader;  // Coefficients on y in leader objective
    public double[] cfol;      // Coefficients on y in follower objective

    // Right-hand sides
    public double[] b;     // RHS for leader constraints
    public double[] bfol;  // RHS for follower constraints

    // Variable bounds
    // If + or - infinity, use Double.POSITIVE_INFINITY or Double.NEGATIVE_INFINITY.
    // You can also use isInfinite() implemented here to check for infinity.
    public double[] lbLeader;
    public double[] ubLeader;
    public double[] lbFollower;
    public double[] ubFollower;

    // String representations of the types
    // Types can be 'bin', 'int', and 'cont' (continuous)
    public String[] typeLeaderStr;      // Type of leader variables
    public String[] typeFollowerStr;    // Type of follower variables

    // Variable names in the original MPS file
    public String[] leaderVarNames;
    public String[] followerVarNames;
 
    /**
     * Reads a `.dd` file and returns the corresponding instance object.
     *
     * @param filename path to the .dd file
     * @return populated instance
     * @throws IOException if reading fails
     */
    public static Instance readFromFile(String filename) throws IOException {
    	Instance d = new Instance();

        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            // First line: dimensions
            String[] header = reader.readLine().trim().split(" ");
            d.nVarL = Integer.parseInt(header[0]);
            d.nVarF = Integer.parseInt(header[1]);
            d.nConsL = Integer.parseInt(header[2]);
            d.nConsF = Integer.parseInt(header[3]);

            // Read constraint matrices
            d.Aleader = readMatrix(reader, d.nConsL, d.nVarL);
            d.Bleader = readMatrix(reader, d.nConsL, d.nVarF);
            d.Bfol    = readMatrix(reader, d.nConsF, d.nVarL);
            d.Afol    = readMatrix(reader, d.nConsF, d.nVarF);

            // Read objective coefficients
            d.c1leader = readVector(reader, d.nVarL);
            d.c2leader = readVector(reader, d.nVarF);
            d.cfol     = readVector(reader, d.nVarF);

            // Read RHS vectors
            d.b        = readVector(reader, d.nConsL);
            d.bfol     = readVector(reader, d.nConsF);
            
            // Read variable bounds for leader
            d.lbLeader = new double[d.nVarL];
            d.ubLeader = new double[d.nVarL];
            String[] boundsL = reader.readLine().trim().split(" ");
            for (int i = 0; i < d.nVarL; i++) {
                d.lbLeader[i] = parseBound(boundsL[2 * i]);
                d.ubLeader[i] = parseBound(boundsL[2 * i + 1]);
            }

            // Read variable bounds for follower
            d.lbFollower = new double[d.nVarF];
            d.ubFollower = new double[d.nVarF];
            String[] boundsF = reader.readLine().trim().split(" ");
            for (int i = 0; i < d.nVarF; i++) {
                d.lbFollower[i] = parseBound(boundsF[2 * i]);
                d.ubFollower[i] = parseBound(boundsF[2 * i + 1]);
            }

            // Read variable types
            d.typeLeaderStr = readStringVector(reader, d.nVarL);
            d.typeFollowerStr = readStringVector(reader, d.nVarF);

            // Read original variable names
            d.leaderVarNames = readStringVector(reader, d.nVarL);
            d.followerVarNames = readStringVector(reader, d.nVarF);
        }

        return d;
    }

 
    /**
     * Reads a matrix from the file.
     *
     * @param reader BufferedReader positioned at the matrix start
     * @param rows number of rows
     * @param cols number of columns
     * @return parsed matrix
     * @throws IOException if reading fails
     */
    private static double[][] readMatrix(BufferedReader reader, int rows, int cols) throws IOException {
        double[][] mat = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            String[] line = reader.readLine().trim().split(" ");
            for (int j = 0; j < cols; j++) {
                mat[i][j] = Double.parseDouble(line[j]);
            }
        }
        return mat;
    }

    
    /**
     * Reads a vector from the file.
     *
     * @param reader BufferedReader positioned at the vector start
     * @param len length of the vector
     * @return parsed vector
     * @throws IOException if reading fails
     */
    private static double[] readVector(BufferedReader reader, int len) throws IOException {
        String[] line = reader.readLine().trim().split(" ");
        double[] vec = new double[len];
        for (int i = 0; i < len; i++) {
            vec[i] = Double.parseDouble(line[i]);
        }
        return vec;
    }

    
    /**
     * Prints the raw matrix and vector data of the instance.
     */
    public void printMatrix() {
        System.out.println("=== INSTANCE SUMMARY ===");
        System.out.printf("Leader vars: %d, Follower vars: %d\n", nVarL, nVarF);
        System.out.printf("Leader cons: %d, Follower cons: %d\n", nConsL, nConsF);

        System.out.println("\nLEADER OBJECTIVE: min c1^T x + c2^T y");
        System.out.print("c1leader: ");
        printVector(c1leader);
        System.out.print("c2leader: ");
        printVector(c2leader);

        System.out.println("\nFOLLOWER OBJECTIVE: min cfol^T y");
        System.out.print("cfol: ");
        printVector(cfol);

        System.out.println("\nAleader (A x in leader):");
        printMatrix(Aleader);
        System.out.println("Bleader (B y in leader):");
        printMatrix(Bleader);
        System.out.print("b (leader RHS): ");
        printVector(b);

        System.out.println("\nAfol (B x in follower):");
        printMatrix(Bfol);
        System.out.println("Afol (A y in follower):");
        printMatrix(Afol);
        System.out.print("bfol (follower RHS): ");
        printVector(bfol);
    }

    /**
     * Prints the instance summary.
     */
    public void print() {
        System.out.println("=== INSTANCE SUMMARY ===");
        System.out.printf("Leader variables: %d, Follower variables: %d\n", nVarL, nVarF);
        System.out.printf("Leader constraints: %d, Follower constraints: %d\n", nConsL, nConsF);

        System.out.println("\n=== LEADER PROBLEM ===");
        System.out.println("Objective: min " + linearExpr(c1leader, "x") + " + " + linearExpr(c2leader, "y"));

        System.out.println("Subject to:");
        for (int i = 0; i < nConsL; i++) {
            String lhs = linearExpr(Aleader[i], "x") + " + " + linearExpr(Bleader[i], "y");
            System.out.println("  " + lhs + " <= " + b[i]);
        }

        System.out.println("\n=== FOLLOWER PROBLEM ===");
        System.out.println("Objective: min " + linearExpr(cfol, "y"));

        System.out.println("Subject to:");
        for (int i = 0; i < nConsF; i++) {
            String lhs = linearExpr(Bfol[i], "x") + " + " + linearExpr(Afol[i], "y");
            System.out.println("  " + lhs + " <= " + bfol[i]);
        }

/*        System.out.println("\n=== VARIABLE TYPES ===");
        for (int i = 0; i < nVarL; i++) {
            System.out.println("x" + i + " (" + leaderVarNames[i] + "): " + typeLeaderStr[i] + " ");
        }
        System.out.println();
        for (int i = 0; i < nVarF; i++) {
            System.out.println("y" + i + " (" + followerVarNames[i] + "): " + typeFollowerStr[i] + " ");
        }
        System.out.println();
*/
        System.out.println("\n=== VARIABLE BOUNDS ===");
        System.out.println("Leader variable bounds:");
        for (int i = 0; i < nVarL; i++) {
            String lbStr = Double.isInfinite(lbLeader[i]) ? (lbLeader[i] < 0 ? "-inf" : "+inf") : String.format("%.2f", lbLeader[i]);
            String ubStr = Double.isInfinite(ubLeader[i]) ? (ubLeader[i] < 0 ? "-inf" : "+inf") : String.format("%.2f", ubLeader[i]);
            System.out.printf("x%d: [%s, %s] (%s)\n", i, lbStr, ubStr, typeLeaderStr[i]);
        }

        System.out.println("Follower variable bounds:");
        for (int i = 0; i < nVarF; i++) {
            String lbStr = Double.isInfinite(lbFollower[i]) ? (lbFollower[i] < 0 ? "-inf" : "+inf") : String.format("%.2f", lbFollower[i]);
            String ubStr = Double.isInfinite(ubFollower[i]) ? (ubFollower[i] < 0 ? "-inf" : "+inf") : String.format("%.2f", ubFollower[i]);
            System.out.printf("y%d: [%s, %s] (%s)\n", i, lbStr, ubStr, typeFollowerStr[i]);
        }
/*
        System.out.println("\n=== VARIABLE NAMES ===");
        System.out.println("Leader variable names:");
        for (int i = 0; i < nVarL; i++) {
            System.out.println("x" + i + ": " + leaderVarNames[i]);
        }

        System.out.println("Follower variable names:");
        for (int i = 0; i < nVarF; i++) {
            System.out.println("y" + i + ": " + followerVarNames[i]);
        }
*/        
        System.out.println("========================");
    }

    // Helper method to create a linear expression string
    private String linearExpr(double[] coeffs, String varPrefix) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;

        for (int i = 0; i < coeffs.length; i++) {
            double coeff = coeffs[i];
            if (coeff == 0) continue;

            if (!first) {
                sb.append(coeff > 0 ? " + " : " - ");
            } else {
                if (coeff < 0) sb.append("-");
            }

            double abs = Math.abs(coeff);
            if (abs != 1) sb.append(abs);
            sb.append(varPrefix).append(i);
            first = false;
        }

        if (first) return "0";
        return sb.toString();
    }

    // Helper method to print a vector
    private void printVector(double[] vec) {
        for (double val : vec) System.out.print(val + " ");
        System.out.println();
    }

    
    // Helper method to print a matrix
    private void printMatrix(double[][] mat) {
        for (double[] row : mat) {
            for (double val : row) System.out.print(val + " ");
            System.out.println();
        }
    }
    
    // Helper method to read a string vector
    private static String[] readStringVector(BufferedReader reader, int len) throws IOException {
        String line;
        do {
            line = reader.readLine();
            if (line == null)
                throw new IOException("Unexpected end of file while reading string vector.");
            line = line.trim();
        } while (line.isEmpty());

        String[] tokens = line.split("\\s+");
        if (tokens.length != len)
            throw new IOException("Expected " + len + " type entries, got " + tokens.length);

        return tokens;
    }

    // Helper method to parse bounds from a string
    private static double parseBound(String token) {
        if (token.equalsIgnoreCase("inf")) return Double.POSITIVE_INFINITY;
        if (token.equalsIgnoreCase("-inf")) return Double.NEGATIVE_INFINITY;
        return Double.parseDouble(token);
    }

    // Helper method to check if a double value is infinite
    public static boolean isInfinity(double val) {
        return Double.isInfinite(val);
    }
}
