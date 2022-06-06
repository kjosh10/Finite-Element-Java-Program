package fem_ThermoMech;

public class Constraint {
	private boolean[] free = new boolean[3];
	private String[] free_string = new String[3]; 
	// To save constraints as either fixed or free.
	private double[] nonHomDBC; // Non Homogenous Dirichlet Boundary condition
	private double prescNodalTemp = -274; 
	/* This temperature is to differentiate while this temperature 
	 * is not possible between prescribed and unprescribed nodes */
	// Constructor with homogeneous DBC
	public Constraint(boolean u1, boolean u2, boolean u3) {
		/* Three constraints in 3 directions where True is free 
		 * and False is fixed */
		// TODO Auto-generated constructor stub
		free[0] = u1;
		if (u1) {
			free_string[0] = "free";
		} else {
			free_string[0] = "fixed";
		}
		free[1] = u2;
		if (u2) {
			free_string[1] = "free";
		} else {
			free_string[1] = "fixed";
		}
		free[2] = u3;
		if (u3) {
			free_string[2] = "free";
		} else {
			free_string[2] = "fixed";
		}
		nonHomDBC = null;
	}
	// Constructor with inHom DBC
	public Constraint(boolean u1, boolean u2, boolean u3, double[] nonHomDBC) {
		/* Here array is of prescribed displacements and 0 where
		   their is no prescribed displacement or fixed */
		// TODO Auto-generated constructor stub
		this.nonHomDBC = nonHomDBC;
		free[0] = u1;
		if (u1 == false) { // False == Fixed
			free_string[0] = "fixed";
		} else if (u1 == true && nonHomDBC[0] != 0) {
			// If displacement is prescribed and the direction is free 
			free_string[0] = "presc";
		} else {
			free_string[0] = "free";
		}

		free[1] = u2;
		if (u2 == false) {
			free_string[1] = "fixed";
		} else if (u2 == true && nonHomDBC[1] != 0) {
			free_string[1] = "presc";
		} else {
			free_string[1] = "free";
		}

		free[2] = u3;
		if (u3 == false) {
			free_string[2] = "fixed";
		} else if (u3 == true && nonHomDBC[2] != 0) {
			free_string[2] = "presc";
		} else {
			free_string[2] = "free";
		}
	}
	// Check whether is free mechanically in direction X_(c+1).
	public boolean isFree(int c) {
		return this.free[c];
	}
	
	// To set the nodal temperature
	public void setNodalTemp(double temp) {
		this.prescNodalTemp = temp;
	}
	
	// To get the nodal temperature
	public double getNodalTemp() {
		return this.prescNodalTemp;
	}
	
	// To get the string of presc, fixed and free.
	public String[] getStringNonHomDBC() {
		return this.free_string;
	}
	
	// To get the array of non Homogeneous DBC
	public double[] getNonHomDBC() {
		return this.nonHomDBC;
	}
	
	// To print accordingly
	public void print() {
		if (this.prescNodalTemp == -274) {
			System.out.printf("%15s%15s%15s%15s\n", free_string[0],
					free_string[1], free_string[2], "unprescTemp");
		} else {
			System.out.printf("%15s%15s%15s%15s\n", free_string[0],
					free_string[1], free_string[2], "  prescTemp");
		}
		
	}
}
