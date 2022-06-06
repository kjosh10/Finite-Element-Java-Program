package fem_ThermoMech;


import iceb.jnumerics.MatrixFormat;
import iceb.jnumerics.Vector3D;

public class Node {
	private int[] dofNumbers = new int[3]; // Array of enumerated DoFs
	private double[] reactionForces; // Array of the reaction forces
	private Vector3D position, displacement; 
	// Vector of position and displacement of the node
	private Constraint constraint = new Constraint(true, true, true); 
	/* Constraint on the node and default is free node if no constrained is defined 
	 on the node */
	private double nodalTemp; // Nodal Temperature
	private Force force; // Force on the node
	private Flux nodalFlux; // Flux on the node
	private double[] nonHomDBC; 
	// non Homogenous Dirichlet Boundary condition (DBC) on the node
	
	public Node(double x1, double x2, double x3) {
		// TODO Auto-generated constructor stub
		position = new Vector3D(x1, x2, x3); // set the position
	}
	// To get position of the node
	public Vector3D getPosition() {
		return position;
	}
	// To set the constraint
	public void setConstraint(Constraint c) {
		this . constraint = c;
		this . nonHomDBC = c.getNonHomDBC();
		this . nodalTemp = c.getNodalTemp();
	}
	
	// To get non homogenous DBC
	public double[] getNonHomDBC() {
		return this.nonHomDBC;
	}
	
	// To get Constraint
	public Constraint getConstraint() {
		return this.constraint;
	}
	
	// To set Flux
	public void setFlux(Flux f) {
		this.nodalFlux = f;
	}
	
	// To get force
	public Flux getFlux() {
		return this.nodalFlux;
	}
	
	// To set Force
	public void setForce(Force f) {
		this.force = f;
	}
	
	// To get force
	public Force getForce() {
		return this.force;
	}
	// To set Reaction Forces
	public void setReactionForces(double[] reactionForces) {
		this.reactionForces = new double[] {reactionForces[0],
				reactionForces[1], reactionForces[2]};
	}
	// To get Reaction Forces
	public double[] getReactionForces() {
		return this.reactionForces;
	}
	
	// To enumerate DOFs
	public int enumerateDOFs(int start) {
		int end = start;
		if (this.constraint == null) {
			dofNumbers[0] = start; 
			dofNumbers[1] = start + 1;
			dofNumbers[2] = start + 2;
			end = start + 3;
		} else {
			if (this.constraint.isFree(0) == true) {
				dofNumbers[0] = end;++end;
				} else {dofNumbers[0] = -1;}
			if (this.constraint.isFree(1) == true) {
				dofNumbers[1] = end;++end;
				} else {dofNumbers[1] = -1;}
			if (this.constraint.isFree(2) == true) {
				dofNumbers[2] = end;++end;
				} else {dofNumbers[2] = -1;}
		}
		return end;
	}
	// To get enumerated DoF numbers
	public int[] getDofNumbers() {
		return dofNumbers;
	}
	// To set Displacement
	public void setDisplacement(double[] u) {
		this.displacement = new Vector3D(u[0], u[1], u[2]);
	}
	// To get Displacement
	public Vector3D getDisplacement() {
		return this.displacement;
	}
	// To get Nodal Temperature
	public double getNodalTemp() {
		return this.nodalTemp;
	}
	// To print
	public void print() {
		System.out.println ( MatrixFormat.format(this.position) );	
	}
}
