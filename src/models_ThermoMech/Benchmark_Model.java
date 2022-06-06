package models_ThermoMech;
/* This Benchmark model was just used to check whether the Matrices and 
	results are correct */
import fem_ThermoMech.Constraint;
import fem_ThermoMech.Force;
import fem_ThermoMech.Node;
import fem_ThermoMech.Structure;

public class Benchmark_Model {
	
	public static Structure createStructure () {
		Structure twoDTruss = new Structure("ThermoMech");
		Node n1 = twoDTruss.addNode(0, 0, 0);
		Node n2 = twoDTruss.addNode(2, 0, 0);
		Node n3 = twoDTruss.addNode(1, 1, 0);
		Constraint c1 = new Constraint(false, false, false);
		c1.setNodalTemp(40);
		n1.setConstraint(c1);
		
		Constraint c2 = new Constraint(true, false, false);
		c2.setNodalTemp(80);
		/*Flux f1 = new Flux(10);
		n2.setFlux(f1);*/
		n2.setConstraint(c2);
		Force f1 = new Force(2.5, 2.5, 0);
		n1.setForce(f1);
		Force f2 = new Force(0, 5, 0);
		n3.setForce(f2);
		Force f3 = new Force(-2.5, 2.5, 0);
		n2.setForce(f3);
		Constraint c3 = new Constraint(true, true, false);
		n3.setConstraint(c3);
		double e1 = 1e3, a = 0.1, e2 = 1414.21, e3 = 2121.32;
		double l1 = 4, l2 = 1.41, l3 = 0.71;
		double tetaRef = 20; double alpha = 5e-4;
		twoDTruss.setTetaRef(tetaRef);
		twoDTruss.addElement(e1, a, l1, alpha, 0, 1);
		twoDTruss.addElement(e2, a, l2, alpha, 0, 2);
		twoDTruss.addElement(e3, a, l3, alpha, 1, 2);
//		twoDTruss.addElement(e1, a, l1, alpha, 0, 1);
//		twoDTruss.addElement(e2, a, l2, alpha, 0, 2);
		// return the new structure
		return twoDTruss ;
	}
	
	public static void main(String[] args) {
		Structure struct = createStructure ();
		// TODO Auto-generated method stub
		struct.printStructure();
		struct.solve();
		struct.printResults();
	
	}

}
