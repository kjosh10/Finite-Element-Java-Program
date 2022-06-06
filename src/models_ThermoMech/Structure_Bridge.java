package models_ThermoMech;

import fem_ThermoMech.Constraint;
import fem_ThermoMech.Force;
import fem_ThermoMech.Node;
import fem_ThermoMech.Structure;
import fem_ThermoMech.Visualizer;
import inf.v3d.view.Viewer;

public class Structure_Bridge {

	public static void main(String[] args)
	{
		String state = "ThermoMech";
		// State can be Therm, Mech, ThermoMech
		Structure S = new Structure(state);
//		double a = Math.PI*(Math.pow(r,2)-Math.pow(r-t, 2));
		double a = 0.05;
		double e = 2.1e11;
		double lambda = 60; //4
		double tetaRef = 10; double alpha = 1.1e-6; // 5e-4
		
		Constraint c1 = new Constraint(false, false, false);
		Constraint c2 = new Constraint(true, false, false);
		Constraint c3 = new Constraint(true, true, false);
//		double[] nonHOMBDC1 = new double[] {1e-7, 0, 0};
//		Constraint c4 = new Constraint ( true , false , false, nonHOMBDC1);
		
		Force f1 = new Force(0, -2000, 0);
		
		
		// create nodes
		
		Node n1 = S.addNode(0, 0, 0);
		Node n2 = S.addNode(20, 30, 0);
		Node n3 = S.addNode(20, 0, 0);
		Node n4 = S.addNode(40, 30, 0);
		Node n5 = S.addNode(40, 0, 0);
		Node n6 = S.addNode(60, 30, 0);
		Node n7 = S.addNode(60, 0, 0);
		Node n8 = S.addNode(80, 30, 0);
		Node n9 = S.addNode(80, 0, 0);
		Node n10 = S.addNode(100, 30, 0);
		Node n11 = S.addNode(100, 0, 0);
		Node n12 = S.addNode(120, 0, 0);
		
		Node n13 = S.addNode(0, 0, -30);
		Node n14 = S.addNode(20, 30, -30);
		Node n15 = S.addNode(20, 0, -30);
		Node n16 = S.addNode(40, 30, -30);
		Node n17 = S.addNode(40, 0, -30);
		Node n18 = S.addNode(60, 30, -30);
		Node n19 = S.addNode(60, 0, -30);
		Node n20 = S.addNode(80, 30, -30);
		Node n21 = S.addNode(80, 0, -30);
		Node n22 = S.addNode(100, 30, -30);
		Node n23 = S.addNode(100, 0, -30);
		Node n24 = S.addNode(120, 00, -30);
		
		n3.setForce(f1);
		n15.setForce(f1);
		n5.setForce(f1);
		n17.setForce(f1);
		n7.setForce(f1);
		n19.setForce(f1);
		n9.setForce(f1);
		n21.setForce(f1);
		n11.setForce(f1);
		n23.setForce(f1);
		
		n1.setConstraint(c1);
		n3.setConstraint(c3);
		n5.setConstraint(c3);
		n7.setConstraint(c3);
		n9.setConstraint(c3);
		n11.setConstraint(c3);
		n12.setConstraint(c2); 
		n2.setConstraint(c3);
		n4.setConstraint(c3);
		n6.setConstraint(c3);
		n8.setConstraint(c3);
		n10.setConstraint(c3);
		
		n13.setConstraint(c1);
		n14.setConstraint(c3);
		n15.setConstraint(c3);
		n17.setConstraint(c3);
		n16.setConstraint(c3);
		n18.setConstraint(c3);
		n19.setConstraint(c3);
		n20.setConstraint(c3);
		n21.setConstraint(c3);
		n22.setConstraint(c3);
		n23.setConstraint(c3);
		n24.setConstraint(c2); 
		
		if ( state == "Mech" ) {
			S.addElement(e, a, 0, 1);
			S.addElement(e, a, 0, 2);
			S.addElement(e, a, 1, 2);
			S.addElement(e, a, 1, 3);
			S.addElement(e, a, 2, 4);
			S.addElement(e, a, 1, 4);
			S.addElement(e, a, 3, 4);
			S.addElement(e, a, 4, 5);
			S.addElement(e, a, 3, 5);
			S.addElement(e, a, 4, 6);
			S.addElement(e, a, 5, 6);
			S.addElement(e, a, 5, 8);
			S.addElement(e, a, 5, 7);
			S.addElement(e, a, 6, 8);
			S.addElement(e, a, 7, 8);
			S.addElement(e, a, 7, 9);
			S.addElement(e, a, 8, 9);
			S.addElement(e, a, 8, 10);
			S.addElement(e, a, 9, 10);
			S.addElement(e, a, 9, 11);
			S.addElement(e, a, 10, 11);
			S.addElement(e, a, 2, 3);
			S.addElement(e, a, 3, 6);
			S.addElement(e, a, 6, 7);
			S.addElement(e, a, 7, 10);
			
			S.addElement(e, a, 12, 13);
			S.addElement(e, a, 12, 14);
			S.addElement(e, a, 13, 14);
			S.addElement(e, a, 13, 15);
			S.addElement(e, a, 13, 16);
			S.addElement(e, a, 14, 16);
			S.addElement(e, a, 15, 16);
			S.addElement(e, a, 15, 17);
			S.addElement(e, a, 16, 17);
			S.addElement(e, a, 16, 18);
			S.addElement(e, a, 17, 18);
			S.addElement(e, a, 17, 19);
			S.addElement(e, a, 17, 20);
			S.addElement(e, a, 18, 20);
			S.addElement(e, a, 19, 20);
			S.addElement(e, a, 19, 21);
			S.addElement(e, a, 20, 21);
			S.addElement(e, a, 20, 22);
			S.addElement(e, a, 21, 22);
			S.addElement(e, a, 21, 23);
			S.addElement(e, a, 22, 23);
			S.addElement(e, a, 14, 15);
			S.addElement(e, a, 15, 18);
			S.addElement(e, a, 18, 19);
			S.addElement(e, a, 19, 22);
			
			S.addElement(e, a, 0, 12);
			S.addElement(e, a, 1, 13);
			S.addElement(e, a, 2, 14);
			S.addElement(e, a, 3, 15);
			S.addElement(e, a, 4, 16);
			S.addElement(e, a, 5, 17);
			S.addElement(e, a, 6, 18);
			S.addElement(e, a, 7, 19);
			S.addElement(e, a, 8, 20);
			S.addElement(e, a, 9, 21);
			S.addElement(e, a, 10, 22);
			S.addElement(e, a, 11, 23);
			
			S.addElement(e, a, 1, 15);
			S.addElement(e, a, 3, 17);
			S.addElement(e, a, 5, 19);
			S.addElement(e, a, 7, 21);
		
			S.printStructure();
			S.solve();
			System.out.print("\n");
			S.printResults();
			Viewer v = new Viewer();
			//System.out.print(force);
			Visualizer p = new Visualizer(S, v, state);
			
			p . drawElements();
			p . setConstraintSymbolScale (3);
			p . setForceSymbolScale  (6e-3);
			p . setForceSymbolRadius (0.5);
			p . setElementNormalForcesScale (0.0003);
			p . setNonHomDBCScale(1e8);
			p . drawElements();
			p . drawConstraints();
			p . drawElementForces();
			p . drawElementNormalForces();
			p . drawNonHomDBC();
			p . setDisplacementScale(3e4);
			p . drawDisplacements();
			p . drawElementInternalForceColored();
			p . drawLegendElemForcesColored();
			v . setVisible(true);
			
		} else {  

			S . setTetaRef(tetaRef);
			Constraint c5 = new Constraint(false, false, false);
			// Prescribing Temperature
			c5.setNodalTemp(20);
			double[] nonHOMBDC2 = new double[] {1e-7, 0, 0};
			Constraint c6 = new Constraint(true, false, false, nonHOMBDC2);
			
			c6.setNodalTemp(30);
			n1.setConstraint(c5);
			n13.setConstraint(c5);
			n12.setConstraint(c6);
			n24.setConstraint(c6);
			
			S.addElement(e, a, lambda, alpha, 0, 1);
			S.addElement(e, a, lambda, alpha, 0, 2);
			S.addElement(e, a, lambda, alpha, 1, 2);
			S.addElement(e, a, lambda, alpha, 1, 3);
			S.addElement(e, a, lambda, alpha, 2, 4);
			S.addElement(e, a, lambda, alpha, 1, 4);
			S.addElement(e, a, lambda, alpha, 3, 4);
			S.addElement(e, a, lambda, alpha, 4, 5);
			S.addElement(e, a, lambda, alpha, 3, 5);
			S.addElement(e, a, lambda, alpha, 4, 6);
			S.addElement(e, a, lambda, alpha, 5, 6);
			S.addElement(e, a, lambda, alpha, 5, 8);
			S.addElement(e, a, lambda, alpha, 5, 7);
			S.addElement(e, a, lambda, alpha, 6, 8);
			S.addElement(e, a, lambda, alpha, 7, 8);
			S.addElement(e, a, lambda, alpha, 7, 9);
			S.addElement(e, a, lambda, alpha, 8, 9);
			S.addElement(e, a, lambda, alpha, 8, 10);
			S.addElement(e, a, lambda, alpha, 9, 10);
			S.addElement(e, a, lambda, alpha, 9, 11);
			S.addElement(e, a, lambda, alpha, 10, 11);
			S.addElement(e, a, lambda, alpha, 2, 3);
			S.addElement(e, a, lambda, alpha, 3, 6);
			S.addElement(e, a, lambda, alpha, 6, 7);
			S.addElement(e, a, lambda, alpha, 7, 10);
			
			S.addElement(e, a, lambda, alpha, 12, 13);
			S.addElement(e, a, lambda, alpha, 12, 14);
			S.addElement(e, a, lambda, alpha, 13, 14);
			S.addElement(e, a, lambda, alpha, 13, 15);
			S.addElement(e, a, lambda, alpha, 13, 16);
			S.addElement(e, a, lambda, alpha, 14, 16);
			S.addElement(e, a, lambda, alpha, 15, 16);
			S.addElement(e, a, lambda, alpha, 15, 17);
			S.addElement(e, a, lambda, alpha, 16, 17);
			S.addElement(e, a, lambda, alpha, 16, 18);
			S.addElement(e, a, lambda, alpha, 17, 18);
			S.addElement(e, a, lambda, alpha, 17, 19);
			S.addElement(e, a, lambda, alpha, 17, 20);
			S.addElement(e, a, lambda, alpha, 18, 20);
			S.addElement(e, a, lambda, alpha, 19, 20);
			S.addElement(e, a, lambda, alpha, 19, 21);
			S.addElement(e, a, lambda, alpha, 20, 21);
			S.addElement(e, a, lambda, alpha, 20, 22);
			S.addElement(e, a, lambda, alpha, 21, 22);
			S.addElement(e, a, lambda, alpha, 21, 23);
			S.addElement(e, a, lambda, alpha, 22, 23);
			S.addElement(e, a, lambda, alpha, 14, 15);
			S.addElement(e, a, lambda, alpha, 15, 18);
			S.addElement(e, a, lambda, alpha, 18, 19);
			S.addElement(e, a, lambda, alpha, 19, 22);
			
			S.addElement(e, a, lambda, alpha, 0, 12);
			S.addElement(e, a, lambda, alpha, 1, 13);
			S.addElement(e, a, lambda, alpha, 2, 14);
			S.addElement(e, a, lambda, alpha, 3, 15);
			S.addElement(e, a, lambda, alpha, 4, 16);
			S.addElement(e, a, lambda, alpha, 5, 17);
			S.addElement(e, a, lambda, alpha, 6, 18);
			S.addElement(e, a, lambda, alpha, 7, 19);
			S.addElement(e, a, lambda, alpha, 8, 20);
			S.addElement(e, a, lambda, alpha, 9, 21);
			S.addElement(e, a, lambda, alpha, 10, 22);
			S.addElement(e, a, lambda, alpha, 11, 23);
			
			S.addElement(e, a, lambda, alpha, 1, 15);
			S.addElement(e, a, lambda, alpha, 3, 17);
			S.addElement(e, a, lambda, alpha, 5, 19);
			S.addElement(e, a, lambda, alpha, 7, 21);
			
			S.printStructure();
			S.solve();
			System.out.print("\n");
			S.printResults();
			Viewer v = new Viewer();
			Visualizer p = new Visualizer(S, v, state);
			
			p . setConstraintSymbolScale (3);
			p . setDisplacementScale (3e3); // 6
			p . setForceSymbolScale  (6e-3);
			p . setForceSymbolRadius (0.5);
			p . setElementNormalForcesScale (0.0003);
			p . setNonHomDBCScale(1e8);
			p . drawElements();
			p . drawConstraints();
			if (state != "Therm") {
				p . drawElementInternalForceColored();
				p . drawLegendElemForcesColored();
				p . drawNonHomDBC();
				p . drawDisplacements();
			}
			
			p . drawElementForces();
			p . drawTemperatures();
			p . drawLegendNodalTemp();
			v . setVisible(true);
		}
		
	}
}

