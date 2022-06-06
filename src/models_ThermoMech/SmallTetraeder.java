package models_ThermoMech;


import fem_ThermoMech.Constraint;
import fem_ThermoMech.Force;
import fem_ThermoMech.Node;
import fem_ThermoMech.Structure;
import fem_ThermoMech.Visualizer;
import inf.v3d.view.Viewer;

public class SmallTetraeder {
	public static Structure createStructure () {
		Structure struct = new Structure ("ThermoMech");
		double lb = 15.0;
		double r = 457.2 / 2000;
		double t = 10.0 / 1000;
		double a = Math . PI * ( Math . pow (r , 2) - Math . pow ( r - t , 2));
		double e = 2.1e11 ;
		
		Constraint c1 = new Constraint ( false , false , false );
		c1.setNodalTemp(40);
		double[] nonHOMBDC1 = new double[] {-1e-4, 1e-4, 0};
		Constraint c2 = new Constraint ( true , true , false, nonHOMBDC1);
		Force f = new Force (0 , -20e3 , -100e3 );
		c2.setNodalTemp(80);
		// create nodes
		Node n1 = struct . addNode (0.0 , 0.0 , lb * Math . sqrt (2.0 / 3.0));
		Node n2 = struct . addNode (0.0 , lb / Math . sqrt (3) , 0);
		Node n3 = struct . addNode (- lb / 2, -lb / Math . sqrt (12.0) , 0);
		Node n4 = struct . addNode ( lb / 2, - lb / Math . sqrt (12.0) , 0);
		// apply BCs
		n1 . setForce ( f );
		n2 . setConstraint ( c1 );
		n3 . setConstraint ( c1 );
		n4 . setConstraint ( c2 );
		// create elements
		struct . addElement (e , a , 4 , 5e-4 , 0 , 1);
		struct . addElement (e , a , 4 , 5e-4 , 0 , 2);
		struct . addElement (e , a , 4 , 5e-4 , 0 , 3);
		struct . addElement (e , a , 4 , 5e-4 , 1 , 2);
		struct . addElement (e , a , 4 , 5e-4 , 2 , 3);
		struct . addElement (e , a , 4 , 5e-4 , 3 , 1);
		
		
		// return the new structure
		return struct ;
	}
	public static void main ( String [] args ) {
	Structure struct = createStructure ();
	// solve
	struct.printStructure();
	struct.solve();
	struct.printResults();
	Viewer viewer = new Viewer ();
	Visualizer viz = new Visualizer ( struct , viewer , "ThermoMech");
	viz . setConstraintSymbolScale (1);
	viz . setDisplacementScale (9); 	
	viz . setForceSymbolScale (9e-5);	
	viz . setForceSymbolRadius (0.09); 
	viz . setElementNormalForcesScale (1e-7); 
	viz . setNonHomDBCScale(2e4);
	viz . drawElements ();
	viz . drawConstraints ();
	viz . drawElementForces ();
	viz . drawNonHomDBC();
	viz . drawElementInternalForceColored ();
	viz . drawDisplacements ();
	viz . drawLegendElemForcesColored ();
	viz . drawLegendNodalTemp ();
	viz . drawTemperatures();
	viewer . setVisible ( true );
	}
}
