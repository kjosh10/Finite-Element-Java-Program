package models_ThermoMech;

import fem_ThermoMech.Constraint;
import fem_ThermoMech.Force;
import fem_ThermoMech.Node;
import fem_ThermoMech.Structure;
import fem_ThermoMech.Visualizer;
import inf.v3d.view.Viewer;

public class SmallTetraeder_Original_Test_Program {
	public static Structure createStructure () {
		Structure struct = new Structure ("Mech");
		double lb = 15.0;
		double r = 457.2 / 2000;
		double t = 10.0 / 1000;
		double a = Math . PI * ( Math . pow (r , 2) - Math . pow ( r - t , 2));
		double e = 2.1e11 ;
		Constraint c1 = new Constraint ( false , false , false );
		Constraint c2 = new Constraint ( true , true , false );
		Force f = new Force (0 , -20e3 , -100e3 );
		
		
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
		struct . addElement (e , a , 0, 1);
		struct . addElement (e , a , 0, 2);
		struct . addElement (e , a , 0, 3);
		struct . addElement (e , a , 1, 2);
		struct . addElement (e , a , 2, 3);
		struct . addElement (e , a , 3, 1);
		// return the new structure
		return struct ;
	}
	public static void main ( String [] args ) {
	Structure struct = createStructure ();
	// solve
	struct.solve();
	struct.printStructure();
	struct.printResults();
	Viewer viewer = new Viewer ();
	Visualizer viz = new Visualizer ( struct , viewer, "Mech" );
	viz . setConstraintSymbolScale (1);
	viz . setDisplacementScale (1e4);
	viz . setForceSymbolScale (4e-5);
	viz . setForceSymbolRadius (0.075);
	viz . setElementNormalForcesScale (0.00003);
	viz . drawElements ();
	viz . drawConstraints ();
	viz . drawElementForces ();
	viz . drawElementInternalForceColored ();
	viz . drawDisplacements ();
	viz . drawLegendElemForcesColored();
	viewer . setVisible ( true );
	}
}
