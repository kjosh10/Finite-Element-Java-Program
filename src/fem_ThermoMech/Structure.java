package fem_ThermoMech;

import java.util.ArrayList;

import iceb.jnumerics.Array2DMatrix;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.MatrixFormat;
import iceb.jnumerics.QuadraticMatrixInfo;
import iceb.jnumerics.SolveFailedException;
import iceb.jnumerics.lse.GeneralMatrixLSESolver;
import iceb.jnumerics.lse.ILSESolver;
import inf.text.ArrayFormat;

public class Structure {
	private ArrayList<Node> nodes = new ArrayList<Node>();
	private ArrayList<Node> nodes_Ordered = new ArrayList<Node>();
	private ArrayList<Element> elements = new ArrayList<Element>();
	private double tetaRef = 0; // Default value is 0
	private String state;

	// To set state of the model which can be Mech, Therm or ThermoMech
	public Structure(String stateModel) {
		this.state = stateModel;
	}

	// Set Teta Reference
	public void setTetaRef(double tetaReference) {
		this.tetaRef = tetaReference;
	}

	// Get Teta Reference
	public double getTetaRef() {
		return this.tetaRef;
	}

	// To add the node in the structure
	public Node addNode(double x1, double x2, double x3) {
		Node n = new Node(x1, x2, x3);
		nodes.add(n);
		return n;
	}

	// To add the element for Thermal or Thermo-Mechanical Problem
	public Element addElement(double e, double a, double lambda, double alpha, int n1, int n2) {
		Element elem = new Element(e, a, lambda, alpha, nodes.get(n1), nodes.get(n2));
		elements.add(elem);
		return elem;
	}

	// To add the element for Mechanical problem
		public Element addElement(double e, double a, int n1, int n2) {
			Element elem = new Element(e, a, nodes.get(n1), nodes.get(n2));
			elements.add(elem);
			return elem;
		}
	// To get the total number of Nodes
	public int getNumberOfNodes() {
		return this.nodes.size();
	}

	// To get the ith Node
	public Node getNode(int id) {
		return this.nodes.get(id);
	}

	// To get the total number of Elements
	public int getNumberOfElements() {
		return elements.size();
	}

	// To get the ith Element
	public Element getElement(int id) {
		return elements.get(id);
	}

	// To print the Structure
	public void printStructure() {
		System.out.println("------------------------------------------------------------"
				+ "------------------------------------------------");
		System.out.println(" Listing Structure ");
		System.out.println("------------------------------------------------------------"
				+ "------------------------------------------------");
		System.out.print("\n");
		
		// For printing Nodes
		System.out.println(" Nodes ");
		System.out.printf(" idx%13s%15s%15s", "x1", "x2", "x3");
		System.out.print("\n");
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			System.out.printf("%5d", i);
			this.getNode(i).print();
		}
		
		// For printing constraints
		System.out.print("\n");
		System.out.println(" Constraints ");
		System.out.printf(" node%12s%15s%15s", "u1", "u2", "u3");
		System.out.print("\n");
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getConstraint() != null) {
				System.out.printf("%5d", i);
				this.getNode(i).getConstraint().print();
			}

		}
		
		// For printing Forces
		if (this.state == "Mech" || this.state == "ThermoMech") {
			System.out.print("\n");
			System.out.println(" Forces ");
			System.out.printf(" node%12s%15s%15s", "r1", "r2", "r3");
			System.out.print("\n");
			for (int i = 0; i < this.getNumberOfNodes(); i++) {
				if (this.getNode(i).getForce() != null) {
					System.out.printf("%5d", i);
					this.getNode(i).getForce().print();
				}

			}

		}
		
		// For printing Elements
		System.out.print("\n");
		System.out.println(" Elements ");
		System.out.printf(" idx%13s%15s%15s%18s%15s", "E", "A", "Lambda", 
				"Conductivity", "length");
		System.out.print("\n");
		for (int i = 0; i < this.getNumberOfElements(); i++) {

			System.out.printf("%5d", i);
			this.getElement(i).print();
		}
		System.out.print("\n");
	}

	// To solve the system
	public void solve() {
		// To order nodes according to the elements defined
		this.orderingNodes(); 
		if (this.state == "Mech") { // Use Mechanical solver
			this.mechSolver();
		} else if (this.state == "Therm") { // Use Thermal solver
			this.thermalSolver();
		} else if (this.state == "ThermoMech") { // Use ThermoMech solver
			this.thermoMechSolver();
		}
	}

	// Mechanical solver
	public void mechSolver() {
		// size of our matrix
		int neq = this.enumerateDOFs();
		boolean nonHomFlag = false; 
		// Boolean which is true if their is atleast one Non Homogeneous DBC
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getNonHomDBC() != null) {
				nonHomFlag = true; // Their is atleast one NonHomDBC
				break;
			}
		}

		if (nonHomFlag == false) { // Solver for Homogeneous Model
			// create the solver object
			ILSESolver solver = new GeneralMatrixLSESolver();
			// info object for coefficient matrix
			QuadraticMatrixInfo aInfo = solver.getAInfo();
			// get coefficient matrix
			IMatrix kGlobal = solver.getA();
			// right hand side
			double[] rGlobal = new double[neq];
			double[] uGlobal = rGlobal; // As solution will replace right hand side.
			// initialize solver
			aInfo.setSize(neq);
			solver.initialize();
			double[] dispPresc = new double[neq];
			assembleStiffnessMatrix(kGlobal, dispPresc);

			assembleLoadVector(uGlobal);

			try { // To catch solver failed error
				solver.solve(uGlobal);
				// System.out.println(ArrayFormat.format(solution));
			} catch (SolveFailedException e) {
				System.out.println(" Solve failed : " + e.getMessage());
			}
			this.selectDisplacements(uGlobal);
			this.selectReactionForces(this.calcReactions(kGlobal, dispPresc, uGlobal));
		
		} else { // If their is a NonHom DBC, then do static condensation

			IMatrix kGlobal = new Array2DMatrix(neq, neq);
			double[] rGlobal = new double[neq];
			double[] uGlobal = rGlobal; // As solution will replace right hand side.

			double[] dispPresc = new double[neq];
			assembleStiffnessMatrix(kGlobal, dispPresc);
			assembleLoadVector(uGlobal);
			// create the solver object
			ILSESolver solver = new GeneralMatrixLSESolver();
			// info object for coefficient matrix
			QuadraticMatrixInfo aInfo = solver.getAInfo();
			// get coefficient matrix
			IMatrix kGlobalReduced = solver.getA();
			// right hand side
			double[] rGlobalReduced = new double[neq - this.getTotalInHomDBC(this.state)];

			// initialize solver
			aInfo.setSize(neq - this.getTotalInHomDBC(this.state));
			solver.initialize();
			// Static condensation
			this.staticCondensation(kGlobal, rGlobal, dispPresc, kGlobalReduced, rGlobalReduced);
			/*System.out.print("\n");
			System.out.println(MatrixFormat.format(kGlobal));
			System.out.print("\n");
			System.out.println(MatrixFormat.format(kGlobalReduced));
			System.out.print("\n");
			System.out.println(ArrayFormat.format(rGlobal));
			System.out.print("\n");
			System.out.println(ArrayFormat.format(rGlobalReduced));*/
			double[] uGlobalReduced = rGlobalReduced; // As solution will replace right hand side.
			try {
				solver.solve(uGlobalReduced);
				// System.out.println(ArrayFormat.format(solution));
			} catch (SolveFailedException e) {
				System.out.println(" Solve failed : " + e.getMessage());
			}

			/*System.out.println(ArrayFormat.format(uGlobalReduced));
			System.out.print("\n");
			System.out.println(ArrayFormat.format(this.totalDispWithNHDBC(uGlobalReduced)));
			System.out.print("\n");*/
			this.selectDisplacements(this.totalDispWithNHDBC(uGlobalReduced));
			this.selectReactionForces(this.calcReactions(kGlobal, dispPresc, this.totalDispWithNHDBC(uGlobalReduced)));
		}
	}

	public void thermalSolver() {
		boolean isThermoNonHom = false; // is their any prescribed temperature?
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getConstraint().getNodalTemp() != -274) {
				isThermoNonHom = true; // Their is atleast one Temperature constraint
				break;
			}
		}

		if (isThermoNonHom == false) { // Solver when no temperature is prescribed
			// create the solver object
			ILSESolver solver = new GeneralMatrixLSESolver();
			// info object for coefficient matrix
			QuadraticMatrixInfo aInfo = solver.getAInfo();
			// get coefficient matrix
			IMatrix kGlobal = solver.getA();
			// right hand side
			double[] rGlobal = new double[this.nodes.size()];
			double[] tetaGlobal = rGlobal; // As solution will replace right hand side.
			// initialize solver
			aInfo.setSize(this.nodes.size());
			solver.initialize();
			double[] prescTemp = new double[this.nodes.size()];
			this.assembleTempMatrix(kGlobal, prescTemp);

			this.assembleFluxVector(tetaGlobal);

			try {
				solver.solve(tetaGlobal);
				// System.out.println(ArrayFormat.format(solution));
			} catch (SolveFailedException e) {
				System.out.println(" Solve failed : " + e.getMessage());
			}

			this.selectTemperatures(tetaGlobal);
			this.selectReactionFluxes(this.calcReactions(kGlobal, prescTemp,
					tetaGlobal));
			
		} else { // Else do static condensation

			IMatrix kGlobal = new Array2DMatrix(this.nodes.size(), this.nodes.size());
			double[] rGlobal = new double[this.nodes.size()];
			double[] tetaGlobal = rGlobal; // As solution will replace right hand side.

			double[] tempPresc = new double[this.nodes.size()];
			this.assembleTempMatrix(kGlobal, tempPresc);
			this.assembleFluxVector(tetaGlobal);
			// create the solver object
			ILSESolver solver = new GeneralMatrixLSESolver();
			// info object for coefficient matrix
			QuadraticMatrixInfo aInfo = solver.getAInfo();
			// get coefficient matrix
			IMatrix kGlobalReduced = solver.getA();
			// right hand side
			double[] rGlobalReduced = new double[this.nodes.size() - this.getTotalInHomDBC("Therm")];

			// initialize solver
			aInfo.setSize(this.nodes.size() - this.getTotalInHomDBC("Therm"));
			solver.initialize();
			this.staticCondensation(kGlobal, rGlobal, tempPresc, kGlobalReduced, rGlobalReduced);
			
			double[] tetaGlobalReduced = rGlobalReduced; // As solution will replace right hand side.
			try {
				solver.solve(tetaGlobalReduced); // Was rGlobalReduced
				// System.out.println(ArrayFormat.format(solution));
			} catch (SolveFailedException e) {
				System.out.println(" Solve failed : " + e.getMessage());
			}
			
			this.selectTemperatures(tetaGlobalReduced);
			this.selectReactionFluxes(this.calcReactions(kGlobal, tempPresc,
					this.totalTempWithNHDBC()));
		}
		
	}

	// Thermo-Mechanical Solver
	public void thermoMechSolver() {
		//// ----------Thermal part of the solver---------///////////
		IMatrix kTetaTetaGlobal = new Array2DMatrix(this.nodes.size(), this.nodes.size());
		double[] rTetaGlobal = new double[this.nodes.size()];
		double[] prescTemp = new double[this.nodes.size()];
		this.assembleTempMatrix(kTetaTetaGlobal, prescTemp);
		this.assembleFluxVector(rTetaGlobal);
		//// --------Mechanical part of the solver--------///////////
		// size of our matrix
		int neq = this.enumerateDOFs();
		// Here u represents displacements i.e. Mechanical part
		IMatrix kUUGlobal = new Array2DMatrix(neq, neq);
		double[] rUGlobal = new double[neq];
		double[] prescDisp = new double[neq];
		this.assembleStiffnessMatrix(kUUGlobal, prescDisp);
		this.assembleLoadVector(rUGlobal);
		
		//// -------ThermoMech part of the solver----------/////////
		IMatrix kUTetaGlobal = new Array2DMatrix(neq, this.nodes.size());
		this.assembleDispTempMatrix(kUTetaGlobal);
		
		// Solve thermal part
		this.thermalSolver();
		double[] teta = new double[prescTemp.length];
		double[] tetaRefVector = new double[prescTemp.length];
		// Solution of the thermal part
		for (int i = 0; i < prescTemp.length; i++) {
			teta[i] = this.getNode(i).getConstraint().getNodalTemp();
			// To make a Teta vector and a vector with tetaRef everywhere
			tetaRefVector[i] = this.tetaRef;
		}
		double[] rFinal = new double[rUGlobal.length];
		for (int i = 0; i < rUGlobal.length; i++) {
			rFinal[i] = rUGlobal[i];
			for (int j = 0; j < rTetaGlobal.length; j++) {
				rFinal[i] += kUTetaGlobal.get(i, j) * (teta[j] - tetaRefVector[j]);
			}
		}
		
		// create the solver object
		ILSESolver solver = new GeneralMatrixLSESolver();
		// info object for coefficient matrix
		QuadraticMatrixInfo aInfo = solver.getAInfo();
		// get coefficient matrix
		IMatrix kGlobalReduced = solver.getA();
		// right hand side
		double[] rGlobalReduced = new double[prescDisp.length - this.getTotalInHomDBC("Mech")];
		// While we are now only solving for mechanical DOF's

		// initialize solver
		aInfo.setSize(prescDisp.length - this.getTotalInHomDBC("Mech"));
		solver.initialize();
		
		this.staticCondensation(kUUGlobal, rFinal, prescDisp, kGlobalReduced, rGlobalReduced);
		
		double[] uGlobalReduced = rGlobalReduced; 
		// As solution will replace right hand side.
		try {
			solver.solve(uGlobalReduced);
		} catch (SolveFailedException e) {
			System.out.println(" Solve failed : " + e.getMessage());
		}

		this.selectDisplacements(this.totalDispWithNHDBC(uGlobalReduced));
		this.selectReactionForces(this.calcReactions(kUUGlobal, prescDisp, 
				this.totalDispWithNHDBC(uGlobalReduced)));
		
	}

	// Method to do static condensation
	public void staticCondensation(IMatrix kGlobal, double[] rGlobal, double[] presc, IMatrix kGlobalReduced,
			double[] rGlobalReduced) {
		int reducedCounterColumn = 0, reducedCounterRow = 0;
		double prescValue = 0;

		for (int i = 0; i < kGlobal.getRowCount(); i++) {
			reducedCounterColumn = 0;
			prescValue = 0;
			for (int j = 0; j < kGlobal.getColumnCount(); j++) {
				if (presc[j] != 0) {
					prescValue += presc[j] * kGlobal.get(i, j);
					if (i == j) {
						continue;
					}
				} else {
					if (presc[i] == 0) {
						kGlobalReduced.add(reducedCounterRow, reducedCounterColumn, kGlobal.get(i, j));
						reducedCounterColumn++;
					}

				}
			}
			if (presc[i] == 0) {
				rGlobalReduced[reducedCounterRow] = rGlobal[i] - prescValue;
			}
			reducedCounterRow++;
			if (presc[i] != 0) {
				reducedCounterRow--;
			}
		}
	}

	// Get total number of non homogenous Dirichlet Boundary Condition
	public int getTotalInHomDBC(String stateModel) {
		int countDisp = 0, countTemp = 0;
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			//System.out.println(this.getNode(i).getConstraint().getNodalTemp());
			if (this.getNode(i).getConstraint().getNodalTemp() != -274) {
				countTemp++; // Number of prescribed nodal Temperatures
			}
			for (int j = 0; j < this.getNode(i).getDofNumbers().length; j++) {
				if (this.getNode(i).getNonHomDBC() != null) {
					if (this.getNode(i).getConstraint().getStringNonHomDBC()[j] == "presc") {
						countDisp++; // Number of prescribed displacement
					}
				}
			}
		}

		if (stateModel == "Mech") {
			return countDisp;
		} else if (stateModel == "Therm") {
			return countTemp;
		} else {
			return countDisp + countTemp;
		}

	}

	// To order nodes in the way they are defined while defining elements
	public void orderingNodes() {
		for (int i = 0; i < this.getNumberOfElements(); i++) {
			if (this.nodes_Ordered.contains(this.getElement(i).getNode1()) == false) {
				this.nodes_Ordered.add(this.getElement(i).getNode1());
			}
			if (this.nodes_Ordered.contains(this.getElement(i).getNode2()) == false) {
				this.nodes_Ordered.add(this.getElement(i).getNode2());
			}
		}
		this.nodes = this.nodes_Ordered; // To make nodes ordered as well.
	}

	// To enumerate all the Nodes
	public int enumerateDOFs() {
		int start = 0, neq; // neq -> No. of Equations
		this.orderingNodes();
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			Node n = this.nodes.get(i);
			start = n.enumerateDOFs(start);
		}
		for (int i = 0; i < this.getNumberOfElements(); i++) {
			Element elem = this.getElement(i);
			elem.enumerateDOFs();
		}
		neq = start;
		return neq;
	}

	// To assemble Stiffness Matrix
	public void assembleStiffnessMatrix(IMatrix kGlobal, double[] dispPresc) {
		IMatrix kElement;
		int[] dofNumber;
		for (int i = 0; i < this.getNumberOfElements(); i++) {
			dofNumber = this.getElement(i).getDofNumbers();
			kElement = this.getElement(i).computeStiffnessMatrix();
			for (int m = 0; m < kElement.getRowCount(); m++) {
				for (int n = 0; n < kElement.getColumnCount(); n++) {

					if (dofNumber[m] != -1 && dofNumber[n] != -1) {
						kGlobal.add(dofNumber[m], dofNumber[n], kElement.get(m, n));
					}
				}
			}
		} // End Of Assembling
		/* dispPresc is a matrix with size of all degrees of freedom and
		only prescribed displacements in it for rest DOF displacement is zero.*/
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.nodes_Ordered.get(i).getNonHomDBC() != null) {
				for (int j = 0; j < this.nodes_Ordered.get(i).getNonHomDBC().length; j++) {
					// HEREEEEEEE
					if (this.nodes_Ordered.get(i).getNonHomDBC()[j] != 0) {
						dispPresc[this.nodes_Ordered.get(i).getDofNumbers()[j]] = this.nodes_Ordered.get(i)
								.getNonHomDBC()[j];
					}
				}
			}
		}
	}

	// To assemble Temperature Matrix
	public void assembleTempMatrix(IMatrix tempStiffnessGlobal, double[] prescTemp) {
		IMatrix tempStiffnessElement;
		for (int i = 0; i < this.getNumberOfElements(); i++) {

			tempStiffnessElement = this.getElement(i).computeTempStiffnessMatrix();

			int node1Index = this.nodes.indexOf(this.getElement(i).getNode1());
			int node2Index = this.nodes.indexOf(this.getElement(i).getNode2());
		
			tempStiffnessGlobal.add(node1Index, node1Index, tempStiffnessElement.get(0, 0));
			tempStiffnessGlobal.add(node1Index, node2Index, tempStiffnessElement.get(0, 1));
			tempStiffnessGlobal.add(node2Index, node1Index, tempStiffnessElement.get(1, 0));
			tempStiffnessGlobal.add(node2Index, node2Index, tempStiffnessElement.get(1, 1));
			
		} // End Of Assembling
		// dispTemp is a matrix with size of total Nodes and
		// only prescribed temperatures in it for rest Nodes temperature is zero.
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			// this.orderingNodes();
			if (this.nodes.get(i).getConstraint().getNodalTemp() != -274) {
				prescTemp[i] = this.nodes.get(i).getConstraint().getNodalTemp();
			}
		}
	
	}

	// To assemble combined Matrix
	public void assembleDispTempMatrix(IMatrix tempDispStiffnessGlobal) {
		IMatrix tempDispStiffnessElement;
		int[] dofNumber;
		// System.out.println(MatrixFormat.format(tempDispStiffnessGlobal));
		for (int i = 0; i < this.getNumberOfElements(); i++) {

			tempDispStiffnessElement = this.getElement(i).computeDispTempStiffnessMatrix();

			int node1Index = this.nodes_Ordered.indexOf(this.getElement(i).getNode1());
			int node2Index = this.nodes_Ordered.indexOf(this.getElement(i).getNode2());
			dofNumber = this.getElement(i).getDofNumbers();
			tempDispStiffnessElement = this.getElement(i).computeDispTempStiffnessMatrix();
			
			for (int m = 0; m < tempDispStiffnessElement.getRowCount(); m++) {

				if (dofNumber[m] != -1) {
					tempDispStiffnessGlobal.add(dofNumber[m], node1Index, tempDispStiffnessElement.get(m, 0));
					tempDispStiffnessGlobal.add(dofNumber[m], node2Index, tempDispStiffnessElement.get(m, 1));
				}

			}
		}
	}

	// To assemble Thermal stiffness Matrix
	public void assembleThermoMechVector(IMatrix kFinalGlobal, IMatrix kUUGlobal, IMatrix kUTetaGlobal,
			IMatrix kTetaTetaGlobal) {
		int neq = this.enumerateDOFs();
		for (int i = 0; i < kFinalGlobal.getRowCount(); i++) {
			for (int j = 0; j < kFinalGlobal.getColumnCount(); j++) {
				if (i < neq && j < neq) {
					kFinalGlobal.add(i, j, kUUGlobal.get(i, j));
				} else if (i < neq && j >= neq) {
					kFinalGlobal.add(i, j, kUTetaGlobal.get(i, j - neq) * -1);
				} else if (i >= neq && j >= neq) {
					kFinalGlobal.add(i, j, kTetaTetaGlobal.get(i - neq, j - neq));
				}
			}
		}
	}

	// To assemble Load Vector
	public void assembleLoadVector(double[] rGlobal) {
		// this.orderingNodes();
		int[] dofNumber;
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			dofNumber = this.nodes.get(i).getDofNumbers(); // Was nodes_ordered
			for (int j = 0; j < dofNumber.length; j++) {
				if (dofNumber[j] != -1) {
					if (this.getNode(i).getForce() == null) {
						rGlobal[dofNumber[j]] = 0;
					} else {
						rGlobal[dofNumber[j]] = this.getNode(i).getForce().getComponent(j);
					}
				}
			}
		}
	}

	// To assemble Flux Vector
	public void assembleFluxVector(double[] rGlobal) {
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getFlux() != null) {
				if (this.getNode(i).getFlux().getNodalFlux() != 0) {
					rGlobal[i] = this.getNode(i).getFlux().getNodalFlux();
				}
			}
		}
	}

	// To calculate R_u - KuTeta* (Teta - TetaRef) vector
	public void assembleFinalLoadVector(double[] rFinalGlobal, double[] rUGlobal, double[] rTetaGlobal, IMatrix kUTeta,
			double tetaReference) {
		double[] tetaRefVector = new double[rTetaGlobal.length];
		for (int i = 0; i < tetaRefVector.length; i++) {
			tetaRefVector[i] = tetaRefVector[i] + tetaReference;
		}
		double[] rReduced = new double[rUGlobal.length];
		for (int i = 0; i < rUGlobal.length; i++) {
			rReduced[i] = rUGlobal[i];
			for (int j = 0; j < rTetaGlobal.length; j++) {
				rReduced[i] -= kUTeta.get(i, j) * tetaRefVector[j];
			}
		}
		for (int i = 0; i < rFinalGlobal.length; i++) {
			if (i < rUGlobal.length) {
				rFinalGlobal[i] = rReduced[i];
			} else {
				rFinalGlobal[i] = rTetaGlobal[i - rUGlobal.length];
			}
		}
	}

	// To store all the reaction in their right places
	public void selectReactionForces(double[] reactionForces) {
		double[] reactionForceNodal = new double[3];
		int counter = 0;
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			for (int j = 0; j < this.nodes.get(i).getDofNumbers().length; j++) {
				//System.out.println(counter);
				if (this.getNode(i).getDofNumbers()[j] == -1) { // If constrained
					reactionForceNodal[j] = 0;
				} else {
					reactionForceNodal[j] = reactionForces[counter];
					counter++;
				}
			}
			this.getNode(i).setReactionForces(reactionForceNodal);
		}
	}
	
	// To store all the reaction Fluxes in their right places
	public void selectReactionFluxes(double[] reactionFluxes) {
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getFlux() == null) { // If flux is null
				Flux fNodal = new Flux(reactionFluxes[i]);
				this.getNode(i).setFlux(fNodal);
			}
		}
	}
	
	// To store all the displacements in their right places
	public void selectDisplacements(double[] uGlobal) {
		double[] uNodal = new double[3];
		int counter = 0;
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			for (int j = 0; j < this.nodes_Ordered.get(i).getDofNumbers().length; j++) {
				if (this.getNode(i).getDofNumbers()[j] == -1) { // If constrained
					uNodal[j] = 0;
				} else {
					uNodal[j] = uGlobal[counter];
					counter++;
				}
			}
			this.getNode(i).setDisplacement(uNodal);
		}
	}

	// To store all the Temperatures in their right places
	public void selectTemperatures(double[] tetaGlobal) {
		int counter = 0;
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			if (this.getNode(i).getConstraint().getNodalTemp() == -274) {
				this.getNode(i).getConstraint().setNodalTemp(tetaGlobal[counter]);
				counter++;
			}
		}
	}
	
	// To get total number of non Homogeneous DBC (displacement) in the model
	public double[] totalDispWithNHDBC(double[] uGlobalReduced) {
		int[] dofNumber;
		int counter = 0;
		double[] uGlobal = new double[this.enumerateDOFs()];
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			dofNumber = this.nodes_Ordered.get(i).getDofNumbers();
			for (int j = 0; j < dofNumber.length; j++) {
				if (dofNumber[j] != -1) {
					if (this.getNode(i).getNonHomDBC() == null) {
						uGlobal[dofNumber[j]] = uGlobalReduced[counter];
						counter++;
					} else {
						if (this.getNode(i).getNonHomDBC()[j] != 0) {
							uGlobal[dofNumber[j]] = this.getNode(i).getNonHomDBC()[j];
							
						} else {
							uGlobal[dofNumber[j]] = uGlobalReduced[counter];
							
						}
					}
				}
			}
		}
		return uGlobal;
	}

	// To get total number of non Homogeneous DBC (temperature) in the model
	public double[] totalTempWithNHDBC() {
		double[] tetaGlobal = new double[this.getNumberOfNodes()];
		for (int i = 0; i < this.getNumberOfNodes(); i++) {
			tetaGlobal[i] = this.getNode(i).getConstraint().getNodalTemp();
		}
		return tetaGlobal;
	}

	// To calculate Reaction Forces
	public double[] calcReactionForces(IMatrix kGlobal, double[] dispPresc) {
		double reactionForce[] = new double[dispPresc.length];
		for (int i = 0; i < kGlobal.getRowCount(); i++) {
			int sum = 0;
			for (int j = 0; j < kGlobal.getColumnCount(); j++) {
				if (dispPresc[i] != 0) {
					sum += kGlobal.get(i, j)*dispPresc[j];
				}
			}
			reactionForce[i] = sum;
		}
		return reactionForce;
	}
	
	// To calculate Reaction Forces
		public double[] calcReactions(IMatrix kGlobal, double[] prescDof,
				double[] sol) {
			double reaction[] = new double[prescDof.length];
			for (int i = 0; i < kGlobal.getRowCount(); i++) {
				double sum = 0;
				for (int j = 0; j < kGlobal.getColumnCount(); j++) {
					if (prescDof[i] != 0) {
						sum += kGlobal.get(i, j)*sol[j];
					}
				}
				reaction[i] = sum;
			}
			return reaction;
		}
	
	// To print all the results
	public void printResults() {
		System.out.println("------------------------------------------------------------"
				+ "------------------------------------------------");
		System.out.println(" Results ");
		System.out.println("------------------------------------------------------------"
				+ "------------------------------------------------");
		System.out.print("\n");
		if (this.state == "Mech" || this.state == "ThermoMech") {
			// For printing Displacements
			System.out.println(" Displacements ");
			System.out.printf(" node%14s%15s%15s", "u1", "u2", "u3");
			System.out.print("\n");
			for (int i = 0; i < this.getNumberOfNodes(); i++) {
				System.out.printf("%5d  ", i);
				System.out.println(MatrixFormat.format(this.getNode(i).getDisplacement()));
			}
			// For printing Element Forces
			System.out.print("\n");
			System.out.println(" Element Forces ");
			System.out.printf(" %s%15s", "elem", "force");
			System.out.print("\n");
			for (int i = 0; i < this.getNumberOfElements(); i++) {
				System.out.printf("%5d", i);
				System.out.println(ArrayFormat.format(this.getElement(i).computeForce()));
			}
			System.out.print("\n");
			// Printing Reaction Forces
			System.out.println(" Reaction Forces ");
			System.out.printf("node%14s%15s%15s", "r1", "r2", "r3");
			System.out.print("\n");
			for (int i = 0; i < this.getNumberOfNodes(); i++) {
				System.out.printf("%5d", i);
				System.out.println(ArrayFormat.format(this.getNode(i)
							.getReactionForces()));
				}
			System.out.print("\n");
		}
			
		if (this.state == "Therm" || this.state == "ThermoMech") {
			System.out.println(" Temperatures ");
			System.out.printf(" node%14s", "temperature");
			System.out.print("\n");
			for (int i = 0; i < this.getNumberOfNodes(); i++) {
				System.out.printf("%5d%15f ", i, this.getNode(i).getConstraint().getNodalTemp());
				System.out.print("\n");
			}
			System.out.print("\n");
			System.out.println(" Reaction Fluxes ");
			System.out.printf(" node%12s", "q");
			System.out.print("\n");
			for (int i = 0; i < this.getNumberOfNodes(); i++) {
				System.out.printf("%5d%15f ", i, this.getNode(i).getFlux().getNodalFlux());
				System.out.print("\n");
			}
			/*System.out.print("\n");
			for (int i = 0; i < this.getNumberOfNodes(); i++) {
				System.out.printf("%5d", i);
				System.out.print(ArrayFormat.format(this.getNode(i).getDofNumbers()));
				System.out.print("\n");
			}*/
		}
	}
}
