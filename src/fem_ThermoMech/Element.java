package fem_ThermoMech;


import iceb.jnumerics.Array2DMatrix;
import iceb.jnumerics.BLAM;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.Vector3D;
import inf.text.ArrayFormat;

public class Element {
	private double area;
	private Node n1, n2;
	private double eModulus, lambda, thermalConductivity;
	private double[] elements_param = new double[5];
	private int[] dofNumbers = new int[6];
	
	// Constructor for Thermal and Thermo-Mechanical Problem
	public Element(double e, double a, double lamda_parameter, double thermalConductivity,
			Node n1, Node n2) {
		// TODO Auto-generated constructor stub
		// Set parameters in the constructor
		this . eModulus = e;
		this . lambda = lamda_parameter;
		this . area = a;
		this . thermalConductivity = thermalConductivity;
		this . n1 = n1;
		this . n2 = n2;
		elements_param[0] = e;
		elements_param[1] = a;
		elements_param[2] = lamda_parameter;
		elements_param[3] = thermalConductivity;
		elements_param[4] = this.getLength();
	}
	
	// Constructor for Mechanical Problems
	public Element(double e, double a, Node n1, Node n2) {
		// TODO Auto-generated constructor stub
		// Set parameters in the constructor
		this . eModulus = e;
		this . lambda = 0; // Default value is zero
		this . area = a;
		this . thermalConductivity = 0; // Default value is zero
		this . n1 = n1;
		this . n2 = n2;
		elements_param[0] = e;
		elements_param[1] = a;
		elements_param[2] = lambda;
		elements_param[3] = thermalConductivity;
		elements_param[4] = this.getLength();
	}
	
	// Compute stiffness matrix for a particular element
	public IMatrix computeStiffnessMatrix() {
		Array2DMatrix kLocalElement = new Array2DMatrix(2, 2, 1, -1, -1, 1);
		kLocalElement = (Array2DMatrix) kLocalElement.multiply(this.eModulus*this.area/this.getLength());
		
		IMatrix T = new Array2DMatrix(2, 6);
		IMatrix tmp = new Array2DMatrix(6, 2);
		IMatrix kGlobalElement = new Array2DMatrix(6, 6);
		// Set T matrix
		T . set(0, 0, this . getE1() . c1);
		T . set(0, 1, this . getE1() . c2);
		T . set(0, 2, this . getE1() . c3);
		T . set(1, 3, this . getE1() . c1);
		T . set(1, 4, this . getE1() . c2);
		T . set(1, 5, this . getE1() . c3);
		
		//System.out.println(MatrixFormat.format(kLocalElement));
		//System.out.println(MatrixFormat.format(T));
		// Calculate T'*kLocalElement*T
		BLAM.multiply(1.0, BLAM.TRANSPOSE, T, BLAM.NO_TRANSPOSE, kLocalElement, 0.0, tmp);
		BLAM.multiply(1.0, BLAM.NO_TRANSPOSE, tmp, BLAM.NO_TRANSPOSE, T, 0.0, kGlobalElement);
		return kGlobalElement;
	}
	
	// Compute Element wise Temperature Stiffness Matrix
	public IMatrix computeTempStiffnessMatrix() {
		Array2DMatrix kTetaTetaElement = new Array2DMatrix(2, 2, 1, -1, -1, 1);
		kTetaTetaElement = (Array2DMatrix) kTetaTetaElement.multiply(this.lambda*this.area/this.getLength());
		return kTetaTetaElement;
	}
	
	// Compute combined stiffness Matrix K_UTeta
	public IMatrix computeDispTempStiffnessMatrix() {
		IMatrix T = new Array2DMatrix(6, 2);
		IMatrix kUTetaElement = new Array2DMatrix(6, 2);
		// Set T matrix
		T . set(0, 0,  this.getE1().c1); T . set(0, 1,  this.getE1().c1);
		T . set(1, 0,  this.getE1().c2); T . set(1, 1,  this.getE1().c2);
		T . set(2, 0,  this.getE1().c3); T . set(2, 1,  this.getE1().c3);
		T . set(3, 0, -this.getE1().c1); T . set(3, 1, -this.getE1().c1);
		T . set(4, 0, -this.getE1().c2); T . set(4, 1, -this.getE1().c2);
		T . set(5, 0, -this.getE1().c3); T . set(5, 1, -this.getE1().c3);
		
		kUTetaElement = (IMatrix) T.multiply(this.eModulus*this.area*
				this.thermalConductivity/2);
		//System.out.println(MatrixFormat.format(kUTetaElement));
		return kUTetaElement;
	}
	
	// Enumerate elements
	public void enumerateDOFs() {
		for (int i = 0; i < this.dofNumbers.length; i++) {
			if (i < 3) {
				dofNumbers[i] = this.n1.getDofNumbers()[i];
			} else {
				dofNumbers[i] = this.n2.getDofNumbers()[i - 3];
			}
		}
	}
	
	// Get element enumerated DoFs
	public int[] getDofNumbers() {
		return this.dofNumbers;
	}
	
	// To compute the element internal force
	public double computeForce() {
		Vector3D disp_vect = this.getE1();
		double u1 = this.getNode1().getDisplacement().scalarProduct(disp_vect);
		double u2 = this.getNode2().getDisplacement().scalarProduct(disp_vect);
		double force = (this.eModulus*this.getArea()/this.getLength())*(u2-u1);
		return force;
	}
	
	// To compute unit vector along the direction of the element
	public Vector3D getE1() {
		return (this.n2.getPosition().subtract(this.n1.getPosition())).multiply(1/this.getLength());
	}

	// To compute the length of the element
	public double getLength() {
		Vector3D n1_pos = n1.getPosition();
		Vector3D n2_pos = n2.getPosition();
		return n1_pos.subtract(n2_pos).normTwo();
	}

	// To get node 1 of the element
	public Node getNode1() {
		return this.n1;
	}

	// To get node 2 of the element
	public Node getNode2() {
		return this.n2;
	}

	// To get the Area of the element
	public double getArea() {
		return this.area;
	}

	// To get the Young's Modulus of the element
	public double getEModulus() {
		return this.eModulus;
	}

	// To print the element
	public void print() {
		System.out.println(ArrayFormat.format(this.elements_param));
	}
}
