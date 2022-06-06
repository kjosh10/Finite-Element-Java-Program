package fem_ThermoMech;

import inf.text.ArrayFormat;

public class Force {
	private double [] components = new double [3]; 
	// Array for prescribing components of Force
	public Force ( double r1 , double r2 , double r3 ) {
		// Setting components
		this . components [0] = r1 ;
		this . components [1] = r2 ;
		this . components [2] = r3 ;
	}
	// To get Force component c+1, e.g. when c = 1, we get F_y component
	public double getComponent ( int c ) {
		return this . components [c ];
	}
	// To print the Force
	public void print () {
		System . out . println ( ArrayFormat . format ( this . components ));
	}
}
