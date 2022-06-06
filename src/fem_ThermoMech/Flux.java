package fem_ThermoMech;

import inf.text.ArrayFormat;

public class Flux {
	private double flux = 0; // Default value is 0
	public Flux ( double fluxNodal ) {
		// Setting components
		this.flux = fluxNodal;
	}
	// To get the Nodal Flux
	public double getNodalFlux () {
		return this . flux ;
	}
	// To print the Flux
	public void print () {
		System . out . println ( ArrayFormat . format ( this . flux ));
	}
}
