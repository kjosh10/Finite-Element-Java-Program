package fem_ThermoMech;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import iceb.jnumerics.Vector3D;
import inf.v3d.obj.Arrow;
import inf.v3d.obj.Cone;
import inf.v3d.obj.Cylinder;
import inf.v3d.obj.CylinderSet;
import inf.v3d.obj.PolygonSet;
import inf.v3d.obj.Sphere;
import inf.v3d.obj.Text;
import inf.v3d.view.Viewer;

public class Visualizer {
	private double displacementScale = 1e4;
	private double forceSymbolScale = 4e-5, forceSymbolRadius = 0.075, 
			constraintSymbolScale = 1;
	private double elementNormalForcesScale = 0.00003;
	private double nonHomDBCScale = 4e-5;
	private Structure struct;
	private Viewer v;
	private double[] elementNormalForce;

	public Visualizer(Structure struct, Viewer v, String state) {
		// TODO Auto-generated constructor stub
		this.struct = struct;
		this.v = v;
		if (state == "Mech" || state == "ThermoMech") {
			this.setElementNormalForce();
		}
	}
	
	// To set element Normal forces in the vector
	public void setElementNormalForce() {
		this.elementNormalForce = new double[this.struct.getNumberOfElements()];
		for (int i = 0; i <this.struct.getNumberOfElements(); i++) {
			this.elementNormalForce[i] = this.struct.getElement(i).computeForce();
		}
	}

	// To draw Elements
	public void drawElements() {
		CylinderSet elem_visualizer = new CylinderSet();
		for (int i = 0; i < this.struct.getNumberOfElements(); i++) {
			Element elem = this.struct.getElement(i);
			double area = elem.getArea();
			double r = Math.sqrt(4 * area / Math.PI); 
			// 4 Multiplied for scaling
			elem_visualizer.addCylinder(this.struct.getElement(i).getNode1().getPosition().toArray(),
					this.struct.getElement(i).getNode2().getPosition().toArray(), r);
			v.addObject3D(elem_visualizer);
		}
	}

	// To draw the constraints
	public void drawConstraints() {
		double height = this.constraintSymbolScale; // this.symbolScaleConstraint;
		for (int i = 0; i < this.struct.getNumberOfNodes(); i++) {
			if (this.struct.getNode(i).getConstraint() != null) {
				for (int j = 0; j < this.struct.getNode(i).getDofNumbers().length; j++) {
					Cone coneConstraint = new Cone();
					coneConstraint.setHeight(height);
					if (this.struct.getNode(i).getConstraint().isFree(j) == false) {

						if (j == 0) { // to switch directions
							coneConstraint.setCenter(this.struct.getNode(i).getPosition().getX1() - height,
									this.struct.getNode(i).getPosition().getX2(),
									this.struct.getNode(i).getPosition().getX3());
							coneConstraint.setDirection(1, 0, 0);
							coneConstraint.setColor("black");
							v.addObject3D(coneConstraint);
						} else if (j == 1) {
							coneConstraint.setCenter(this.struct.getNode(i).getPosition().getX1(),
									this.struct.getNode(i).getPosition().getX2() - height,
									this.struct.getNode(i).getPosition().getX3());
							coneConstraint.setDirection(0, 1, 0);
							coneConstraint.setColor("black");
							v.addObject3D(coneConstraint);
						} else if (j == 2) {
							coneConstraint.setCenter(this.struct.getNode(i).getPosition().getX1(),
									this.struct.getNode(i).getPosition().getX2(),
									this.struct.getNode(i).getPosition().getX3() - height);
							coneConstraint.setDirection(0, 0, 1);
							coneConstraint.setColor("black");
							v.addObject3D(coneConstraint);
						}
					}

				}
			}
		}

	}

	// To draw Element Forces
	public void drawElementForces() {
		for (int i = 0; i < this.struct.getNumberOfElements(); i++) {
			// For Node 1
			if (this.struct.getElement(i).getNode1().getForce() != null) {
				for (int j = 0; j < this.struct.getElement(i).getNode1().getDofNumbers().length; j++) {
					Arrow forceArrow = new Arrow();
					forceArrow.setRadius(this.forceSymbolRadius);
					forceArrow.setPoint1(this.struct.getElement(i).getNode1().getPosition().toArray());
					double[] point2 = new double[3];
					if (this.struct.getElement(i).getNode1().getForce().getComponent(j) != 0) {
						if (j == 0) { // to switch directions
							point2[0] = this.struct.getElement(i).getNode1().getPosition().getX1()
									+ this.struct.getElement(i).getNode1().getForce().getComponent(j)
											* this.forceSymbolScale;
							point2[1] = this.struct.getElement(i).getNode1().getPosition().getX2();
							point2[2] = this.struct.getElement(i).getNode1().getPosition().getX3();
						} else if (j == 1) {
							point2[0] = this.struct.getElement(i).getNode1().getPosition().getX1();
							point2[1] = this.struct.getElement(i).getNode1().getPosition().getX2()
									+ this.struct.getElement(i).getNode1().getForce().getComponent(j)
											* this.forceSymbolScale;
							point2[2] = this.struct.getElement(i).getNode1().getPosition().getX3();
						} else if (j == 2) {
							point2[0] = this.struct.getElement(i).getNode1().getPosition().getX1();
							point2[1] = this.struct.getElement(i).getNode1().getPosition().getX2();
							point2[2] = this.struct.getElement(i).getNode1().getPosition().getX3()
									+ this.struct.getElement(i).getNode1().getForce().getComponent(j)
											* this.forceSymbolScale;
						}
						forceArrow.setPoint2(point2);
						forceArrow.setColor("red");
						v.addObject3D(forceArrow);
					}
				}

			}

			// For Node 2
			if (this.struct.getElement(i).getNode2().getForce() != null) {
				for (int j = 0; j < this.struct.getElement(i).getNode1().getDofNumbers().length; j++) {
					Arrow forceArrow = new Arrow();
					forceArrow.setRadius(this.forceSymbolRadius);
					forceArrow.setPoint1(this.struct.getElement(i).getNode2().getPosition().toArray());
					double[] point2 = new double[3];
					if (this.struct.getElement(i).getNode2().getForce().getComponent(j) != 0) {
						if (j == 0) { // to switch directions
							point2[0] = this.struct.getElement(i).getNode2().getPosition().getX1()
									+ this.struct.getElement(i).getNode2().getForce().getComponent(j)
											* this.forceSymbolScale;
							point2[1] = this.struct.getElement(i).getNode2().getPosition().getX2();
							point2[2] = this.struct.getElement(i).getNode2().getPosition().getX3();
						} else if (j == 1) {
							point2[0] = this.struct.getElement(i).getNode2().getPosition().getX1();
							point2[1] = this.struct.getElement(i).getNode2().getPosition().getX2()
									+ this.struct.getElement(i).getNode2().getForce().getComponent(j)
											* this.forceSymbolScale;
							point2[2] = this.struct.getElement(i).getNode2().getPosition().getX3();
						} else if (j == 2) {
							point2[0] = this.struct.getElement(i).getNode2().getPosition().getX1();
							point2[1] = this.struct.getElement(i).getNode2().getPosition().getX2();
							point2[2] = this.struct.getElement(i).getNode2().getPosition().getX3()
									+ this.struct.getElement(i).getNode2().getForce().getComponent(j)
											* this.forceSymbolScale;
						}
						forceArrow.setPoint2(point2);
						forceArrow.setColor("red");
						v.addObject3D(forceArrow);
					}
				}

			}

		}
	}

	// To draw Displacements
	public void drawDisplacements() {
		for (int i = 0; i < this.struct.getNumberOfElements(); i++) {
			Node n1_deform = new Node(
					this.struct.getElement(i).getNode1().getPosition().c1
							+ this.struct.getElement(i).getNode1().getDisplacement().c1 * this.displacementScale,
					this.struct.getElement(i).getNode1().getPosition().c2
							+ this.struct.getElement(i).getNode1().getDisplacement().c2 * this.displacementScale,
					this.struct.getElement(i).getNode1().getPosition().c3
							+ this.struct.getElement(i).getNode1().getDisplacement().c3 * this.displacementScale);
			Node n2_deform = new Node(
					this.struct.getElement(i).getNode2().getPosition().c1
							+ this.struct.getElement(i).getNode2().getDisplacement().c1 * this.displacementScale,
					this.struct.getElement(i).getNode2().getPosition().c2
							+ this.struct.getElement(i).getNode2().getDisplacement().c2 * this.displacementScale,
					this.struct.getElement(i).getNode2().getPosition().c3
							+ this.struct.getElement(i).getNode2().getDisplacement().c3 * this.displacementScale);
			CylinderSet display_visualizer = new CylinderSet();
			/* lambda and conductivity is 0 here because they are not 
				used for visualizing displacements*/
			Element elem_deform = new Element(this.struct.getElement(i).getEModulus(),
					this.struct.getElement(i).getArea(), 0, 0, n1_deform, n2_deform); // this.struct.getElement(i);
			double area = elem_deform.getArea();
			double r = Math.sqrt(area / Math.PI);
			display_visualizer.addCylinder(elem_deform.getNode1().getPosition().toArray(),
					elem_deform.getNode2().getPosition().toArray(), r);
			display_visualizer.setColor("blue");
			v.addObject3D(display_visualizer);

		}
		
	}

	// To draw temperatures
	public void drawTemperatures() {
		// Before drawing the elements, we set some values for coloring
	    // Create an array of all Nodal Temperatures to know the maximum and minimum
	    ArrayList<Double> nodalTemp = new ArrayList<Double>();
	    for (int i = 0; i < this.struct.getNumberOfNodes(); i++) { 
			// Loop over all elements
			nodalTemp.add(this.struct.getNode(i).getConstraint().getNodalTemp());
		}
		// Knowing the max and min, we create a scale array of 4 intervals
		double max = Collections.max(nodalTemp);
		double min = Collections.min(nodalTemp);
		double step = (max - min) / 4;
		double[] colorIntervals = { min, min + step, min + 2 * step, min + 3 * step, max };
		Sphere temperature;
		
		for (int i = 0; i < this.struct.getNumberOfElements(); i++) {
			double area = this.struct.getElement(i).getArea();
			double r = Math.sqrt(5 * 4 * area / Math.PI);
			// 4 and 5 Multiplied for scaling
			// For each node create a sphere object
			temperature = new Sphere(this.struct.getElement(i).getNode1().getPosition().getX1(),
					this.struct.getElement(i).getNode1().getPosition().getX2(),
					this.struct.getElement(i).getNode1().getPosition().getX3());
			temperature.setRadius(r);
	        
			double nodalTemp1 = this.struct.getElement(i).getNode1()
					.getConstraint().getNodalTemp(); 
			// Get the nodal Temperature of the Node 1 of the current element

	        /* Check in which interval is it located to set the color
	         these intervals are based on the RGB color convention */
	        // RGB (R, G, B) (each value from 1 to 0, where 1 means 255)
	          if (nodalTemp1 < colorIntervals[1]) {
	            temperature.setColor(0, Math.abs(nodalTemp1 - colorIntervals[0]) / step, 1);
	            // in this interval, the color changes between (0, 0, 255) to (0, 255, 255)
	            // [blue to cyan]
	          } else if (nodalTemp1 < colorIntervals[2]) {
	        	  temperature.setColor(0, 1, 1 - Math.abs(nodalTemp1 - colorIntervals[1]) / step);
	            // (0, 255, 255) to (0, 255, 0)
	          } else if (nodalTemp1 < colorIntervals[3]) {
	        	  temperature.setColor(Math.abs(nodalTemp1 - colorIntervals[2]) / step, 1, 0);
	            // (0, 255, 0) to (255, 255, 0)
	          } else {
	        	  temperature.setColor(1, 1 - Math.abs(nodalTemp1 - colorIntervals[3]) / step, 0);
	            // (255, 255, 0) to (255, 0, 0) yellow to red
	          }
	          temperature.setOpacity(1);
	          v.addObject3D(temperature);
	          
	          // For each node create a sphere object
	          temperature = new Sphere(this.struct.getElement(i).getNode2().getPosition().getX1(),
	        		  this.struct.getElement(i).getNode2().getPosition().getX2(),
						this.struct.getElement(i).getNode2().getPosition().getX3());
	          temperature.setRadius(r*1.15);
	          double nodalTemp2 = this.struct.getElement(i).getNode2()
	        		  .getConstraint().getNodalTemp(); 
	          // Get the nodal Temperature of the Node 2 of the current element

		      // Check in which interval is it located to set the color
		      // these intervals are based on the RGB color convention
		      // RGB (R, G, B) (each value from 1 to 0, where 1 means 255)
		      if (nodalTemp2 < colorIntervals[1]) {
		          temperature.setColor(0, Math.abs(nodalTemp2 - colorIntervals[0]) / step, 1);
		          // in this interval, the color changes between (0, 0, 255) to (0, 255, 255)
		          // [blue to cyan]
		      } else if (nodalTemp2 < colorIntervals[2]) {
		        	  temperature.setColor(0, 1, 1 - Math.abs(nodalTemp2 - colorIntervals[1]) / step);
		            // (0, 255, 255) to (0, 255, 0)
		      } else if (nodalTemp2 < colorIntervals[3]) {
		        	  temperature.setColor(Math.abs(nodalTemp2 - colorIntervals[2]) / step, 1, 0);
		            // (0, 255, 0) to (255, 255, 0)
		      } else {
		        	  temperature.setColor(1, 1 - Math.abs(nodalTemp2 - colorIntervals[3]) / step, 0);
		            // (255, 255, 0) to (255, 0, 0) yellow to red
		      }
		      temperature.setOpacity(1);
		      v.addObject3D(temperature); 
		}
		
	}
	
	// To draw Internal Element Normal Forces
	public void drawElementNormalForces() {
		double alphaScale = this.elementNormalForcesScale;
		Vector3D n1 = new Vector3D(0, 0, 1); // Unit Vector along Z-axis
		Vector3D n2 = new Vector3D(0, 1, 0); // Unit Vector along Y-axis
		Vector3D p, d, s1, s2;
		PolygonSet elementNormalForce = new PolygonSet();
		for (int i = 0; i < this.struct.getNumberOfElements(); i++) {
			d = this.struct.getElement(i).getE1(); // Unit vector along element
			if (n1.vectorProduct(d).normTwo() != 0) {  // <= 1e-2
				p = n1.vectorProduct(d);
				s1 = this.struct.getElement(i).getNode1().getPosition()
						.add(p.multiply(alphaScale* this.struct.getElement(i).computeForce())); //* this.struct.getElement(i).computeForce()
				s2 = this.struct.getElement(i).getNode2().getPosition()
						.add(p.multiply(alphaScale* this.struct.getElement(i).computeForce()));
			} else {
				p = n2.vectorProduct(d);
				s1 = this.struct.getElement(i).getNode1().getPosition()
						.add(p.multiply(alphaScale* this.struct.getElement(i).computeForce()));
				s2 = this.struct.getElement(i).getNode2().getPosition()
						.add(p.multiply(alphaScale* this.struct.getElement(i).computeForce()));
			}
			
			elementNormalForce.insertVertex(this.struct.getElement(i).getNode1().getPosition().toArray(), this.struct.getElement(i).computeForce());
			elementNormalForce.insertVertex(this.struct.getElement(i).getNode2().getPosition().toArray(), this.struct.getElement(i).computeForce());
			elementNormalForce.insertVertex(s2.toArray(), this.struct.getElement(i).computeForce());
			elementNormalForce.insertVertex(s1.toArray(), this.struct.getElement(i).computeForce());
			elementNormalForce.polygonComplete();
			
			elementNormalForce.setVisible(true);
			elementNormalForce.setColoringByData(true);
			elementNormalForce.setOutlinesVisible(true);
			elementNormalForce.setContourLinesVisible(false);
			elementNormalForce.createColors();
			elementNormalForce.setOpacity(1);
			
			
			v.addObject3D(elementNormalForce);

		}

	}
	
	// To Draw Element Normal Forces in the elements itself
	public void drawElementInternalForceColored() {
		// Before drawing the elements, we set some values for coloring
		double absMax = this.maxAbsNormalForce();
		double step = this.maxAbsNormalForce() * 2 / 4;
		double[] colorIntervals = { -absMax, -absMax + step, -absMax + 2 * step,
				-absMax + 3 * step, -absMax + 4 * step};
		Cylinder elem_visualizer;
		for (int i = 0; i < this.struct.getNumberOfElements(); i++) {
			Element elem = this.struct.getElement(i);
			double area = elem.getArea();
			double r = Math.sqrt(4 * area / Math.PI);
			// 4 Multiplied for scaling
			elem_visualizer = new Cylinder(this.struct.getElement(i).getNode1().getPosition().toArray(),
					this.struct.getElement(i).getNode2().getPosition().toArray());
			elem_visualizer.setRadius(r);
			
			// Coloring Elements
			double internalForce = this.struct.getElement(i).computeForce(); 
			// Get the internal force of the current element
			// Check in which interval is it located to set the color
			// these intervals are based on the RGB color convention
			// RGB (R, G, B) (each value from 1 to 0, where 1 means 255)
			if (internalForce == 0) {
				// Green for elements with no internal force
				elem_visualizer.setColor("green");
			} else if (internalForce < colorIntervals[1]) {
				elem_visualizer.setColor(0, Math.abs(internalForce - colorIntervals[0]) / step, 1);
				// in this interval, the color changes between (0, 0, 255) to (0, 255, 255)
				// [blue to cyan]
			} else if (internalForce < colorIntervals[2]) {
				elem_visualizer.setColor(0, 1, 1 - Math.abs(internalForce - colorIntervals[1]) / step);
				// (0, 255, 255) to (0, 255, 0)
			} else if (internalForce < colorIntervals[3]) {
				elem_visualizer.setColor(Math.abs(internalForce - colorIntervals[2]) / step, 1, 0);
				// (0, 255, 0) to (255, 255, 0)
			} else {
				elem_visualizer.setColor(1, 1 - Math.abs(internalForce - colorIntervals[3]) / step, 0);
				// (255, 255, 0) to (255, 0, 0) yellow to red
			}
			v.addObject3D(elem_visualizer);
		}
	}

	// To draw Non Homogeneous Prescribed Displacement
	public void drawNonHomDBC() {
		boolean nonHomFlag = false;
		for (int i = 0; i < this.struct.getNumberOfNodes(); i++) {
			if (this.struct.getNode(i).getNonHomDBC() != null) {
				nonHomFlag = true; // Their is atleast one NHDBC
				break;
			}
		}
		if (nonHomFlag) { // Draw Non-Hom DBC
			for (int i = 0; i < this.struct.getNumberOfNodes(); i++) {
				if (this.struct.getNode(i).getNonHomDBC() != null) {
					double[] point2 = new double[3];
					for (int j = 0; j < this.struct.getNode(i).getNonHomDBC().length;
							j++) {
						Arrow nonHomDBCArrow = new Arrow();
						if (this.struct.getNode(i).getNonHomDBC()[j] != 0) {
							double[] prescDispVector = new double[3];
							prescDispVector[j] = this.struct.getNode(i).getNonHomDBC()[j];
							point2[0] = this.struct.getNode(i).getPosition().getX1()
									+ prescDispVector[0]* this.nonHomDBCScale;
            				point2[1] = this.struct.getNode(i).getPosition().getX2()
									+ prescDispVector[1]* this.nonHomDBCScale;
							point2[2] = this.struct.getNode(i).getPosition().getX3()
									+ prescDispVector[2]* this.nonHomDBCScale;
						}
						nonHomDBCArrow.setRadius(this.forceSymbolRadius);
						// Same Radius used as in Force Symbol
						nonHomDBCArrow.setPoint1(this.struct.getNode(i).getPosition().toArray());
						nonHomDBCArrow.setPoint2(point2);
						nonHomDBCArrow.setColor("blue");
						v.addObject3D(nonHomDBCArrow);			
					}
					
				}
			}
		}
	}
	
	// To draw Legend for Element Normal Forces 
	public void drawLegendElemForces() {
		double force;
		Viewer legendViewer = new Viewer(); // create a new viewer for viewing the legend
		
		PolygonSet boxColor = new PolygonSet(); // use this object to show respective 
		// colors in a box and stack them one above the other
		
		for (int i = 0; i < 11; i++) { // here 10 is the number of colors which 
									   // we want to show
			// force for current iteration
			
			force = -this.maxAbsNormalForce() + i*
					2*this.maxAbsNormalForce() / 10;
			// insert the box for the current iteration
			boxColor.insertVertex(0, i, 0, force);
			boxColor.insertVertex(1, i, 0, force);
			boxColor.insertVertex(1, i+1, 0, force);
			boxColor.insertVertex(0, i+1, 0, force);
			boxColor.polygonComplete();
			legendViewer.addObject3D(boxColor);
			// create a text field to show the values of the color in the box
			Text text = new Text(1.5, 0.5 + i, 0, 0.2, String.valueOf(force));
			text.setColor("black");
			legendViewer.addObject3D(text);
		}
		boxColor.setVisible(true);
		boxColor.setColoringByData(true);
		boxColor.setOutlinesVisible(true);
		boxColor.setContourLinesVisible(false);
		boxColor.createColors();
		boxColor.setOpacity(1);
		legendViewer.setVisible(true);

	}
	
	// To draw Legend for Element Normal Forces with coloring in the elements 
	public void drawLegendElemForcesColored() {
		double force;
		Viewer legendViewer = new Viewer(); // create a new viewer for viewing the legend
		
		PolygonSet boxColor = new PolygonSet(); // use this object to show respective 
		// colors in a box and stack them one above the other
		
		for (int i = 0; i < 6; i++) { // here 5 is the number of colors which 
									   // we want to show and one Caption
			// force for current iteration
			
			if (i < 5) {
				force = -this.maxAbsNormalForce() + i*
						2*this.maxAbsNormalForce() / 4;
				
				// insert the box for the current iteration
				boxColor.insertVertex(0, i, 0, force);
				boxColor.insertVertex(1, i, 0, force);
				boxColor.insertVertex(1, i+1, 0, force);
				boxColor.insertVertex(0, i+1, 0, force);
				boxColor.polygonComplete();
				legendViewer.addObject3D(boxColor);
				// create a text field to show the values of the color in the box
				Text text = new Text(1.5, 0.5 + i, 0, 0.2, String.valueOf(force));
				text.setColor("black");
				legendViewer.addObject3D(text);
			} else {
				// create a text field to show the Caption
				Text textCaption = new Text(1.5, 0.5 + i, 0, 0.2, "Internal Force Legend");
				textCaption.setColor("black");
				legendViewer.addObject3D(textCaption);
			}
		}
		boxColor.setVisible(true);
		boxColor.setColoringByData(true);
		boxColor.setOutlinesVisible(true);
		boxColor.setContourLinesVisible(false);
		boxColor.createColors();
		boxColor.setOpacity(1);
		legendViewer.setVisible(true);

	}
	
	// To draw legend for Nodal Temperatures
	public void drawLegendNodalTemp() {
		double nodalTemp;
		Viewer legendViewer = new Viewer(); // create a new viewer for viewing the legend
		
		ArrayList<Double> nodalTemperatures = new ArrayList<Double>();
		for (int i = 0; i < this.struct.getNumberOfNodes(); i++) { 
			// Loop over all elements
			nodalTemperatures.add(this.struct.getNode(i).getConstraint().getNodalTemp());
		}
		// Knowing the max and min, we create a scale array of 4 intervals
		double max = Collections.max(nodalTemperatures);
		double min = Collections.min(nodalTemperatures);
	
		PolygonSet boxColor = new PolygonSet(); // use this object to show respective 
		// colors in a box and stack them one above the other
		
		for (int i = 0; i < 6; i++) { // here 5 is the number of colors which 
									   // we want to show and a Caption
			if (i < 5) {
				// nodalTemp for current iteration
				nodalTemp = min + i * (max - min) / 4;
				// insert the box for the current iteration
				boxColor.insertVertex(0, i, 0, nodalTemp);
				boxColor.insertVertex(1, i, 0, nodalTemp);
				boxColor.insertVertex(1, i+1, 0, nodalTemp);
				boxColor.insertVertex(0, i+1, 0, nodalTemp);
				boxColor.polygonComplete();
				legendViewer.addObject3D(boxColor);
				// create a text field to show the values of the color in the box
				Text text = new Text(1.5, 0.5 + i, 0, 0.2, String.valueOf(nodalTemp));
				text.setColor("black");
				legendViewer.addObject3D(text);
			} else {
				// create a text field to show the Caption
				Text textCaption = new Text(1.5, 0.5 + i, 0, 0.2, "Nodal Temperature Legend");
				textCaption.setColor("black");
				legendViewer.addObject3D(textCaption);
			}
			
		}
		boxColor.setVisible(true);
		boxColor.setColoringByData(true);
		boxColor.setOutlinesVisible(true);
		boxColor.setContourLinesVisible(false);
		boxColor.createColors();
		boxColor.setOpacity(1);
		legendViewer.setVisible(true);

	}
	
	// To calculate Minimum and Maximum Internal Normal Force
	public double[] minMaxNormalForce() { 
		double[] elementNormalForcesCopy = 
				Arrays.copyOf(this.elementNormalForce, this.elementNormalForce.length);
		Arrays.sort(elementNormalForcesCopy);
		return new double[] {elementNormalForcesCopy[0], 
				elementNormalForcesCopy[this.struct.getNumberOfElements() - 1]};
	}
	
	// To calculate maximum absolute normal internal force
	public double maxAbsNormalForce() {
		double[] elementNormalForcesCopy = 
				Arrays.copyOf(this.elementNormalForce, this.elementNormalForce.length);
		Arrays.sort(elementNormalForcesCopy);
		if (Math.abs(elementNormalForcesCopy[0]) > 
			Math.abs(elementNormalForcesCopy[this.struct.getNumberOfElements() - 1])){
			return Math.abs(elementNormalForcesCopy[0]);
		} else {
			return Math.abs(elementNormalForcesCopy[this.struct.getNumberOfElements() - 1]);
		}
	}
	
	// To set Constraint symbol scale
	public void setConstraintSymbolScale(double scale) {
		this.constraintSymbolScale = scale;
	}

	// To set Displacement scale
	public void setDisplacementScale(double scale) {
		this.displacementScale = scale;
	}

	// To set Force symbol scale
	public void setForceSymbolScale(double scale) {
		this.forceSymbolScale = scale;
	}

	// To set Force symbol arrow radius 
	public void setForceSymbolRadius(double radius) {
		this.forceSymbolRadius = radius;
	}
	
	// To set Element Normal Forces scale
	public void setElementNormalForcesScale(double scale) {
		this.elementNormalForcesScale = scale;
	}
	
	// To set non homogeneous arrow symbol scale
	public void setNonHomDBCScale(double scale) {
		this.nonHomDBCScale = scale;
	}
}
