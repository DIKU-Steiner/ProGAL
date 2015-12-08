package ProGAL.geom2d;

import java.awt.Color;

import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.math.Functions;

/* 
 * A hyperbola is the set of all points P = (x,y) for which the absolute value of the difference of the distances to 
 * two focal points F_1 and F_2 is a constant 
 */
public class Hyperbola {
	Vector tVector;
	double alpha;
	Point F1;        // focal point
	Point F2;        // focal point
	
	public Hyperbola(Point F1, Point F2) {
		Point midPoint = Point.midPoint(F1,  F2);
		tVector = new Vector(midPoint);
		this.F1 = F1.subtract(tVector);
		this.F2 = F2.subtract(tVector);
		alpha = this.F1.polarAngle(); System.out.println(Functions.toDeg(alpha));
		this.F1.rotation(-alpha);
		this.F2.rotation(-alpha);
	}
	/** returns the pair o points on the two branches that are closest to each other */
	public Point[] getVertices() {
		return null;  // to be implemented
	}
	
	/** returns the smallest radius of curvature (at the vertices) */
	public double getSmallestRadiusOfCurvature() {
		return -1;   // to be implemented;
	}
	
	/** returns the transverse axis (through the foci and vertices) */
	public LineSegment getTransverseAxis() {
		return null;   // to be implemented
	}
	
	/** returns the semi-major axis (perpendicular to the transverse axis) */
	public LineSegment getSemiMajorAxis() {
		return null; // to be implemented
	}
	
	/** returns the center of the hyperbola */
	public Point getCenter() {
		return null;   // to be implemented
	}
	
	/** returns the two asymptotes (intersecting in the center) */
	public Line[] getAsymptotes() {
		return null; // to be implemented
	}
	
	/** translates and rotates the hyperbola so that its center is 
	 * in origo and its transverse axis is the x-axis */
	public void alignXAxis() {
		// to be implemented
	}
	
	public void alignYAxis() {
		// to be implemented
	}
	
	public static void main(String[] args) {
		J2DScene scene  = J2DScene.createJ2DSceneInFrame();
		Point R1 = new Point(1, 0); R1.toScene(scene, 0.1, Color.blue);
		Point R2 = new Point(2, 4); R2.toScene(scene, 0.1, Color.blue);
		
		Hyperbola h = new Hyperbola(R1, R2);
		h.F1.toScene(scene, 0.1, Color.red);
		h.F2.toScene(scene, 0.1, Color.red);
	}
}
