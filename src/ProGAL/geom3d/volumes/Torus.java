package ProGAL.geom3d.volumes;

import java.awt.Color;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;

public class Torus  {
	protected Point center;
	protected Vector normal;
	protected double R;
	protected double r;


	public Torus(Point center, Vector normal, double majorRadius, double minorRadius) {
		this.center = center;
		this.normal = normal;
		this.R = majorRadius;
		this.r = minorRadius;
	}
	
	public Circle getToroidalCircle() {
		return new Circle(center, R, normal);
	}
	
	public Circle getPoloidalCircle() {
		return null;
	}
	
	/** Two circles of major radius R, tilted to the center plane
		by slopes of r/R and -r/R and offset from the center by minor radius distance r */
	public Circle[] getVillarceauCircles () {
		return null;
	}
	
	
	public void toScene(J3DScene scene, Color clr) {
		int iMax = 200;
		double delta = 2*Math.PI/iMax;
		
		Vector v = normal.getOrthonormal();
		Point p = center.add(v.scaleToLength(R));
		for (int i = 0; i < iMax; i++) {
			new Circle(p, r, v.cross(normal)).toScene(scene, 0.01, 32);
			p.rotation(normal, delta, center);
		}
				
	}
	
	public static void main(String[] args) {
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		Point center = new Point(0, 0, 0);
		Vector normal = new Vector(0, 0, 1).normalize();
		double R = 0.5;
		double r = 0.1;
		Torus torus = new Torus(center, normal, R, r);
		torus.toScene(scene, Color.blue);
	}
}
