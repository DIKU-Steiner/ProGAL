package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;

/**
 * An infinite cone is the union of all half-lines starting at an <em>apex</em>-point whose <em>angle</em> to a specified 
 * <em>axis</em>-vector is less than a specified constant. This class represents an infinitely extended cone 
 * defined by these three elements.
 * @author R.Fonseca
 */
public class InfCone {
	protected Point apex;
	protected Vector axis;
	protected double angle;
	
	/**
	 * Construct the cone by specifying the apex-point, direction-vector and angle.
	 * @param apex
	 * @param axis
	 * @param angle
	 */
	public InfCone(Point apex, Vector axis, double angle){
		this.apex = apex;
		this.axis = axis;
		this.angle = angle;
	}
	
	/**
	 * Construct the cone with the specified apex such that the points <code>p1</code>, <code>p2</code> and 
	 * <code>p3</code> are all on the boundary of the cone. 
	 * @param apex
	 * @param p1
	 * @param p2
	 * @param p3
	 */
	public InfCone(Point apex, Point p1, Point p2, Point p3){
		Vector v1 = apex.vectorTo(p1); v1.normalizeThis();
		Vector v2 = apex.vectorTo(p2); v2.normalizeThis();
		Vector v3 = apex.vectorTo(p3); v3.normalizeThis();

		Circle c = new Circle( apex.add(v1), apex.add(v2), apex.add(v3) );
		axis = apex.vectorTo(c.getCenter()).normalizeThis();
		angle = 2*Math.atan(c.getRadius()/apex.distance(c.getCenter()));
	}

	public Point getApex(){ return apex; }
	public Vector getAxis(){ return axis; }
	public double getAngle(){ return angle; }
	
	public String toString(int dec){
		return String.format("Cone[apex:%s,axis:%s,angle:%."+dec+"f]", apex,axis,angle);
	}
	public String toString(){ return toString(2); }
	
	public static void main(String[] args){
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		Point p0 = new Point(Point.getRandomPoint(3, -1.0, 1.0).getCoords());
		Point p1 = new Point(Point.getRandomPoint(3, -1.0, 1.0).getCoords());
		Point p2 = new Point(Point.getRandomPoint(3, -1.0, 1.0).getCoords());
		Point p3 = new Point(Point.getRandomPoint(3, -1.0, 1.0).getCoords());
		InfCone n = new InfCone(p0, p1, p2, p3);
		System.out.println(p0);
		System.out.println(n);
		Cone c = new Cone( n.getApex().add(n.getAxis()),n.getApex(), Math.tan(n.getAngle()/2) );
//		System.out.println(c);
//		Cone c = new Cone(new Point(0,0,1),new Point(1,1,0), 0.1);
		scene.addShape(new Sphere(p0,0.05), java.awt.Color.BLACK);
		scene.addShape(new Sphere(p1,0.05));
		scene.addShape(new Sphere(p2,0.05));
		scene.addShape(new Sphere(p3,0.05));
		scene.addShape(c, new java.awt.Color(0,0,200,240), 100);
		scene.setAxisEnabled(true);
	}
}
