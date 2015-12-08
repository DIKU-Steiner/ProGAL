package ProGAL.geom3d.volumes;

import java.awt.Color;

import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;

/** 
 * A cone consists of all the points in an infinite cone that are less than a certain distance from 
 * the apex along the axis. The circular area furthest along the axis from the apex is called the base 
 * of the cone.
 * @author R.Fonseca
 */
public class Cone extends InfCone implements Volume{
//	private Point p1, p2;
//	private double rad;
	private double axisLength, baseRad;

	public Cone(Point apex, Vector axis, double angle, double axisLength){
		super(apex, axis, angle);
		this.axisLength = axisLength;
		this.baseRad = Math.tan(angle)*axisLength;
	}

	/**
	 * Construct a cone from the apex-point, the center of the base and the radius of the base.
	 * @param apex
	 * @param baseCenter
	 * @param baseRad
	 */
	public Cone(Point apex, Point baseCenter, double baseRad){
		super(apex,apex.vectorTo(baseCenter).normalizeThis(), Math.atan(baseRad/apex.distance(baseCenter)));
		this.axisLength = apex.distance(baseCenter);
		this.baseRad = baseRad;
	}
	
	public double getAxisLength(){
		return axisLength;
	}

	public Point getBaseCenter(){
		return apex.add(axis.multiply(axisLength));
	}
	
	public double getBaseRadius() {
		return baseRad;
	}
	
	@Override
	public double getVolume() {
		return Math.PI*baseRad*baseRad*axisLength/3.0;
	}

	@Override
	public Point getCenter() {
		return apex.add(axis.multiply(axisLength/2));
	}

	@Override
	public Cone clone(){
		return new Cone(apex, axis, angle, axisLength);
	}

	@Override
	public boolean overlaps(Volume vol) {
		throw new RuntimeException("Overlaps not implemented");
	}
	
	public static void main(String[] args){
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		Cone c1 = new Cone(new Point(0,0,0), new Vector(-1,0,0), Math.PI*2.0/3.0, 1);
		Cone c2 = new Cone(new Point(0,0,0), new Vector(0.5, Math.sqrt(3)/2,0), Math.PI*2.0/3.0, 1);
		Cone c3 = new Cone(new Point(0,0,0), new Vector(0.5,-Math.sqrt(3)/2,0), Math.PI*2.0/3.0, 1);
		scene.addShape(c1, new Color(200,0,0,100), 100);
//		scene.addShape(c2, Color.GREEN, 100);
//		scene.addShape(c3, Color.BLUE, 100);
		scene.addShape(new Sphere(new Point(-0.5,0,0), 0.01), Color.BLACK);
		scene.addShape(new LSS(new Point(-0.5,0,0), new Point(-0.5,1,0), 0.01));
		scene.addShape(new LSS(new Point(-0.5,0,0), new Point(-0.5,-0.5,Math.sqrt(3)/2), 0.01));
		scene.addShape(new LSS(new Point(-0.5,0,0), new Point(-0.5,-0.5,-Math.sqrt(3)/2), 0.01));
		Point p1 = new Point(-0.2,0,Math.tan(Math.PI/3.0)*0.2);
		double r = Math.tan(Math.PI/3.0)*0.6;
		Point p2 = new Point(-0.6,Math.cos(Math.PI/3)*r, Math.sin(Math.PI/3)*r);
		Point p3 = new Point(-0.6,Math.cos(-Math.PI/3)*r, Math.sin(-Math.PI/3)*r);
		scene.addShape(new Sphere(p1, 0.02));
		scene.addShape(new Sphere(p2, 0.02));
		scene.addShape(new Sphere(p3, 0.02));
		scene.addShape(new Triangle(p1, p2, p3), Color.GREEN);
		scene.addShape(new LSS(new Point(0,0,0),p1, 0.01));
		scene.addShape(new LSS(new Point(0,0,0),p2, 0.01));
		scene.addShape(new LSS(new Point(0,0,0),p3, 0.01));
//		scene.setAxisEnabled(true);
	}


}

