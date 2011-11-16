package ProGAL.geom3d.volumes;

import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Vector;
import ProGAL.math.Matrix3x3;


public class Cylinder implements Volume{
	protected double rad;
	protected LineSegment segment;
	
	public Cylinder(Point p1, Point p2, double r) {
		this(new LineSegment(p1, p2), r);
	}
	
	public Cylinder(LineSegment sgm, double r) {
		this.rad = r;
		this.segment = sgm;
	}

	
	public LineSegment getSegment() { return segment; }
	
	public double getLength() { return 2*segment.getLength(); }
	public double getRadius() { return rad; }
	public double getVolume() { return Math.PI * rad * rad * getLength(); }
	public double getSurfaceArea() { return 2*Math.PI * rad * (rad + getLength()); }
	
	public boolean inCylinder(Point p) { return segment.getDistance(p) < rad; }
	
	public String toString(){ return toString(2); }
	
	public String toString(int dec) {
		return String.format("Cylinder[%s,rad=%."+dec+"f]", segment.toString(dec), rad);
	}

	
	public void toConsole(int dec) { System.out.println(toString(dec)) ; }
		
	public Point getCenter() {
		return segment.getMidPoint();
	}

	public static Cylinder createBoundingCylinder_CovarianceFit(PointList points) {
		if(points.size()<=0)	throw new Error("Cannot create cylinder enclosing 0 points");
		if(points.size()==1)	return new Cylinder(points.get(0).clone(), points.get(0).clone(), 0);
		if(points.size()==2)	return new Cylinder(points.get(0).clone(), points.get(1).clone(), 0);

		Matrix3x3 covMatr = points.getCovariance();
		Vector[] eigenVecs = covMatr.getEigenvectors();
		
		Vector dir = eigenVecs[0];
		if(eigenVecs[1].length()>dir.length()) dir = eigenVecs[1];
		if(eigenVecs[2].length()>dir.length()) dir = eigenVecs[2];

		InfCylinder iCyl = InfCylinder.createMinRadCylinderFromDirection(points, dir);
		Cylinder ret = iCyl.capWithDiscs(points);
		return ret;
	}

	public boolean overlaps(Volume vol) {
		throw new Error("Not implemented");
	}
	
	public Cylinder clone(){ 
		return new Cylinder(segment.clone(), rad);
	}

}
