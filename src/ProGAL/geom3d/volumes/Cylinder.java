package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Line;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants;
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
	
	public void setSegment(LineSegment sgm) { segment = sgm; }
	
	public boolean inCylinder(Point p) { return segment.getDistance(p) < rad; }
	
	public Double intersectionParameter(Line l){
		double t;
		Vector d = segment.getAToB();
		Vector m = segment.getA().vectorTo(l.getP());
		Vector n = l.getDir();
		double md = m.dot(d);
		double nd = n.dot(d);
		double dd = d.dot(d);
		// Test if segment fully outside either endcap of cylinder
		if (md < 0.0f && md + nd < 0.0f) return null; // Segment outside ’p’ side of cylinder 
		if (md > dd && md + nd > dd) return null; // Segment outside ’q’ side of cylinder 
		double nn = n.dot(n);
		double mn = m.dot(n);
		double a = dd*nn - nd*nd;
		double k = m.dot(m) - rad*rad;
		double c = dd * k - md * md;
		if (Math.abs(a) < Constants.EPSILON) {
			// Segment runs parallel to cylinder axis
			if (c > 0.0f) return null; // ’a’ and thus the segment lie outside cylinder 
			// Now known that segment intersects cylinder; figure out how it intersects 
			if (md < 0.0f) t = -mn / nn; // Intersect segment against ’p’ endcap
			else if (md > dd) t = (nd - mn) / nn; // Intersect segment against ’q’ endcap 
			else t = 0.0f; // ’a’ lies inside cylinder

			if(t>0) return t;
			else return null;
		}
		double b = dd * mn - nd * md;
		double discr = b * b - a * c;
		if (discr < 0.0f) return null;// No real roots; no intersection

		t = (-b - Math.sqrt(discr)) / a;
		double t0 = t = (-b - Math.sqrt(discr)) / a;
	    if (md + t * nd < 0.0f) {
	        // Intersection outside cylinder on ‘p’ side
	        if (nd <= 0.0f) return null; // Segment pointing away from endcap
	        t = -md / nd;
	        // Keep intersection if Dot(S(t) - p, S(t) - p) <= r^2
	        if (k + t * (2.0 * mn + t * nn) <= 0.0f) return t;
	    } else if (md + t * nd > dd) {
	        // Intersection outside cylinder on ‘q’ side
	        if (nd >= 0.0f) return 0.0; // Segment pointing away from endcap
	        t = (dd - md) / nd;
	        // Keep intersection if Dot(S(t) - q, S(t) - q) <= r^2
	        if (k + dd - 2.0f * md + t * (2.0f * (mn - nd) + t * nn) <= 0.0f) return t;
	    }
	    t = t0;
	    return t;
	}
	
	public static void main(String[] args){
		Line l = new Line(new Point(0,1.01,0), new Vector(1,0,0));
		Cylinder cyl = new Cylinder(new Point(3,1,0), new Point(3,-1,0), 1);
		System.out.println(cyl.intersectionParameter(l));
	}
	
	
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
