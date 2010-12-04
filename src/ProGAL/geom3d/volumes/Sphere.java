package ProGAL.geom3d.volumes;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Line;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Segment;
import ProGAL.geom3d.Vector;

/** 
 * A sphere represented by a center-point and a radius. 
 */
public class Sphere implements Volume{
	protected Point center;
	protected double radius;
	
	/** Constructs a sphere with the specified center and the specified radius. */
	public Sphere(Point center, double radius) {
		this.center = center;
		this.radius = radius;
	}

	/** Get the center */
	public Point getCenter() { return center; }
	/** Get the radius */
	public double getRadius() { return radius; }
	/** Get the squared radius */
	public double getRadiusSquared() { return radius*radius; }
	/** Get the surface area */
	public double getSurfaceArea() { return 4*Math.PI*radius*radius; }
	/** Get the volume */
	public double getVolume() { return getSurfaceArea()*radius/3; }
	/** Returns true if the point is inside this sphere */
	public boolean isInside(Point p) { 
		return center.getDistanceSquared(p) < getRadiusSquared(); 
	}
	
	/** Returns true if this sphere is intersected or touched by another sphere. */
	public boolean isIntersected (Sphere sphere) {	return overlaps(sphere);	}
	
	/** Gets the secant on the line. TODO: Rename.*/
	public Segment getIntersection(Line line) {
		Point p1 = line.getP();
		Point p2 = line.getPoint(1.0);
		double dx = p2.getX() - p1.getX();
		double dy = p2.getY() - p1.getY();
		double dz = p2.getZ() - p1.getZ();
		double ex = p1.getX() - center.getX();
		double ey = p1.getY() - center.getY();
		double ez = p1.getZ() - center.getZ();
		double a = dx*dx + dy*dy + dz*dz;
		double b = 2*(dx*ex + dy*ey + dz*ez);
		double c = center.getX()*center.getX() + center.getY()*center.getY() + center.getZ()*center.getZ() + 
		           p1.getX()*p1.getX() + p1.getY()*p1.getY() + p1.getZ()*p1.getZ() - 
		           2*(center.getX()*p1.getX() + center.getY()*p1.getY() + center.getZ()*p1.getZ()) - radius*radius;
		double delta = b*b - 4*a*c; 
		if (delta < 0) return null;
		double u1, u2;
		if (delta == 0) u1 = u2 = - b/(2*a);
		else {
			double sqr = Math.sqrt(delta);
			u1 = (-b + sqr)/(2*a);
			u2 = (-b - sqr)/(2*a);
		}
		return new Segment(new Point(p1.getX() + u1*dx, p1.getY() + u1*dy, p1.getZ() + u1*dz),
							 new Point(p1.getX() + u2*dx, p1.getY() + u2*dy, p1.getZ() + u2*dz));
	}
	

	
	/** 
	 * Returns the two line-parameters that indicate where <code>line</code> intersects 
	 * this sphere. TODO: Coordinate line-intersection methods (see above). 
	 */
	public double[] intersectionParameters(Line line) {
		Vector l = line.getDir();//.norm();
		Vector c = line.getP().vectorTo(center);
		double lc = l.dot(c);
		double cc = c.dot(c);
		double rr = radius*radius;
		double tmp = lc*lc-cc+rr;
		if(tmp<0) return new double[0];
		else if(tmp==0) return new double[]{lc};
		else {
			double d1 = lc-Math.sqrt(tmp);
			double d2 = lc+Math.sqrt(tmp);
			return new double[]{d1, d2};
		}
	}

	
	
	/** Returns true if none of the given points is in the sphere. */
	public boolean containsNone(PointList points) {
		double rr = radius*radius-0.000000001;
		for(Point p: points)
			if(p.getDistanceSquared(center)<rr) return false;
		return true;
	}
	
	
	/** Gets the squared distance of a point from a sphere surface 
	 * (negative if the point is inside the sphere). */
	public double powerDistance(Point p) {
		return center.getDistanceSquared(p) - radius*radius; 
	}
	
	/** Return a string representation of this sphere. */
	public String toString() {
		return String.format("Sphere3d[%s,%d]",center.toString(), radius);
	}
	
	/** Return a string representation of this sphere with <code>dec</code> decimals precision */
	public String toString(int dec) {
		return String.format("Sphere3d[%s,%"+dec+"d]",center.toString(dec), radius);
	}
	
	/** Writes this sphere to <code>System.out</code>. */
	public void toConsole() { System.out.println(toString()); }
	
	/** Writes this sphere to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	/** Returns true if the sphere overlaps with <code>vol</code>. TODO: Implement for all volumes. */
	public boolean overlaps(Volume vol) {
		if(vol instanceof Sphere) 
			return ((Sphere)vol).center.getDistance(this.center)<=((Sphere)vol).radius+radius;
		throw new IllegalArgumentException();
	}

	/** Returns a deep clone of this sphere. */
	public Sphere clone(){
		return new Sphere(center.clone(), radius);
	}

	/** Get the sphere with the specified circle as equator */
	public static Sphere getMinSphere(Circle c) {
		return new Sphere(c.getCenter(), c.getRadius());
	}
	
	/** Get the smallest sphere through two given points. */
	public static Sphere getMinSphere(Point p1, Point p2) {
		return new Sphere( Point.getMidpoint(p1,p2), p1.getDistance(p2)/2 );
	}
	
	/** Get the smallest sphere through three points. */
	public static Sphere getMinSphere(Point p0, Point p1, Point p2) {
		Point center = new Point((p0.getX()+p1.getX()+p2.getX())/3, (p0.getY()+p1.getY()+p2.getY())/3, (p0.getZ()+p1.getZ()+p2.getZ())/3);
		double radius = p0.getDistance(center);
		return new Sphere(center, radius);
	}

	/** Constructs the smallest sphere through four points. An error is thrown 
	 * if the points are coplanar. */ 
	public static Sphere getMinSphere(Point p0, Point p1, Point p2, Point p3) {
		double x0 = p0.getX(); double y0 = p0.getY(); double z0 = p0.getZ();
		double x1 = p1.getX(); double y1 = p1.getY(); double z1 = p1.getZ();
		double x2 = p2.getX(); double y2 = p2.getY(); double z2 = p2.getZ();
		double x3 = p3.getX(); double y3 = p3.getY(); double z3 = p3.getZ();
		
		double xx0 = x0*x0 + y0*y0 + z0*z0, xx1 = x1*x1 + y1*y1 + z1*z1;
		double xx2 = x2*x2 + y2*y2 + z2*z2, xx3 = x3*x3 + y3*y3 + z3*z3;
		
		double x1y2 = x1*y2, x1y3 = x1*y3, x1z2 = x1*z2, x1z3 = x1*z3;
		double x2y1 = x2*y1, x2y3 = x2*y3, x2z1 = x2*z1, x2z3 = x2*z3; 
		double x3y2 = x3*y2, x3y1 = x3*y1, x3z2 = x3*z2, x3z1 = x3*z1;

		double y1z2 = y1*z2, y1z3 = y1*z3;
		double y2z1 = y2*z1, y2z3 = y2*z3;
		double y3z1 = y3*z1, y3z2 = y3*z2;
		
		
		double m11 =  x0*(y1z2 + y3z1 + y2z3 - y1z3 - y2z1 - y3z2)
		             -y0*(x1z2 + x3z1 + x2z3 - x1z3 - x2z1 - x3z2)
		             +z0*(x1y2 + x3y1 + x2y3 - x1y3 - x2y1 - x3y2)
		             -((x1y2-x2y1)*z3 + (x3y1-x1y3)*z2 + (x2y3-x3y2)*z1);
			
		if (m11 != 0.0) {
			
			double m12 =  xx0*(y1z2 + y3z1 + y2z3 - y1z3 - y2z1 - y3z2)
            -y0*(xx1*(z2-z3)     + xx3*(z1-z2)     + xx2*(z3-z1))
            +z0*(xx1*(y2-y3)     + xx3*(y1-y2)     + xx2*(y3-y1))
               -(xx1*(y2z3-y3z2) + xx3*(y1z2-y2z1) + xx2*(y3z1-y1z3));
		
			double m13 =  xx0*(x1z2 + x3z1 + x2z3 - x1z3 - x2z1 - x3z2)
			-x0*(xx1*(z2-z3)     + xx3*(z1-z2)     + xx2*(z3-z1))
            +z0*(xx1*(x2-x3)     + xx3*(x1-x2)     + xx2*(x3-x1))
               -(xx1*(x2z3-x3z2) + xx3*(x1z2-x2z1) + xx2*(x3z1-x1z3));

			double m14 =  xx0*(x1y2 + x3y1 + x2y3 - x1y3 - x2y1 - x3y2)
            -x0*(xx1*(y2-y3)     + xx3*(y1-y2)     + xx2*(y3-y1))
            +y0*(xx1*(x2-x3)     + xx3*(x1-x2)     + xx2*(x3-x1))
               -(xx1*(x2y3-x3y2) + xx3*(x1y2-x2y1) + xx2*(x3y1-x1y3));

			double m15 =  xx0*(z3*(x1y2-x2y1) + z2*(x3y1-x1y3) + z1*(x2y3-x3y2))
            -x0*(xx1*(y2z3-y3z2) + xx3*(y1z2-y2z1) + xx2*(y3z1-y1z3))
            +y0*(xx1*(x2z3-x3z2) + xx3*(x1z2-x2z1) + xx2*(x3z1-x1z3))
            -z0*(xx1*(x2y3-x3y2) + xx3*(x1y2-x2y1) + xx2*(x3y1-x1y3));

	
		    double x =  0.5*m12/m11;
		    double y = -0.5*m13/m11;
		    double z =  0.5*m14/m11;
		    return new Sphere(new Point(x, y, z), Math.sqrt(x*x + y*y + z*z - m15/m11));
		}
		throw new Error("Points are coplanar");
	}
	

	
	/**
	 * Gets the smallest sphere containing a set of points. Uses a 
	 * randomized, O(n) expected time algorithm. 
	 */
	public static Sphere getMinSphere(PointList points) {
		return getMinSphere(points.getRandomPermutation(), points.size(), new PointList());
	}
	
	private static Sphere getMinSphere(PointList points, int n, PointList boundaryPoints) {
		Sphere sphere = null;
		int k = 0;
		switch (boundaryPoints.size()) {
			case 0: sphere = getMinSphere(points.get(0), points.get(1)); k = 2; break;
			case 1: sphere = getMinSphere(points.get(0), boundaryPoints.get(0)); k = 1; break;
			case 2: sphere = getMinSphere(boundaryPoints.get(0), boundaryPoints.get(1)); break;
			case 3: sphere = getMinSphere(boundaryPoints.get(0), boundaryPoints.get(1), boundaryPoints.get(2)); break;
		}
		
		for (int i = k; i < n + boundaryPoints.size(); i++) {
			Point p = (Point)points.get(i);
			if (!boundaryPoints.contains(p)) {
				if (!sphere.isInside(p)) {
					if (boundaryPoints.size() < 3) {
						boundaryPoints.add(p);
						sphere = getMinSphere(points, i-1, boundaryPoints);
						boundaryPoints.remove(p);
					}
					else sphere = getMinSphere(boundaryPoints.get(0), boundaryPoints.get(1), boundaryPoints.get(2), p);
				}
			}
		}
		return sphere;
	}
	
}

