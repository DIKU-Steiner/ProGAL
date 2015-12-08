package ProGAL.geom3d.volumes;


import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Line;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.PointWeighted;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.kineticDelaunay.Vertex;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.math.Constants;

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

	public Sphere(double x, double y, double z, double radius) {
		center = new Point(x, y, z);
		this.radius = radius;
	}

	
	public Sphere(Point[] ps) {
		this(0,0,0,0);
		radius = computeSphere_fast(ps[0],ps[1],ps[2],ps[3], center);
	}
	
	public Sphere(Point p0, Point p1, Point p2, Point p3) {
		this(0,0,0,0);
		radius = computeSphere_fast(p0,p1,p2,p3, center);
	}
	
	public Sphere(CTetrahedron tetr) {
		this(0,0,0,0);
		Point[] ps = tetr.getCorners();
		radius = computeSphere_fast(ps[0],ps[1],ps[2],ps[3], center);
	}
	
	/** Constructs a sphere with the weighted point as center and a radius with 
	 * the square root of the points weight. */
	public Sphere(PointWeighted p) {
		this.center = p;
		this.radius = Math.sqrt(p.getWeight());
	}

	/** creates a sphere with the specified circle as equator */
	public Sphere(Circle c) {
		center = c.getCenter();
		radius = c.getRadius();
	}
	
	/**
	 * Calculates the radius and center of the sphere touching the four specified points. The radius is 
	 * returned and if <code>center!=null</code> the coordinates of <code>center</code> are overwritten 
	 * with the coordinates of the circumsphere. Uses 46 "heavy" operations (multiplication/division).
	 */
	public static double computeSphere_fast(Point p0, Point p1,  Point p2, Point p3, Point center ) {
		double x1 = p1.x()-p0.x(), y1 = p1.y()-p0.y(), z1 = p1.z()-p0.z();
		double x2 = p2.x()-p0.x(), y2 = p2.y()-p0.y(), z2 = p2.z()-p0.z();
		double x3 = p3.x()-p0.x(), y3 = p3.y()-p0.y(), z3 = p3.z()-p0.z();
		
		double xx1 = x1*x1 + y1*y1 + z1*z1; //3HOp
		double xx2 = x2*x2 + y2*y2 + z2*z2, xx3 = x3*x3 + y3*y3 + z3*z3;  //6HOp
		
		double x1y2 = x1*y2, x1y3 = x1*y3, x1z2 = x1*z2, x1z3 = x1*z3;   //4HOp
		double x2y1 = x2*y1, x2y3 = x2*y3, x2z1 = x2*z1, x2z3 = x2*z3;   //4HOp
		double x3y2 = x3*y2, x3y1 = x3*y1, x3z2 = x3*z2, x3z1 = x3*z1;   //4HOp

		double y1z2 = y1*z2, y1z3 = y1*z3; //2HOp
		double y2z1 = y2*z1, y2z3 = y2*z3; //2HOp
		double y3z1 = y3*z1, y3z2 = y3*z2; //2HOp
		
		double m11 = -((x1y2-x2y1)*z3 + (x3y1-x1y3)*z2 + (x2y3-x3y2)*z1);  //3HOp
			
		if (m11 != 0.0) {
			
			double m12 = -(xx1*(y2z3-y3z2) + xx3*(y1z2-y2z1) + xx2*(y3z1-y1z3)); //3HOp
			double m13 = -(xx1*(x2z3-x3z2) + xx3*(x1z2-x2z1) + xx2*(x3z1-x1z3)); //3HOp
			double m14 = -(xx1*(x2y3-x3y2) + xx3*(x1y2-x2y1) + xx2*(x3y1-x1y3)); //3HOp
	
			double m11x2 = 0.5/m11; //1HOp
			double x = m12*m11x2; //1HOp
			double y =-m13*m11x2; //1HOp
			double z = m14*m11x2; //1HOp
			
			if(center!=null){
				center.setX(x+p0.x());
				center.setY(y+p0.y());
				center.setZ(z+p0.z());
			}
		    return Math.sqrt(x*x + y*y + z*z); //3HOp
		}
		else throw new RuntimeException("Points are coplanar");
	}
	
	/**
	 * 96 HOps
	 */
	private void computeSphere(Point p0, Point p1,  Point p2, Point p3 ) {
		double x0 = p0.x(); double y0 = p0.y(); double z0 = p0.z();
		double x1 = p1.x(); double y1 = p1.y(); double z1 = p1.z();
		double x2 = p2.x(); double y2 = p2.y(); double z2 = p2.z();
		double x3 = p3.x(); double y3 = p3.y(); double z3 = p3.z();
		
		double xx0 = x0*x0 + y0*y0 + z0*z0, xx1 = x1*x1 + y1*y1 + z1*z1; // 6HOps
		double xx2 = x2*x2 + y2*y2 + z2*z2, xx3 = x3*x3 + y3*y3 + z3*z3; // 6HOps
		
		double x1y2 = x1*y2, x1y3 = x1*y3, x1z2 = x1*z2, x1z3 = x1*z3; // 4HOps
		double x2y1 = x2*y1, x2y3 = x2*y3, x2z1 = x2*z1, x2z3 = x2*z3; // 4HOps
		double x3y2 = x3*y2, x3y1 = x3*y1, x3z2 = x3*z2, x3z1 = x3*z1; // 4HOps

		double y1z2 = y1*z2, y1z3 = y1*z3; // 2HOps
		double y2z1 = y2*z1, y2z3 = y2*z3; // 2HOps
		double y3z1 = y3*z1, y3z2 = y3*z2; // 2HOps
		
		
		double m11 =  x0*(y1z2 + y3z1 + y2z3 - y1z3 - y2z1 - y3z2)
		             -y0*(x1z2 + x3z1 + x2z3 - x1z3 - x2z1 - x3z2)
		             +z0*(x1y2 + x3y1 + x2y3 - x1y3 - x2y1 - x3y2)
		             -((x1y2-x2y1)*z3 + (x3y1-x1y3)*z2 + (x2y3-x3y2)*z1); // 6HOps
			
		if (m11 != 0.0) {
			
			double m12 =  xx0*(y1z2 + y3z1 + y2z3 - y1z3 - y2z1 - y3z2)
            -y0*(xx1*(z2-z3)     + xx3*(z1-z2)     + xx2*(z3-z1))
            +z0*(xx1*(y2-y3)     + xx3*(y1-y2)     + xx2*(y3-y1))
               -(xx1*(y2z3-y3z2) + xx3*(y1z2-y2z1) + xx2*(y3z1-y1z3)); // 12HOps
		
			double m13 =  xx0*(x1z2 + x3z1 + x2z3 - x1z3 - x2z1 - x3z2)
			-x0*(xx1*(z2-z3)     + xx3*(z1-z2)     + xx2*(z3-z1))
            +z0*(xx1*(x2-x3)     + xx3*(x1-x2)     + xx2*(x3-x1))
               -(xx1*(x2z3-x3z2) + xx3*(x1z2-x2z1) + xx2*(x3z1-x1z3)); // 12HOps

			double m14 =  xx0*(x1y2 + x3y1 + x2y3 - x1y3 - x2y1 - x3y2)
            -x0*(xx1*(y2-y3)     + xx3*(y1-y2)     + xx2*(y3-y1))
            +y0*(xx1*(x2-x3)     + xx3*(x1-x2)     + xx2*(x3-x1))
               -(xx1*(x2y3-x3y2) + xx3*(x1y2-x2y1) + xx2*(x3y1-x1y3)); // 12HOps

			double m15 =  xx0*(z3*(x1y2-x2y1) + z2*(x3y1-x1y3) + z1*(x2y3-x3y2))
            -x0*(xx1*(y2z3-y3z2) + xx3*(y1z2-y2z1) + xx2*(y3z1-y1z3))
            +y0*(xx1*(x2z3-x3z2) + xx3*(x1z2-x2z1) + xx2*(x3z1-x1z3))
            -z0*(xx1*(x2y3-x3y2) + xx3*(x1y2-x2y1) + xx2*(x3y1-x1y3)); // 16HOps
	
			double m11x2 = 0.5/m11; // 1HOps
			double x = m12*m11x2; // 1HOps
			double y =-m13*m11x2; // 1HOps
			double z = m14*m11x2; // 1HOps
		    center = new Point(x, y, z);
		    radius = Math.sqrt(x*x + y*y + z*z - m15/m11); // 4HOps
		}
		else System.out.println("Points are coplanar");
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
		return center.distanceSquared(p) < getRadiusSquared(); 
	}
	
	/** Returns true iff the squared distance from the point to the sphere center is less than the squared 
	 * radius minus <code>eps</code> 
	 */
	public boolean isInside(Point p, double eps) { return center.distanceSquared(p) < radius*radius - eps; }
	
	public boolean isEmpty(Point[] points, double eps) {
		for (int i = 0; i < points.length; i++) if (isInside(points[i], eps)) return false;
		return true;
	}
	/** returns TRUE if the interior of the sphere (for a given eps reduction of the radius) is empty */
	public boolean isEmpty(List<Vertex> points, double eps) {
		for (Point p : points) if (isInside(p, eps)) return false;
		return true;
	}
	
	public void contains(List<Vertex> points, double eps) {
		for (Vertex p : points) 
			if (isInside(p, eps)) 
				System.out.print(p.getId() + " " + (radius*radius) + " " + center.distanceSquared(p) + ", ");
		System.out.println();
	}

	public void setCenter(Point center) { this.center = center; }
	
	public void setCenter(Point p0, Point p1, Point p2, Point p3) { computeSphere(p0, p1, p2, p3); }

	public void setRadius(double radius) { this.radius = radius; }
	
	/** Returns true if this sphere is intersected or touched by another sphere. */
	public boolean isIntersected (Sphere sphere) {	return overlaps(sphere);	}

	/** Gets the secant on the line. TODO: Rename.*/
	public LineSegment getIntersection(Line line) {
		Point p1 = line.getP();
		Point p2 = line.getPoint(1.0);
		double dx = p2.x() - p1.x();
		double dy = p2.y() - p1.y();
		double dz = p2.z() - p1.z();
		double ex = p1.x() - center.x();
		double ey = p1.y() - center.y();
		double ez = p1.z() - center.z();
		double a = dx*dx + dy*dy + dz*dz;
		double b = 2*(dx*ex + dy*ey + dz*ez);
		double c = center.x()*center.x() + center.y()*center.y() + center.z()*center.z() + 
		p1.x()*p1.x() + p1.y()*p1.y() + p1.z()*p1.z() - 
		2*(center.x()*p1.x() + center.y()*p1.y() + center.z()*p1.z()) - radius*radius;
		double delta = b*b - 4*a*c; 
		if (delta < 0) return null;
		double u1, u2;
		if (delta == 0) u1 = u2 = - b/(2*a);
		else {
			double sqr = Math.sqrt(delta);
			u1 = (-b + sqr)/(2*a);
			u2 = (-b - sqr)/(2*a);
		}
		return new LineSegment(new Point(p1.x() + u1*dx, p1.y() + u1*dy, p1.z() + u1*dz),
				new Point(p1.x() + u2*dx, p1.y() + u2*dy, p1.z() + u2*dz));
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

	public Point[] getIntersections(Circle c) {
		Plane plane = new Plane(c.getCenter(), c.getNormalVector());
		Circle c2 = plane.getIntersection(this);
		if (c2 != null) return c.getIntersection(c2); else return null;
	}

	public Double getIntersectionAngle(Circle c, Point p, Vector dir) {
		Plane plane = new Plane(c.getCenter(), c.getNormalVector());                           
		Circle c2 = plane.getIntersection(this);
		if (c2 != null) return c.getFirstIntersection(c2, p, dir); else return null;
	}
	
	
	
	
	/** Returns true if none of the given points is in the sphere. */
	public boolean containsNone(List<Point> points) {
		double rr = radius*radius-0.000000001;
		for(Point p: points)
			if(p.distanceSquared(center)<rr) return false;
		return true;
	}

	/** Returns number of points inside the sphere */
	public int containsNumber(List<Vertex> points) {
		int count = 0;
		double rr = radius*radius-0.000000001;
		for(Point p: points)
			if(p.distanceSquared(center)<rr) count++;
		return count;
	}
	
	public boolean containsNoneButAtMostOne(Vertex v, List<Vertex> points) {
		double rr = radius*radius-0.000000001;
		for(Vertex p: points)
			if ((p.distanceSquared(center) < rr) && (p != v)) {
				System.out.println("Contains vertex " + p.getId());
				return false;
			}
		return true;
	}
	
	/** Returns true if the given point is in the sphere. */
	public boolean contains(Point p) {
		double rr = radius*radius-0.000000001;
		return p.distanceSquared(center)<rr;
	}


	/** Gets the squared distance of a point from a sphere surface 
	 * (negative if the point is inside the sphere). */
	public double powerDistance(Point p) {
		return center.distanceSquared(p) - radius*radius; 
	}

	/** Return a string representation of this sphere. */
	public String toString() {
		return toString(2);
	}

	/** Return a string representation of this sphere with <code>dec</code> decimals precision */
	public String toString(int dec) {
		return String.format("Sphere3d[%s,%"+dec+"f]",center.toString(dec), radius);
	}

	/** Writes this sphere to <code>System.out</code>. */
	public void toConsole() { System.out.println(toString()); }

	/** Writes this sphere to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	public void toScene(J3DScene scene, Color clr) {
		scene.addShape(this, clr);
	}

	public void toScene(J3DScene scene, Color clr, int res) {
		scene.addShape(this, clr, res);
	}

	public void toSceneSpiral(J3DScene scene, Color clr, int s, int step, double width) {
		double alpha;
		double delta;
		double t;
		double iPI;
		double tStep = 2.0/step;
		double r2 = radius * radius;
		for (int i = 1; i <= s; i++) {
			iPI = i*Math.PI;
			t = -1;
			for (int j = 0; j < step; j++) {
				delta = Math.sqrt(r2 - r2 * t * t);
				alpha = t * iPI;
				new Point(delta*Math.cos(alpha) + center.x(), delta*Math.sin(alpha) + center.y(), radius*t + center.z()).toScene(scene, width, clr);

				t += tStep;
			}
		}
		
	}
	
	/** Returns true if the sphere overlaps with <code>vol</code>. TODO: Implement for all volumes. */
	public boolean overlaps(Volume vol) {
		if(vol instanceof Sphere) 
			return ((Sphere)vol).center.distance(this.center)<=((Sphere)vol).radius+radius;
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
		return new Sphere( Point.getMidpoint(p1,p2), p1.distance(p2)/2 );
	}

	/** Get the smallest sphere through three points. */
	public static Sphere getMinSphere(Point p0, Point p1, Point p2) {
		Point center = new Point((p0.x()+p1.x()+p2.x())/3, (p0.y()+p1.y()+p2.y())/3, (p0.z()+p1.z()+p2.z())/3);
		double radius = p0.distance(center);
		return new Sphere(center, radius);
	}
	
	/** Constructs the sphere through four points. An error is thrown 
	 * if the points are coplanar. */ 
	public static Sphere getMinSphere(Point p0, Point p1, Point p2, Point p3) {
		Sphere ret = new Sphere(0,0,0,0);
		ret.radius = computeSphere_fast(p0, p1, p2, p3, ret.center);
		return ret;
//		double x0 = p0.x(); double y0 = p0.y(); double z0 = p0.z();
//		double x1 = p1.x(); double y1 = p1.y(); double z1 = p1.z();
//		double x2 = p2.x(); double y2 = p2.y(); double z2 = p2.z();
//		double x3 = p3.x(); double y3 = p3.y(); double z3 = p3.z();
//
//		double xx0 = x0*x0 + y0*y0 + z0*z0, xx1 = x1*x1 + y1*y1 + z1*z1;
//		double xx2 = x2*x2 + y2*y2 + z2*z2, xx3 = x3*x3 + y3*y3 + z3*z3;
//
//		double x1y2 = x1*y2, x1y3 = x1*y3, x1z2 = x1*z2, x1z3 = x1*z3;
//		double x2y1 = x2*y1, x2y3 = x2*y3, x2z1 = x2*z1, x2z3 = x2*z3; 
//		double x3y2 = x3*y2, x3y1 = x3*y1, x3z2 = x3*z2, x3z1 = x3*z1;
//
//		double y1z2 = y1*z2, y1z3 = y1*z3;
//		double y2z1 = y2*z1, y2z3 = y2*z3;
//		double y3z1 = y3*z1, y3z2 = y3*z2;
//
//
//		double m11 =  x0*(y1z2 + y3z1 + y2z3 - y1z3 - y2z1 - y3z2)
//		-y0*(x1z2 + x3z1 + x2z3 - x1z3 - x2z1 - x3z2)
//		+z0*(x1y2 + x3y1 + x2y3 - x1y3 - x2y1 - x3y2)
//		-((x1y2-x2y1)*z3 + (x3y1-x1y3)*z2 + (x2y3-x3y2)*z1);
//
//		if (m11 != 0.0) {
//
//			double m12 =  xx0*(y1z2 + y3z1 + y2z3 - y1z3 - y2z1 - y3z2)
//			-y0*(xx1*(z2-z3)     + xx3*(z1-z2)     + xx2*(z3-z1))
//			+z0*(xx1*(y2-y3)     + xx3*(y1-y2)     + xx2*(y3-y1))
//			-(xx1*(y2z3-y3z2) + xx3*(y1z2-y2z1) + xx2*(y3z1-y1z3));
//
//			double m13 =  xx0*(x1z2 + x3z1 + x2z3 - x1z3 - x2z1 - x3z2)
//			-x0*(xx1*(z2-z3)     + xx3*(z1-z2)     + xx2*(z3-z1))
//			+z0*(xx1*(x2-x3)     + xx3*(x1-x2)     + xx2*(x3-x1))
//			-(xx1*(x2z3-x3z2) + xx3*(x1z2-x2z1) + xx2*(x3z1-x1z3));
//
//			double m14 =  xx0*(x1y2 + x3y1 + x2y3 - x1y3 - x2y1 - x3y2)
//			-x0*(xx1*(y2-y3)     + xx3*(y1-y2)     + xx2*(y3-y1))
//			+y0*(xx1*(x2-x3)     + xx3*(x1-x2)     + xx2*(x3-x1))
//			-(xx1*(x2y3-x3y2) + xx3*(x1y2-x2y1) + xx2*(x3y1-x1y3));
//
//			double m15 =  xx0*(z3*(x1y2-x2y1) + z2*(x3y1-x1y3) + z1*(x2y3-x3y2))
//			-x0*(xx1*(y2z3-y3z2) + xx3*(y1z2-y2z1) + xx2*(y3z1-y1z3))
//			+y0*(xx1*(x2z3-x3z2) + xx3*(x1z2-x2z1) + xx2*(x3z1-x1z3))
//			-z0*(xx1*(x2y3-x3y2) + xx3*(x1y2-x2y1) + xx2*(x3y1-x1y3));
//
//
//			double x =  0.5*m12/m11;
//			double y = -0.5*m13/m11;
//			double z =  0.5*m14/m11;
//			return new Sphere(new Point(x, y, z), Math.sqrt(x*x + y*y + z*z - m15/m11));
//		}
//		throw new Error("Points are coplanar");
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
	
	
	public static Circle getIntersection(Sphere s1, Sphere s2){
		double r1 = s1.radius;
		double r2 = s2.radius;
		double d = s1.center.distance(s2.center);
		if(d>r1+r2) return null;
		
		double h = ProGAL.geom2d.Triangle.calculateHeight(r1,r2,d);
		double d1 = Math.sqrt(r1*r1 - h*h);
		double d2 = Math.sqrt(r2*r2 - h*h);
		if(d2>d) d1*=-1;
		Vector normal = s1.center.vectorTo(s2.center).normalizeThis();
		Point center = s1.center.add(normal.multiply(d1));
		return new Circle(center, h, normal);
	}
	
	/**
	 * Find the two, one or zero points that is at the intersection of the three sphere shells.
	 * If the three spheres intersect in more than two points or one sphere contains another, 
	 * the result of this method not specified.  
	 * @hops 57-70
	 */
	public static Point[] getIntersections(Sphere s1, Sphere s2, Sphere s3){
		
		Circle i12 = getIntersection(s1,s2);
		if(i12==null) return new Point[]{};
		
		return s3.getIntersections(i12);
		
		//Inspired by http://mathforum.org/library/drmath/view/63138.html
		//Theres a bug. Above works
//		double x1 = s1.center.x();
//		double y1 = s1.center.y();
//		double z1 = s1.center.z();
//		double c1Sq = s1.center.dot(s1.center);								//3HOp
//		double c2Sq = s2.center.dot(s2.center);								//3HOp
//		double r1 = s1.radius,							r1Sq = r1*r1;		//2HOp
//		Vector v12 = s2.center.vectorTo(s1.center);
//		Vector v23 = s3.center.vectorTo(s2.center);
//		double c12 = c1Sq-c2Sq;	
//		double c23 = c2Sq-s3.center.dot(s3.center);							//3HOp
//		double r2Sq = s2.radius*s2.radius;									//1HOp (12)
//		double r12 = r2Sq-r1Sq + c12;
//		double r23 = s3.radius*s3.radius - r2Sq + c23;						//3HOp
//		double Dyz = v12.y()*v23.z()-v23.y()*v12.z(),	DyzSq = Dyz*Dyz;	//3HOp
//		double Dry = r12*v23.y()-r23*v12.y(), 			DrySq = Dry*Dry;	//3HOp
//		double Dzx = v12.z()*v23.x()-v23.z()*v12.x(), 	DzxSq = Dzx*Dzx;	//3HOp
//		double Dxr = v12.x()*r23-v23.x()*r12,			DxrSq = Dxr*Dxr;	//3HOp
//		double Dxy = v12.x()*v23.y()-v23.x()*v12.y(), 	DxySq = Dxy*Dxy;	//3HOp (30)
//				
//		double k2 = (DyzSq+DzxSq)/DxySq + 1;								//1HOp
//		double k1 = (Dyz*Dry+Dzx*Dxr)/DxySq - 4*(x1*Dyz+y1*Dzx)/Dxy - z1;	//7HOp
//		double k0 = (DrySq+DxrSq)/(4*DxySq) + c1Sq-r1Sq - 2*(x1*Dry+y1*Dxr)/Dxy;	//6HOp
//				
//		double d = k1*k1-4*k2*k0;											//3HOp
//		if(d<-Constants.EPSILON) return new Point[]{};
//		if(d<Constants.EPSILON) { //One intersection-point
//			double zi = -k1/(2*k2);											//2HOp (49)
//			return new Point[]{
//					new Point( (2*zi*Dyz+Dry)/(2*Dxy) , (2*zi*Dzx+Dxr)/(2*Dxy) , zi ) 	//8HOp (57)
//			};
//		}
//		double dRt = Math.sqrt(d);											//1HOp (50)
//		double zi0 = (-k1-dRt)/(2*k2);										//2HOp
//		double zi1 = (-k1+dRt)/(2*k2);										//2HOp 
//				
//		return new Point[]{
//				new Point( (2*zi0*Dyz+Dry)/(2*Dxy) , (2*zi0*Dzx+Dxr)/(2*Dxy) , zi0 ),	//8HOp
//				new Point( (2*zi1*Dyz+Dry)/(2*Dxy) , (2*zi1*Dzx+Dxr)/(2*Dxy) , zi1 )	//8HOp (70)
//		};
	}
	
	
	//Daisy
	private double angle(Vector v1, Vector v2) {
		double ret = Math.atan2(v1.y(), v1.x()) - Math.atan2(v2.y(), v2.x());
		if (ret<0) ret = Constants.TAU+ret;
		return ret;
	}
	
	// Daisy : numerical method
	private Double findAngle(int start, int end, Vertex A, Vertex B, Vertex C, Vertex D, int dir, double alphaVal) {
		Point Cnew = C.clone();
		Point Dnew = D.clone();
		if (dir == 0) {
			Cnew.rotationCW(new Vector(0,0,1), Math.toRadians(start));
			Dnew.rotationCW(new Vector(0,0,1), Math.toRadians(start));
		} else {
			Cnew.rotationCCW(new Vector(0,0,1), Math.toRadians(start));
			Dnew.rotationCCW(new Vector(0,0,1), Math.toRadians(start));
		}
		double radius = new Sphere(A, B, Cnew, Dnew).getRadius();
		if (Math.abs(radius-alphaVal)<Math.pow(10, -9)) return start*(1.0);
		System.out.println("Radius-alphaVal = "+(radius-alphaVal));
		double sign = Math.signum(radius-alphaVal);
		double newStart = start;
		double newEnd = end;
		double angle;
		double newSign;
		
		for (int i = 0 ; i<99 ; i++) {
			angle = (newEnd+newStart)/2.0;
			System.out.println("Angle = "+angle);
			Cnew = C.clone();
			Dnew = D.clone();
			if (dir == 0) {
				Cnew.rotationCW(new Vector(0,0,1), Math.toRadians(angle));
				Dnew.rotationCW(new Vector(0,0,1), Math.toRadians(angle));
			} else {
				Cnew.rotationCCW(new Vector(0,0,1), Math.toRadians(angle));
				Dnew.rotationCCW(new Vector(0,0,1), Math.toRadians(angle));
			}
			radius = new Sphere(A, B, Cnew, Dnew).getRadius();
			System.out.println("Radius-alphaVal = "+(radius-alphaVal));
			if (Math.abs(radius-alphaVal)<Math.pow(10, -9)) {
				
				return angle;
			}
			newSign = Math.signum(radius-alphaVal);
			if (newSign!=sign) {
				newEnd = angle;
			} else newStart = angle;
		}
		return null;
	}
	
	public void plotRadius(Vertex A, Vertex B, Vertex C, Vertex D, J3DScene scene) {
		Sphere s = new Sphere(A, B, C, D);
		System.out.println("Initial radius2 = "+s.getRadius());
		double R;
		Sphere p;
		Point Cnew;
		Point Dnew;
		for (int j = 1 ; j<360 ; j++) {
/*			double R2 = getCircumRadiusSquared(A, B, C, D, Math.toRadians(j), radiusNoAngle);
			System.out.println("radius^2 = "+R2);
			R = Math.sqrt(R2);
			System.out.println("radius = "+R);
			p = new Sphere(new Point(j*1.0, R, 0.0), 0.01);
			p.toScene(scene, Color.BLACK);*/
			Cnew = C.clone();
			Dnew = D.clone();
			Cnew.rotationCW(new Vector(0,0,1), Math.toRadians(j));
			Dnew.rotationCW(new Vector(0,0,1), Math.toRadians(j));
			Sphere sphere = new Sphere(A, B, Cnew, Dnew);
			double sphereRadius = sphere.getRadius();
			System.out.println("plot : radius = "+sphereRadius+" at j = "+j);
			if (sphereRadius>Math.pow(10, 40)) continue;
			Color c = Color.BLACK;
			if (sphereRadius<1.0) {
				c = Color.ORANGE;
			}
			p = new Sphere(new Point(j*1.0, sphereRadius, 0.0), 1.0);
//			System.out.println("p = "+p.toString());
			p.toScene(scene, c);
		}
	}
		
	public static void main(String[] args){
		Circle c = new Circle(new Point(61.608,-6.951,5.080), 1.480, new Vector(0.548,0.501,-0.670));
		Sphere s = new Sphere(new Point(61.648,-9.430,4.846), 1.425);
		System.out.println(s.getIntersections(c)[0]);
		System.out.println(s.getIntersections(c)[1]);
		if(true) return;
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
//		Sphere s1 = new Sphere( new Point(46.29, 7.24,79.23),31.929843);
//		Sphere s2 = new Sphere( new Point(44.33,31.27,72.11),44.430195);
//		Sphere s3 = new Sphere( new Point(17.73,16.28,56.71),35.639535);
		Sphere s1 = new Sphere( new Point(4.5,0.7,8.0),3.0);
		Sphere s2 = new Sphere( new Point(4.0,3.0,7.0),4.0);
		Sphere s3 = new Sphere( new Point(1.5,1.6,5.5),3.5);
//		Sphere s1 = new Sphere( new Point(2,0,0), 1);
//		Sphere s2 = new Sphere( new Point(0,0,0), 1.2);
//		Sphere s3 = new Sphere( new Point(1,1,0), 1);
		scene.addShape(s1, new Color(200,0,0), 50);
		scene.addShape(s2, new Color(0,200,0), 50);
		scene.addShape(s3, new Color(0,0,200), 50);
		Point[] intersections = Sphere.getIntersections(s1, s2, s3);
		for(Point p: intersections){
			System.out.println(p);
			scene.addShape(new Sphere(p,0.1), Color.GRAY.darker());
		}
	}

	/**
	 * Estimate the volume of the union of a set of spheres. A grid is placed
	 * around the spheres and the volume of the cells are used to compute an
	 * upper and a lower bound. The average value of the upper and lower bound 
	 * is used. Typically the result is accurate to 1/100 (not to 1/1000) and 
	 * takes around 150ms to compute.  
	 */
	public static double unionVolume_Grid(Collection<Sphere> spheres){
		Point[] centers = new Point[spheres.size()];
		double[] radSqs = new double[spheres.size()];
		double[] rads = new double[spheres.size()];
		Point minPoint = new Point(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
		Point maxPoint = new Point(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
		int c=0;
		for(Sphere s: spheres){
			centers[c] = s.center;
			radSqs[c] = s.getRadiusSquared();
			rads[c] = s.radius;
			c++;
			for(int i=0;i<3;i++){
				if(s.center.get(i)-s.radius<minPoint.get(i)) minPoint.set(i, s.center.get(i)-s.radius);
				if(s.center.get(i)+s.radius>maxPoint.get(i)) maxPoint.set(i, s.center.get(i)+s.radius);
			}
		}

		//Each dimension is divided into <code>cells</code> cells.
		int cells = 80;
		double[] delta = {
				(maxPoint.x()-minPoint.x())/cells,
				(maxPoint.y()-minPoint.y())/cells,
				(maxPoint.z()-minPoint.z())/cells
		};
		//Determine which grid-vertices are inside at least one sphere
		boolean[][][] bits = new boolean[cells+1][cells+1][cells+1];
		int xC=0, yC=0, zC=0;
		for(double x=minPoint.x();x<=maxPoint.x();x+=delta[0]){
			yC=0;
			for(double y=minPoint.y();y<=maxPoint.y();y+=delta[1]){
				zC=0;
				for(double z=minPoint.z();z<=maxPoint.z();z+=delta[2]){
					//Determine if (x,y,z) is inside a sphere
					for(int i=0;i<centers.length;i++){
						double dX = Math.abs(x-centers[i].x());
						double dY = Math.abs(y-centers[i].y());
						double dZ = Math.abs(z-centers[i].z());
						if(dX>rads[i] || dY>rads[i] || dZ>rads[i]) continue;
						if(dX*dX+dY*dY+dZ*dZ<radSqs[i]) {bits[xC][yC][zC] = true; break; } 
					}
					zC++;
				}
				yC++;
			}
			xC++;
		}

		//Determine how many cells are completely inside a sphere and how many are on the border
		int borderCells = 0;
		int insideCells = 0;
		for(int x=0;x<cells;x++){
			for(int y=0;y<cells;y++){
				for(int z=0;z<cells;z++){
					if(		bits[x][y][z] || 
							bits[x][y][z+1] || 
							bits[x][y+1][z] || 
							bits[x][y+1][z+1] || 
							bits[x+1][y][z] ||
							bits[x+1][y][z+1] ||  
							bits[x+1][y+1][z] ||  
							bits[x+1][y+1][z+1]){
						if(		bits[x][y][z] && 
								bits[x][y][z+1] && 
								bits[x][y+1][z] && 
								bits[x][y+1][z+1] && 
								bits[x+1][y][z] &&  
								bits[x+1][y][z+1] &&  
								bits[x+1][y+1][z] &&  
								bits[x+1][y+1][z+1]){
							insideCells++;
						}else{
							borderCells++;
						}
					}
				}
			}
		}

		double cellVol = delta[0]*delta[1]*delta[2];
		double insideVol = cellVol*insideCells;
		double borderVol = cellVol*borderCells;
		return insideVol+borderVol/2;
	}

	/** TODO: Comment, move up and test */
	public PointList generateRandomPointsOnSphere(int n){
		PointList ret = PointList.generateRandomPointsOnSphere(n);
		
		for(Point p: ret){
			p.scaleThis(radius);
			p.addThis(center.toVector());
		}
		return ret;
	}
	
	/** TODO: Comment, move up and test */
	public PointList generatePointsOnSphere(int n){
		PointList ret = PointList.generatePointsOnSphere(n);
		
		for(Point p: ret){
			p.scaleThis(radius);
			p.addThis(center.toVector());
		}
		return ret;
	}
}

