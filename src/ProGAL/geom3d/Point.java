package ProGAL.geom3d;

import java.awt.Color;

import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;

import ProGAL.math.Constants;

/** 
 *  A point in (x,y,z)-space represented with double precision. 
 */
public class Point extends ProGAL.geomNd.Point implements Simplex{
	private static final long serialVersionUID = -2120468832687547475L;

	/** Construct a point with the specified coordinates. */
	public Point(double x, double y, double z) { 
		super(new double[]{x,y,z});
	}
	/** Construct a point with the specified coordinates. */
	public Point(double[] coords) { 
		super(coords);
	}
	/** Construct a point that is a clone of p. */
	public Point(ProGAL.geomNd.Point p) { 
		super(p);
	}
	/** Construct a point at the coordinates of v. */
	public Point(Vector v) { 
		super(new double[]{v.x(),v.y(),v.z()});
	}

	
//	/** Get the i'th coordinate. */
//	public double getCoord(int i) {
//		switch(i){
//		case 0: return coords[0];
//		case 1: return coords[1];
//		case 2: return coords[2];
//		}
//		throw new Error("Trying to get invalid coordinate");
//	}
	
//	/** Get the i'th coordinate. */
//	public double get(int i) { return getCoord(i); }

	/** Get the first coordinate. */
	public final double x() { return coords[0]; }
	/** Get the second coordinate. */
	public final double y() { return coords[1]; }
	/** Get the third coordinate. */
	public final double z() { return coords[2]; }
	
	/** Set the first coordinate */
	public final void setX(double x) { this.coords[0] = x; }
	/** Set the second coordinate */
	public final void setY(double y) { this.coords[1] = y; }
	/** Set the third coordinate */
	public final void setZ(double z) { this.coords[2] = z; }
	
	/** 
	 * Return the 'dimension' of this object. Required by the interface Simplex.
	 * Beware not to confuse this method with getDimensions from geomNd.Point.  
	 */
	public final int getDimension() { return 0; }
	
	/** Get the vector that points from this point to p */
	public Vector vectorTo(Point p){
		return new Vector(p.coords[0]-coords[0], p.coords[1]-coords[1], p.coords[2]-coords[2]);
	}

	/** 
	 * Returns true if three points are on the same line. This implies that
	 * overlapping points are considered collinear.
	 */
	public static boolean collinear(Point p0, Point p1, Point p2) {
		Vector v1v0 = p1.vectorTo(p0);
		Vector v1v2 = p1.vectorTo(p2);
		return v1v0.cross(v1v2).getLengthSquared()<Constants.EPSILON;
	}

	/** Returns true if four specified points are in the same plane	 */
	public static boolean coplanar(Point p0, Point p1, Point p2, Point p3) {
		double ax = p0.coords[0]; double ay = p0.coords[1]; double az = p0.coords[2];
		double bx = p1.coords[0]; double by = p1.coords[1]; double bz = p1.coords[2];
		double cx = p2.coords[0]; double cy = p2.coords[1]; double cz = p2.coords[2];
		double dx = p3.coords[0]; double dy = p3.coords[1]; double dz = p3.coords[2];
		return   Math.abs(
				-az*by*cx + ay*bz*cx + az*bx*cy - ax*bz*cy
				-ay*bx*cz + ax*by*cz + az*by*dx - ay*bz*dx
				-az*cy*dx + bz*cy*dx + ay*cz*dx - by*cz*dx
				-az*bx*dy + ax*bz*dy + az*cx*dy - bz*cx*dy
				-ax*cz*dy + bx*cz*dy + ay*bx*dz - ax*by*dz
				-ay*cx*dz + by*cx*dz + ax*cy*dz - bx*cy*dz	) < Constants.EPSILON; 
	}

	/* returns a positive double if point s is above the plane through p, q, r */
	public static double orientation(Point p, Point q, Point r, Point s) {
		double M1 = q.x()*(r.y()-s.y()) + r.x()*(s.y()-q.y()) + s.x()*(q.y()-r.y());
		double M2 = p.x()*(r.y()-s.y()) + r.x()*(s.y()-p.y()) + s.x()*(p.y()-r.y());
		double M3 = p.x()*(q.y()-s.y()) + q.x()*(s.y()-p.y()) + s.x()*(p.y()-q.y());
		double M4 = p.x()*(q.y()-r.y()) + q.x()*(r.y()-p.y()) + r.x()*(p.y()-q.y());
		return -p.z()*M1 + q.z()*M2 - r.z()*M3 + s.z()*M4;
	}

	/* returns a positive double if points s and t are on the same side of above the plane through p, q, r */
	/* This strange solution is needed as the ordering of p, q, r in tetrahedra is not known. It would be much 
	 * easier if the ordering was counterclockwise w.r.t. the fourth corner t
	 */
	public static double orientation(Point p, Point q, Point r, Point s, Point t) {
		double pqy = p.y()-q.y();
		double pry = p.y()-r.y();
		double psy = p.y()-s.y();
		double pty = p.y()-t.y();
		double qry = q.y()-r.y();
		double qsy = q.y()-s.y();
		double qty = q.y()-t.y();
		double rsy = r.y()-s.y();
		double rty = r.y()-t.y();
		double Ms1 = q.x()*rsy - r.x()*qsy + s.x()*qry;
		double Ms2 = p.x()*rsy - r.x()*psy + s.x()*pry;
		double Ms3 = p.x()*qsy - q.x()*psy + s.x()*pqy;
		double Mt1 = q.x()*rty - r.x()*qty + t.x()*qry;
		double Mt2 = p.x()*rty - r.x()*pty + t.x()*pry;
		double Mt3 = p.x()*qty - q.x()*pty + t.x()*pqy;
		double M4 =  p.x()*qry - q.x()*pry + r.x()*pqy;
		double Ms = -p.z()*Ms1 + q.z()*Ms2 - r.z()*Ms3 + s.z()*M4;
		double Mt = -p.z()*Mt1 + q.z()*Mt2 - r.z()*Mt3 + t.z()*M4;
		if (Ms*Mt >= 0.0) return Math.abs(Ms);
		return -Math.abs(Ms);
	}

	
	/** Returns a positive double if point t is inside the sphere through points p, q, r, s. */
	public static double inSphere(Point p, Point q, Point r, Point s, Point t) {
		double pp = p.x()*p.x() + p.y()*p.y() + p.z()*p.z();
		double qq = q.x()*q.x() + q.y()*q.y() + q.z()*q.z();
		double rr = r.x()*r.x() + r.y()*r.y() + r.z()*r.z();
		double ss = s.x()*s.x() + s.y()*s.y() + s.z()*s.z();
		double tt = t.x()*t.x() + t.y()*t.y() + t.z()*t.z();
		
		double M12 = r.x()*(s.y()-t.y()) + s.x()*(t.y()-r.y()) + t.x()*(r.y()-s.y()); 
		double M13 = q.x()*(s.y()-t.y()) + s.x()*(t.y()-q.y()) + t.x()*(q.y()-s.y());
		double M14 = q.x()*(r.y()-t.y()) + r.x()*(t.y()-q.y()) + t.x()*(q.y()-r.y());
		double M15 = q.x()*(r.y()-s.y()) + r.x()*(s.y()-q.y()) + s.x()*(q.y()-r.y());
		double M23 = p.x()*(s.y()-t.y()) + s.x()*(t.y()-p.y()) + t.x()*(p.y()-s.y());
		double M24 = p.x()*(r.y()-t.y()) + r.x()*(t.y()-p.y()) + t.x()*(p.y()-r.y());
		double M25 = p.x()*(r.y()-s.y()) + r.x()*(s.y()-p.y()) + s.x()*(p.y()-r.y());
		double M34 = p.x()*(q.y()-t.y()) + q.x()*(t.y()-p.y()) + t.x()*(p.y()-q.y());
		double M35 = p.x()*(q.y()-s.y()) + q.x()*(s.y()-p.y()) + s.x()*(p.y()-q.y());
		double M45 = p.x()*(q.y()-r.y()) + q.x()*(r.y()-p.y()) + r.x()*(p.y()-q.y());
		                         
		double M1 = -q.z()*M12 + r.z()*M13 - s.z()*M14 + t.z()*M15;
		double M2 = -p.z()*M12 + r.z()*M23 - s.z()*M24 + t.z()*M25;
		double M3 = -p.z()*M13 + q.z()*M23 - s.z()*M34 + t.z()*M35;
		double M4 = -p.z()*M14 + q.z()*M24 - r.z()*M34 + t.z()*M45;
		double M5 = -p.z()*M15 + q.z()*M25 - r.z()*M35 + s.z()*M45;
		
		return pp*M1 - qq*M2 + rr*M3 - ss*M4 + tt*M5;
	}

	/** Translates this point by (x,y,z). */
	public void translateThis(double dx, double dy, double dz) {
		this.coords[0] += dx;
		this.coords[1] += dy;
		this.coords[2] += dz;
	}
	/** Translates this point by p. */
	public void translateThis(Point p) { translateThis(-p.x(), -p.y(), -p.z()); }
	
	/** Scale this point by a factor s */
	public void scaleThis(double s){
		this.coords[0]*=s;
		this.coords[1]*=s;
		this.coords[2]*=s;
	}
	
	/** Returns p added to this (changing this object). */
	public Point addThis(Vector p) { translateThis(p.x(),p.y(),p.z()); return this; }

	public Point addThis(Point p) { translateThis(p.x(), p.y(), p.z()); return this; }
	
	/** Returns a point translated from this one by p (without changing this object). */
	public Point add(Vector p) { return new Point(coords[0]+p.x(), coords[1]+p.y(), coords[2]+p.z()); }

	/** Returns a point translated from this one by (x,y,z) */
	public Point add(double x, double y, double z) { return new Point(coords[0]+x, coords[1]+y, coords[2]+z); }
	public Point add(Point p) { return new Point(coords[0]+p.x(), coords[1]+p.y(), coords[2]+p.z()); }
	
	/** Returns p subtracted from this (changing this object). */
	public Point subtractThis(Vector p) { coords[0] -= p.x(); coords[1] -= p.y(); coords[2] -= p.z(); return this; }
	public Point subtractThis(Point  p) { coords[0] -= p.x(); coords[1] -= p.y(); coords[2] -= p.z(); return this; }

	/** Returns p subtracted from this (without changing this object). */
	public Point subtract(Vector p) {return new Point(coords[0]-p.x(),coords[1]-p.y(),coords[2]-p.z());	}
	public Point subtract(Point p) {return new Point(coords[0]-p.x(),coords[1]-p.y(),coords[2]-p.z());	}
	
	/** Reflects this point through origo. */
	public Point reflectThroughOrigoThis() { coords[0]*=-1; coords[1]*=-1; coords[2]*=-1; return this; }

	/** rotates (clockwise) the point around the line through the origo with the direction unit vector v.
	 * For counterclockwise rotation change signs within parentheses in non-diagonal terms. */
	public void rotationCW(Vector v, double alpha) {
		double c = Math.cos(alpha);
		double d = 1.0-c;
		double s = Math.sin(alpha);
		double vxyd = v.x()*v.y()*d, vxzd = v.x()*v.z()*d, vyzd = v.y()*v.z()*d;
		double vxs = v.x()*s, vys = v.y()*s, vzs = v.z()*s; 
		double xNew = (v.x()*v.x()*d+c)*coords[0] + (vxyd-vzs)*coords[1] + (vxzd+vys)*coords[2];
		double yNew = (vxyd+vzs)*coords[0] + (v.y()*v.y()*d+c)*coords[1] + (vyzd-vxs)*coords[2];
		double zNew = (vxzd-vys)*coords[0] + (vyzd+vxs)*coords[1] + (v.z()*v.z()*d+c)*coords[2];
		setX(xNew);
		setY(yNew);
		setZ(zNew);
	}
	// Daisy
	/* rotates (counter-clockwise) the point around the line through the origo with the direction unit vector v. */
	public void rotationCCW(Vector v, double alpha) {
		double c = Math.cos(alpha);
		double d = 1.0-c;
		double s = Math.sin(alpha);
		double vxyd = v.x()*v.y()*d, vxzd = v.x()*v.z()*d, vyzd = v.y()*v.z()*d;
		double vxs = v.x()*s, vys = v.y()*s, vzs = v.z()*s; 
		double xNew = (v.x()*v.x()*d+c)*x() + -(vxyd-vzs)*y()    + -(vxzd+vys)*z();
		double yNew = -(vxyd+vzs)*x()    + (v.y()*v.y()*d+c)*y() + -(vyzd-vxs)*z();
		setZ(-(vxzd-vys)*x()    + -(vyzd+vxs)*y()    + (v.z()*v.z()*d+c)*z());
		setX(xNew);
		setY(yNew);
	}

	/** rotates (clockwise) the point around the line through the  point p with the direction unit vector v.
	 * For counterclockwise rotation change signs within parentheses in non-diagonal terms. */
	public void rotation(Vector v, double alpha, Point p) {
		this.translateThis(-p.x(), -p.y(), -p.z());
		rotationCW(v, alpha);
		this.translateThis(p.x(), p.y(), p.z());
	}

	/** Returns the sinus of the polar angle of this point with the z-axis*/
	public double polarAngleSinZ() { return coords[1]/distance(); }

	/** Returns the cosinus of the polar angle of this point with the z-axis*/
	public double polarAngleCosZ() { return coords[0]/distance(); }


	
	/** Get the squared distance from this point to point q. */
	public double distanceSquared(Point q) {
		double dx = coords[0]-q.coords[0];
		double dy = coords[1]-q.coords[1];
		double dz = coords[2]-q.coords[2];
		return dx*dx+dy*dy+dz*dz;
	}

	/** Get the distance from this point to point q */
	public double distance(Point q) { 
		double dx = coords[0]-q.coords[0];
		double dy = coords[1]-q.coords[1];
		double dz = coords[2]-q.coords[2];
		return Math.sqrt(dx*dx+dy*dy+dz*dz); 
	}

	public double dot(Point p) { return x()*p.x() + y()*p.y() + z()*p.z(); }
	public double dot(Vector v) { return x()*v.x() + y()*v.y() + z()*v.z(); }
	
	/** Creates a bisector between points p and q */
	public static Plane getBisector(Point p, Point q) {
		if (!p.equals(q)) 
			return new Plane(getMidpoint(p, q), p.vectorTo(q).normalizeThis()); 
		else return null;
	}
	
	/** Creates the midpoint of two points. */
	public static Point getMidpoint(Point p, Point q) {
		return new Point( (p.coords[0] + q.coords[0])/2,(p.coords[1] + q.coords[1])/2,(p.coords[2] + q.coords[2])/2 );		
	}


	// ANGLE METHODS

	/** Get the angle between the line segments p2->p1 and p2->p3. */
	public static double getAngle(Point p1, Point p2, Point p3) {
		return p2.vectorTo(p1).angle(p2.vectorTo(p3));
	}

	/** Get the dihedral angle defined by the 4 non-collinear points p1, p2, p3, p4. */
	public static double getDihedralAngle(Point p1, Point p2, Point p3, Point p4) {
		return Vector.getDihedralAngle(p1.vectorTo(p2), p2.vectorTo(p3), p3.vectorTo(p4));
	}

  	/*
  	 * returns cosinus of the dihedral angle between 4 non-collinear points p1, p2, p3, p4
  	 * added by pawel 12-11-2011
  	 */
  	public static double getCosDihedralAngle(Point p1, Point p2, Point p3, Point p4) {
 		return Vector.getCosDihedralAngle(new Vector(p1, p2), new Vector(p2, p3), new Vector(p3, p4));
  	}

	
	// COMPARISON METHODS

	/**
	 * Returns true if this point dominates point q. One point is said to dominate another 
	 * if it has a higher x-coordinate. If two points have identical x-coordinates, the 
	 * y-coordinate is considered and so forth.
	 */
	public boolean dominates(Point q) { 
		if (coords[0] > q.coords[0]) return true;
		if (coords[0] < q.coords[0]) return false;
		if (coords[1] > q.coords[1]) return true;
		if (coords[1] < q.coords[1]) return false; 
		return coords[2] > q.coords[2];
	}

	/**
	 * Returns true if this point dominates point q (i=0,1,2 is the most important coordinate,
	 * j=0,1,2 is the second most important coordinate and k=0,1,2 is the least important coordinate).
	 */
	public boolean dominates(Point q, int i, int j, int k) {
		if(i==j || i==k || j==k) 
			throw new Error(String.format("i, j and k must be distinct coordinate indices (%d,%d,%d)",i,j,k));
		if (this.getCoord(i) > q.getCoord(i)) return true;
		if (this.getCoord(i) < q.getCoord(i)) return false;
		if (this.getCoord(j) > q.getCoord(j)) return true;
		if (this.getCoord(j) < q.getCoord(j)) return false; 
		return this.getCoord(k) > q.getCoord(k);
	}

	/** 
	 * Returns a clone of this point. Since a point can be interpreted as a geometric shape 
	 * (a 0-simplex) the Shape interface requires the getCenter method to be implemented.
	 * TODO: Test
	 */
	public Point getCenter() {
		return clone();
	}
	
	public Point getPoint(int i){
		if(i!=0) throw new IllegalArgumentException("Invalid index ("+i+") 0-simplex has one point only");
		return this;
	}
	
	/** Returns true iff o is a point that equals this point. */
	public boolean equals(Object o){
		if(o instanceof Point) return equals((Point)o);
		else return false;
	}
	
	/** Returns true iff this point and point p are overlapping. */
	public boolean equals(Point p) {
		if(Math.abs(coords[0]-p.coords[0])>Constants.EPSILON) return false;
		if(Math.abs(coords[1]-p.coords[1])>Constants.EPSILON) return false;
		if(Math.abs(coords[2]-p.coords[2])>Constants.EPSILON) return false;
		return true;
	}
	
	public static Point getCircumCenter(Point a, Point b, Point c) {
		/*Plane p0 = Point.getBisector(a, b);
		Plane p1 = Point.getBisector(a, c);
		Line l = p0.getIntersection(p1);
		Plane p = new Plane(a, b, c);
		return p.getIntersection(l);*/
		
		Vector ca = new Vector(c, a);
		Vector cb = new Vector(c, b);
		Vector cr = ca.cross(cb);
		Vector v1 = cb.multiply(ca.getLengthSquared());
		Vector v2 = ca.multiply(cb.getLengthSquared());
		v1.subtractThis(v2);
		v1.crossThis(cr);
		v1.divideThis(2.0*cr.getLengthSquared());
		return c.add(v1);
	}
	
	
	/* Returns equilateral point of 2 points a and b in the plane through a, b and c */
	public static Point getEquilateralPoint(Point a, Point b, Point c) {
		Point e = a.clone();
		Vector ba = new Vector(b, a);
		Vector normal = ba.cross(new Vector(b,c));
		normal.normalizeThis();
		e.rotation(normal, -Math.PI/3,b);	
		return e;
	}
 
	/* Returns the equilateral circle of a and b in the plane through a, b, and c */
	public static Circle getEquilateralCircle(Point a, Point b, Point c) {
		return new Circle(a, b, Point.getEquilateralPoint(a, b, c));
	}
	
	public static Circle getEquilateralPoints(Point a, Point b) {
		Point center = Point.getMidpoint(a, b);
		double radius = Math.sqrt(3)*a.distance(b)/2;
		Vector normal = new Vector(a,b).normalizeThis();
		return new Circle(center, radius, normal);
	}
	
	/* Returns Steiner point of 3 points */
	public static Point getSteinerPoint(Point a, Point b, Point c) {
		// if the angle at a is 120 or more, return a. Similarly for b and c
		Vector ab = new Vector(a,b);
		Vector ba = new Vector(b,a);
		Vector bc = new Vector(b,c);
		Vector cb = new Vector(c,b);
		Vector ca = new Vector(c,a);
		Vector ac = new Vector(a,c);

		if (!ab.isSteinerAngle(ac)) return a;
		if (!ba.isSteinerAngle(bc)) return b;
		if (!ca.isSteinerAngle(cb)) return c;
		Point eab = Point.getEquilateralPoint(a, b, c);
		Point center = Point.getCircumCenter(a, b, eab);
		Vector normal = ab.cross(ac).normalize();
		Sphere sphere = new Sphere(new Circle(center, a, normal));
		LineSegment sgm =  sphere.getIntersection(new Line(eab, c));
		if (sgm.a.distanceSquared(eab) > sgm.b.distanceSquared(eab)) return sgm.a; else return sgm.b;
	}

	/** Return a new object that equals this object. */
	public Point clone(){ return new Point(coords[0], coords[1], coords[2]); }
	
	/** Swaps points a and b */
	public static void swap(Point a, Point b) {
		Point temp = a;
		a = b;
		b = temp;
	}
	
	/** Returns the vector from origo to this point. Converts this point to a vector. */
	public Vector toVector() { return new Vector(coords[0], coords[1], coords[2]); }

	/** Returns a string-representation of this point formatted with two decimals precision. */ 
	public String toString() {	return toString(2); }

	/** Returns a string-representation of this point formatted with <code>dec</code> decimals precision. */ 
	public String toString(int dec) {
		return String.format("Point[%."+dec+"f,%."+dec+"f,%."+dec+"f]", coords[0], coords[1], coords[2]); 
	}	
	
	public Sphere toScene(J3DScene scene, double r, Color clr) { 
		Sphere sph = new Sphere(this, r);
		scene.addShape(sph, clr);
		return sph;
	}

	public Sphere toScene(J3DScene scene, double r, Color clr, int divisions) { 
		Sphere sph = new Sphere(this, r);
		scene.addShape(sph, clr, divisions);
		return sph;
	}

//	public void draw(J3DScene scene, double r) { draw(scene, r, Color.BLUE); }
//	public void draw(J3DScene scene) { draw(scene,0.1f, Color.BLUE); }


	/** Writes this point to <code>System.out</code>. */
	public void toConsole() { System.out.println(toString()); }
	
	/** Writes this point to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	public static void main(String[] args) {
		
		J3DScene scene  = J3DScene.createJ3DSceneInFrame();
		Point a = new Point(-1, 0, 0);
		Point b = new Point(1, 0, 0);
		a.toScene(scene, 0.03, Color.blue);
		b.toScene(scene, 0.03, Color.blue);
		Circle cab = Point.getEquilateralPoints(a, b);
		cab.toScene(scene, 0.01, Color.blue);
		Point e = cab.getPoint();
		for (int i = 0; i < 36; i++) {
			e.rotation(cab.getNormal(), Math.PI/18, cab.getCenter());
			Circle eqCircle = new Circle(a, b, e);
			eqCircle.toSceneArc(scene, 0.01, Color.blue, 120, a);
		}
		Point c = new Point(0, 1, 5);
		Point d = new Point(0, 2, 4);
		c.toScene(scene, 0.03, Color.red);
		d.toScene(scene, 0.03, Color.red);
		Circle ccd = Point.getEquilateralPoints(c, d);
		ccd.toScene(scene, 0.01, Color.red);
		Point f = ccd.getPoint();
		for (int i = 0; i < 36; i++) {
			f.rotation(ccd.getNormal(), Math.PI/18, ccd.getCenter());
			Circle fqCircle = new Circle(c, d, f);
			fqCircle.toSceneArc(scene, 0.01, Color.red, 120, c);
			e = cab.getFarthestPoint(f);
			LineSegment ef = new LineSegment(e, f);
			ef.toScene(scene, 0.01, Color.cyan);
		}
	}
	
}



