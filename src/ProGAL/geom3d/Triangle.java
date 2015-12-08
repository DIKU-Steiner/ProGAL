package ProGAL.geom3d;

import java.awt.Color;

import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;

/**
 * A triangle in (x,y,z)-space represented by the three corner-points.
 */
public class Triangle implements Simplex{
	protected Point p1, p2, p3;
	protected Shape[] LSSs = new Shape[3];
	protected Shape face;


	/** Construct a triangle using the three specified points as corners */
	public Triangle(Point p1, Point p2, Point p3) {
		this.p1 = p1;
		this.p2 = p2;
		this.p3 = p3;
	}

	public Triangle(Point[] p) {
		this.p1 = p[0];
		this.p2 = p[1];
		this.p3 = p[2];
	}

	
	/** Get the first corner */
	public Point getP1(){ return p1; }
	/** Get the second corner */
	public Point getP2(){ return p2; }
	/** Get the third corner */
	public Point getP3(){ return p3; }
	/** Get the specified corner of this triangle */
	public Point getCorner(int c){ return getPoint(c); }
	/** Get the specified corner-point of this triangle */
	public Point getPoint(int c){
		switch(c){
		case 0: return p1;
		case 1: return p2;
		case 2: return p3;
		}
		throw new Error("Badly specified point number ("+c+"). Should be between 0 and 2");
	}

	/** Return the 'dimension' of this object. Required by the interface Simplex. */
	public int getDimension() { return 2; }

	public boolean orient(Point p){
		Matrix m = new Matrix(4,4);
		for(int r=0;r<3;r++){
			for(int c=0;c<3;c++){
				m.set(r, c, getCorner(r).getCoord(c));
			}
			m.set(r, 3, 1);
		}
		for(int c=0;c<3;c++)
			m.set(3, c, p.getCoord(c));
		m.set(3,3,1);

		double det = m.determinant();
		if (Math.abs(det)<Constants.EPSILON) return true;
		return det<0;
	}
	
	/** Return the center of the triangle. Here average of the corners is used.*/
	public Point getCenter() { 
		return new Point( 
			(p1.x()+p2.x()+p3.x())/3, 
			(p1.y()+p2.y()+p3.y())/3, 
			(p1.z()+p2.z()+p3.z())/3
		); 
	}

	/** Return the area of one side of the triangle. */
	public double getArea(){
		return 0.5*p1.vectorTo(p2).crossThis(p1.vectorTo(p3)).length();
	}
	
	/** Return a vector that is normal to this triangle. */
	public Vector getNormal() {
		return p1.vectorTo(p2).crossThis(p1.vectorTo(p3)).normalizeThis();
	}
	
	/** 
	 * Return the circumradius of the triangle. If one side has zero length this method returns 
	 * the length of the two remaining sides.
	 */
	public double circumradius(){
		double a = p1.distance(p2);
		double b = p1.distance(p3);
		double c = p2.distance(p3);
		double s = (a+b+c)/2;//Semiperemiter
		return a*b*c/(4*Math.sqrt(s*(a+b-s)*(a+c-s)*(b+c-s)));
	}

	/** Return the circumcenter of the triangle. 
	 * TODO: Test
	 * TODO: Make more efficient (transform to origo with n as z and use 2D formula)	 
	 */
	public Point circumcenter() {
		Vector n = getNormal();
		Point m1 = Point.getMidpoint(p1, p2);
		Point m2 = Point.getMidpoint(p1, p3);
		Line l1 = new Line(m1, p1.vectorTo(p2).crossThis(n));
		Line l2 = new Line(m2, p1.vectorTo(p3).crossThis(n));
		return l1.getIntersection(l2);
	}

	public double inradius(){
		double a = p1.distance(p2);
		double b = p1.distance(p3);
		double c = p2.distance(p3);
		double s = (a+b+c)/2;//Semiperemiter
		return Math.sqrt( ((s-a)*(s-b)*(s-c))/s );
	}
	
	public Point incenter(){
		double a = p1.distance(p2);
		double b = p1.distance(p3);
		double c = p2.distance(p3);
		double P = a+b+c;
		Vector C = p3.toVector().multiplyThis(a);
		C.addThis(p2.toVector().multiplyThis(b));
		C.addThis(p1.toVector().multiplyThis(c));
		C.divideThis(P);
		return C.toPoint();
	}
	
	public Point getIntersection(Point p, Point q) {
		Vector dir = new Vector(p, q);
		//Daisy
/*		Plane plane = new Plane(p1, p2, p3);
		Line line = new Line(p, q);
		Point intersectionPoint = plane.getIntersection(line);*/
		//Rasmus
		Vector u = new Vector(p1, p2);
		Vector v = new Vector(p1, p3);
		Vector n = u.cross(v);
		if (n.isZeroVector())  { System.out.println("Normal is zero"); return null; }    // triangle is degenerated
		Vector w0 = new Vector(p1, p);
		double a = -n.dot(w0);
		double b = n.dot(dir);
		if (Math.abs(b) < Constants.EPSILON) {
			if (a == 0) { System.out.println("a is zero"); return null; } // ray is in the triangle plane
			else { System.out.println("a is not zero"); return null; } // ray is not intersecting the plane
		}
		double r = a/b;
		if (r < 0.0) { /*System.out.println("r is less than zero");*/ return null; } // ray goes away from the triangle plane, no intersection
		// for a segment, also test if r > 1.0. if so, no intersection
		
		Point intersection = new Point(p.x() + r*dir.x(), p.y() + r*dir.y(), p.z() + r*dir.z());
		if (contains(intersection)) return intersection; else { return null; } 
	}
	
	//Daisy
	public boolean containsPoint(Point p) {
		return (p1.equals(p) || p2.equals(p) || p3.equals(p));
	}
	
	public boolean contains(Point p) {
		//Daisy
/*		double[] r1 = {p1.getCoord(0), p2.getCoord(0), p3.getCoord(0)};
		double[] r2 = {p1.getCoord(1), p2.getCoord(1), p3.getCoord(1)};
		double[] r3 = {1, 1, 1};
		double[][] rows = {r1, r2, r3};
		Matrix orientation = new Matrix(rows);
		if (orientation.determinant()<0) {
			System.out.println("Wooh!");
			System.out.println("Triangle = "+this.toString());
			System.out.println("Point = "+p.toString());
			double[] r1_12 = {p1.getCoord(0), p2.getCoord(0), p.getCoord(0)};
			double[] r2_12 = {p1.getCoord(1), p2.getCoord(1), p.getCoord(1)};
			double[] r3_12 = {1, 1, 1};
			double[][] rows12 = {r1_12, r2_12, r3_12};
			Matrix orientation12 = new Matrix(rows12);
			if (orientation12.determinant()<0) {
				System.out.println("Wooh!");
				double[] r1_23 = {p2.getCoord(0), p3.getCoord(0), p.getCoord(0)};
				double[] r2_23 = {p2.getCoord(1), p3.getCoord(1), p.getCoord(1)};
				double[] r3_23 = {1, 1, 1};
				double[][] rows23 = {r1_23, r2_23, r3_23};
				Matrix orientation23 = new Matrix(rows23);
				if (orientation23.determinant()<0) {
					System.out.println("Wooh!");
					double[] r1_31 = {p3.getCoord(0), p1.getCoord(0), p.getCoord(0)};
					double[] r2_31 = {p3.getCoord(1), p1.getCoord(1), p.getCoord(1)};
					double[] r3_31 = {1, 1, 1};
					double[][] rows31 = {r1_31, r2_31, r3_31};
					Matrix orientation31 = new Matrix(rows31);
					if (orientation31.determinant()<0) { System.out.println("Wooh!"); return true; }
				}
			}
		}
		return false;*/
/*		
		List<Double> xs = new ArrayList<Double>();
		xs.add(p1.get(0)); xs.add(p2.get(0)); xs.add(p3.get(0));
		List<Double> ys = new ArrayList<Double>();
		ys.add(p1.get(1)); ys.add(p2.get(1)); ys.add(p3.get(1));
		List<Double> zs = new ArrayList<Double>();
		zs.add(p1.get(2)); zs.add(p2.get(2)); zs.add(p3.get(2));
		
		double xMax = Collections.max(xs);
		double xMin = Collections.min(xs);
		double yMax = Collections.max(ys);
		double yMin = Collections.min(ys);
		double zMax = Collections.max(zs);
		double zMin = Collections.min(zs);
		if (xMax > p.get(0) && xMin < p.get(0) && yMax > p.get(1) && yMin < p.get(1) && zMax > p.get(2) && zMin < p.get(2)) return true;
		else return false;*/
		
		/*double alpha = ((p2.getCoord(1) - p3.getCoord(1))*(p.getCoord(0) - p3.getCoord(0)) + (p3.getCoord(0) - p2.getCoord(0))*(p.getCoord(1) - p3.getCoord(1))) /
		        ((p2.getCoord(1) - p3.getCoord(1))*(p1.getCoord(0) - p3.getCoord(0)) + (p3.getCoord(0) - p2.getCoord(0))*(p1.getCoord(1) - p3.getCoord(1)));
		System.out.println("Alpha = "+alpha);
		if (alpha>-Constants.EPSILON) {
			double beta = ((p3.getCoord(1) - p1.getCoord(1))*(p.getCoord(0) - p3.getCoord(0)) + (p1.getCoord(0) - p3.getCoord(0))*(p.getCoord(1) - p3.getCoord(1))) /
				       ((p2.getCoord(1) - p3.getCoord(1))*(p1.getCoord(0) - p3.getCoord(0)) + (p3.getCoord(0) - p2.getCoord(0))*(p1.getCoord(1) - p3.getCoord(1)));
			System.out.println("Beta = "+beta);
			if (beta>-Constants.EPSILON) {
				double gamma = 1.0 - alpha - beta;
				if (gamma>-Constants.EPSILON) return true;
				System.out.println("Gamma = "+gamma);
			}
		}
		return false;*/
		//Rasmus
		Vector u = new Vector(p1, p2);
		Vector v = new Vector(p1, p3);
		double uu = u.dot(u);
		double uv = u.dot(v);
		double vv = v.dot(v);
		Vector w = new Vector(p1, p);
		double wu = w.dot(u);
		double wv = w.dot(v);
		double D = uv*uv -uu*vv;
		
		double s = (uv*wv -vv*wu)/D;
		if ((s < 0.0) || (s > 1.0)) return false;
		double t = (uv*wu -uu*wv)/D;
		if ((t < 0.0) || (s + t > 1.0)) return false;
		return true;
	}
	
	/** Returns a string-representation of this triangle formatted with two decimals precision. */
	public String toString(){	return toString(2);	}

	/** Returns a string-representation of this triangle formatted with <code>dec</code> decimals precision. */
	public String toString(int dec) {
		return String.format("Triangle[p1=%s,p2=%s,p3=%s]",p1.toString(dec), p2.toString(dec), p3.toString(dec));
	}

	/** Writes this triangle to <code>System.out</code> with 2 decimals precision. */
	public void toConsole() { toConsole(2); }

	/** Writes this triangle to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) {
		System.out.println(toString(dec)); 
	}
	
	public Tetrahedron toScene(J3DScene scene, Color clr) {
		Tetrahedron t = new Tetrahedron(getCorner(0), getCorner(1), getCorner(2), getCorner(2));
		scene.addShape(t, clr);
		return t;
	}
	
	private boolean isBig(Point p) {
		return ((Math.abs(p.x()) > 1000) || (Math.abs(p.y()) > 1000) || (Math.abs(p.z()) > 1000));
	}

	
	/* draws edges of the triangle.*/
	public void toSceneEdges(J3DScene scene, Color clr, double width) {
		if (LSSs[0] == null) 
			if (!isBig(p1) && !isBig(p2)) LSSs[0] = new LSS(p1, p2, width);
		if (LSSs[1] == null) 
			if (!isBig(p2) && !isBig(p3)) LSSs[1] = new LSS(p2, p3, width);
		if (LSSs[2] == null)
			if (!isBig(p3) && !isBig(p1)) LSSs[2] = new LSS(p3, p1, width);
		for (int k = 0; k < 3; k++) if (LSSs[k] != null) scene.addShape(LSSs[k], clr, 3);
	}

	/* deletes edges of the triangle.*/
	public void fromSceneEdges(J3DScene scene) {
		for (int k = 0; k < 3; k++) scene.removeShape(LSSs[k]);
	}

	public Triangle clone(){
		return new Triangle(p1.clone(), p2.clone(), p3.clone());
	}

 
}

