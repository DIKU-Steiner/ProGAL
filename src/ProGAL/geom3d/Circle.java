package ProGAL.geom3d;



import java.awt.Color;

import ProGAL.Function;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Cylinder;
import ProGAL.math.Constants;
import ProGAL.math.Functions;
import ProGAL.math.RootFinding;

/*
 * If the circle has center c, radius r, and unit-length normal vector n, compute unit-length vectors 
u and v so that {u,v,n} are mutually orthogonal. The circle is then parameterized by P(t) = c + rcos(t)u + rsin(t)v for 0 <= t < 2*pi. 
 */

/** 
 * A circle in (x,y,z)-space represented by a center-point, a radius and a normal-vector.
 */
public class Circle implements Shape{
	private Point center;
	private double radius;
	private Vector normal;

	public Circle(Point center, double radius, Vector normal) {
		this.center = center;
		this.radius = radius;
		if (normal!=null) this.normal = normal.normalize();
	}
	
	/** A circle with given center through a given point  and with specified normal vector */
	public Circle(Point center, Point through, Vector normal) {
		this.center = center;
		radius = center.distance(through);
		this.normal = normal.normalize();
	}
	
	/** Circle in the plane through p0, p1, p2 */ 
    public Circle(Point p0, Point p1, Point p2) {
	  center = Point.getCircumCenter(p0,  p1, p2);
	  radius = center.distance(p0);
	  Vector v0 = new Vector(center, p0);
	  Vector v1 = new Vector(center, p1);
	  normal = v0.cross(v1).normalize();
	}		

	
	/*
	 * Given three points p0, p1, p2 and a vector v, find the circle or just the radius of the circle through the 
	 * projection of the points onto the plane with v as its normal vector. Claim: radius is minimized when of the 
	 * plane goes through p0, p1, p2
	 */
	/*public Circle3d(Point3d p0, Point3d p1, Point3d p2, Vector3d v) {
		// create the plane through the origo with v as its normal vector
		Plane3d plane = new Plane3d(v, new Point3d(0,0,0));
		Point3d q0 = plane.projectPoint(p0);
		Point3d q1 = plane.projectPoint(p1);
		Point3d q2 = plane.projectPoint(p2);
		double a = q0.getSquaredDistance(q1);
		double b = q1.getSquaredDistance(q2);
		double c = q2.getSquaredDistance(q0);
		radius = Math.sqrt(a*b*c)/Math.sqrt(2*a*b + 2*b*c + 2*c*a - a*a - b*b - c*c);
		normalVector = Vector3d.crossProduct(new Vector3d(p0,p1), new Vector3d(p0,p2));
	}*/
	

	/*
	 * returns the radius of the circle through 3 given points (without creating the circle)
	 */
	/*public static double getRadius(Point3d p0, Point3d p1, Point3d p2) {

		// get the plane through q0, q1, q2
		
		Point3d q0 = new Point3d(0,0,0);
		Point3d q1 = new Point3d(p1.x-p0.x, p1.y-p0.y, p1.z-p0.z);
		Point3d q2 = new Point3d(p2.x-p0.x, p2.y-p0.y, p2.z-p0.z);
		
		// get the plane through q0, q1, q2
		
		Plane3d plane = new Plane3d(q0,q1,q2);
		Vector3d verticalVector = new Vector3d(0,0,1);
		double angle = Vector3d.getAngle(plane.n, verticalVector);
		Vector3d rotationVector = Vector3d.crossProduct(plane.n, verticalVector);
		rotationVector.normalizeThis();
		q1.rotation(rotationVector, angle);
		q2.rotation(rotationVector, angle);
		Point2d r0 = new Point2d(0,0);
		Point2d r1 = new Point2d(q1.x, q1.y);
		Point2d r2 = new Point2d(q2.x, q2.y);
		Circle2d circle2 = new Circle2d(r0, r1, r2);
		return circle2.getRadius();
	}*/
	
	/** Get the center of the circle. */
	public Point getCenter() { return center; }
	/** Get the radius of the circle. */
	public double getRadius() { return radius; }
	/** Get the normal of the circle. */
	public Vector getNormalVector() { return normal; }
	public Vector getNormal() { return normal; }
	/** Sets the center of the circle */
	public void setCenter(Point p) { center = p; }
	
	/** return a point on the circle */
	public Point getPoint() {
		return center.add(normal.getOrthonormal().scaleToLength(radius));		
	}
	
	/** returns plane through this circle */
	public Plane getPlane() {
		return new Plane(center, normal);
	}
	
	/** returns the point on the circle closets to the query point p. If p is equidistant to all points on the circle, 
	 * an arbitrary point of the circle is returned.
	 */
	public Point getClosestPoint(Point p) {
		Point pr = new Plane(center, normal).projectPoint(p);
		if (center.distanceSquared(pr) <= Constants.EPSILON) return getPoint();
		else return center.add(new Vector(center, pr).scaleToLength(radius));
	}

	/** returns the point on the circle farthest from the query point p. If p is equidistant to all points on the circle, 
	 * an arbitrary point of the circle is returned.
	 */
	public Point getFarthestPoint(Point p) {
		Point pr = new Plane(center, normal).projectPoint(p);
		if (center.distanceSquared(pr) <= Constants.EPSILON) return getPoint();
		else return center.add(new Vector(center, pr).scaleToLength(radius).multiply(-1));
	}

	public double getClosestDistance(Circle cb) {
		Vector u = cb.normal.getOrthonormal();
		Vector v = cb.normal.cross(u);
		double ra2 = radius*radius;
		double ca2 = center.dot(center);
		double cacb = center.dot(cb.center);
		double cau = center.dot(u);
		double cav = center.dot(v);
		double naca = normal.dot(center);
		double nacb = normal.dot(cb.center);
		double nau = normal.dot(u);
		double nav = normal.dot(v);
		double rb2 = cb.radius*cb.radius;
		double cb2 = cb.center.dot(cb.center);
		double cbu = cb.center.dot(u);
		double cbv = cb.center.dot(v);
		double a0 = rb2 + cb2 + ca2 - 2*cacb - (nacb - naca)*(nacb - naca);
		double a1 = -2*rb2*nau*nav;
		double a2 = 2*cb.radius*(cbv - cav - (nacb - naca)*nav);
		double a3 = 2*cb.radius*(cbu - cau - (nacb - naca)*nau);
		double a4 = -rb2*nav*nav;
		double a5 = -rb2*nau*nau;
		double a6 = -2*radius;
		double a7 = ra2 + rb2 + ca2 + cb2 - 2*cacb;
		double a8 = 2*cb.radius*(cbv - cav);
		double a9 = 2*cb.radius*(cbu - cau);
		double[] c = new double[5];
		c[0] = a9*a9-a6*a6*a5 + 2*a7*a9 - a6*a6*a3 + a7*a7 - a6*a6*a0;
		c[1] = 2*(2*a7*a8 - a6*a6*a2) + 2*(2*a8*a9 - a6*a6*a1);
		c[2] = 2*(a6*a6*a5 - a9*a9) + 4*(a8*a8 -a6*a6*a4) + 2*(a7*a7 - a6*a6*a0);
		c[3] = 2*(2*a7*a8 - a6*a6*a2) - 2*(2*a8*a9 - a6*a6*a1);
		c[4] = a9*a9-a6*a6*a5 - 2*a7*a9 + a6*a6*a3 + a7*a7 -a6*a6*a0;
		Function f = new Function(c);
		double root = RootFinding.brent(f, 0.5, 0.8, 0.000001, 0.00001, 1000);
		return root;
	}
	
	
	/** Create the equilateral circle of two points. */
	public static Circle getEquilateralCircle(Point a, Point b) {
		Point center = Point.getMidpoint(a, b);
		Vector ab = a.vectorTo(b);
		double radius = Math.sqrt(3)*ab.length()/2;
		return new Circle(center, radius, ab);
	}
	

	/** Returns a string-representation of this circle formatted with two decimals precision. */
	public String toString(){
		return toString(2);
	}
	
	/** Returns a string-representation of this circle formatted with <code>dec</code> decimals precision. */
	public String toString(int dec) {
		return String.format("Circle[center=%s,radius=%"+dec+"f,normal=%s]",center.toString(dec),radius,normal.toString(dec));
	}

	/** Writes this circle to <code>System.out</code> with 2 decimals precision. */
	public void toConsole() { toConsole(2); }

	/** Writes this circle to <code>System.out</code> with <code>dec</code> decimals precision. */
	public void toConsole(int dec) {
		System.out.println(toString(dec)); 
	}
	
 
	/** Intersection of 2 circles in the same plane */
	public Point[] getIntersection(Circle c) {
		Point[] intersectionPoints = null;
		double dist = center.distance(c.getCenter());
		double sum = radius + c.getRadius();
		if (dist - sum > Constants.EPSILON) return intersectionPoints;
		else {
			double diff = Math.abs(radius - c.getRadius());
			if (dist - diff < -Constants.EPSILON) return intersectionPoints;
			else {
				if (Math.abs(dist - sum) < Constants.EPSILON) { // one intersection point
					Vector v = new Vector(center, c.getCenter());
					double fraction = radius/(radius+c.getRadius());
					intersectionPoints = new Point[1];
					intersectionPoints[0] = center.addThis(v.multiply(fraction));
//					System.out.println("Scale vector correct in Circle class? "+(center.distance(center.addThis(v.multiply(fraction)))==radius));
					return intersectionPoints;
				}
				else {
					if (dist - diff < Constants.EPSILON) {
						Vector v;
						intersectionPoints = new Point[1];
						if (radius > c.radius) {
							v = new Vector(center, c.getCenter());
							v.scaleToLengthThis(radius);
							intersectionPoints[0] = center.add(v);
						}
						else {
							v = new Vector(c.getCenter(), center);
							v.scaleToLengthThis(c.getRadius());
							intersectionPoints[0] = c.getCenter().add(v);
						}
					}
					else {
						double alpha = Math.acos((radius*radius + dist*dist - c.getRadius()*c.getRadius())/(2*radius*dist));
						Vector v = new Vector(center, c.getCenter());
						v.scaleToLengthThis(radius);
						intersectionPoints = new Point[2];
						getNormalVector().rotateIn(v, alpha);
						intersectionPoints[0] = center.add(v);
						getNormalVector().rotateIn(v, -2*alpha);
						intersectionPoints[1] = center.add(v);
					}
				}
			}
		}
		return intersectionPoints;
	}
	
	/** returns the smallest rotation angle (direction determined by vector dir) 
	 * needed to bring point p on this circle to be on the line l as well. Returns null if 
	 * there is no intersection or just one intersection.
	 * @return
	 */
	public Double getFirstIntersection(Line line, Point p, Vector dir, J3DScene scene) {
		double dist = line.getDistance(center);
		if (dist > radius - Constants.EPSILON) return null;
		Vector op = new Vector(center, p);
		double r2 = radius*radius;
		if (dist < 0.001) {
			Vector oco = line.dir.multiply(radius);
			double cosBeta = op.dot(oco)/r2;
			Vector cr  = op.cross(oco);
			cr.divideThis(r2);
			double sinBeta = cr.length();
			double beta = Math.atan2(sinBeta,cosBeta); System.out.println("beta = " + Functions.toDeg(beta));
			if (normal.dot(dir) > 0) {
				if (normal.dot(cr) > 0) return beta; else return 2*Math.PI - beta;
			}
			else { 
				if (normal.dot(cr) > 0) return 2*Math.PI - beta; else return beta;
			}
		}
		Point co = line.orthogonalProjection(center);
		double oco2 = dist*dist;
		double pco2 = co.distanceSquared(p);
//		double cosAlpha = Math.sqrt(1-oco2/r2); 
		double cosAlpha = dist/radius;
		double cosBeta = (r2 + oco2 - pco2)/(2*radius*dist);
		Vector oco = new Vector(center, co);
//		LineSegment seg1 = new LineSegment(center, co);      // seg1.toScene(scene,0.03, Color.red);
//		LineSegment seg2 = new LineSegment(center, p);       // seg2.toScene(scene, 0.03, Color.pink);
		Vector cr = op.cross(oco);
		cr.divideThis(op.length()*oco.length());
		double sinBeta = cr.length();
		double alpha = Math.acos(cosAlpha); System.out.println("alpha = " + Functions.toDeg(cosAlpha));
		double beta = Math.atan2(sinBeta,cosBeta); System.out.println("beta = " + Functions.toDeg(beta));
		if (normal.dot(dir) > 0) {
			if (normal.dot(cr) > 0) return beta-alpha; else return 2*Math.PI - beta - alpha;
		}
		else { 
			if (normal.dot(cr) > 0) return 2*Math.PI - beta - alpha; else return beta - alpha;
		}
	}
	
	//Daisy
	public double getDistanceSquared(Point p) {
		double NPC = normal.dot(p.subtract(center));
		double secondTerm = Math.sqrt(p.distanceSquared(center)-NPC*NPC)-radius;
		return NPC*NPC+secondTerm*secondTerm;
	}
	//Daisy
	public double getDistance(Point p) {
		return Math.sqrt(getDistanceSquared(p));
	}
	
	/** returns the smallest rotation angle (direction determined by vector dir) 
	 * needed to bring point p on this circle to be on the circle c as well. Returns null if 
	 * there is no intersection or just one intersection.
	 * @return
	 */
	public Double getFirstIntersection(Circle c, Point p, Vector dir) {
		double dist = center.distance(c.getCenter()); 
		if (dist > radius + c.getRadius() - Constants.EPSILON) return null;
		double r2d2 = radius*radius + dist*dist;
		double rd2 = 2*radius*dist;
		double cosAlpha = (r2d2- c.getRadius()*c.getRadius())/rd2; System.out.println("cosAlpha =" + cosAlpha + " " + Math.acos(cosAlpha));
		double distP2 = p.distanceSquared(c.getCenter());
		double cosBeta = (r2d2 - distP2)/rd2; System.out.println("cosBeta =" + cosAlpha + " " + Math.acos(cosBeta));
		Vector op = new Vector(center, p);
		Vector oco = new Vector(center, c.getCenter());
		Vector cr = op.cross(oco);
		cr.divideThis(op.length()*oco.length());
		double sinBeta = cr.length();
		double alpha = Math.acos(cosAlpha); System.out.println("alpha = " + Functions.toDeg(cosAlpha));
		double beta = Math.atan2(sinBeta,cosBeta); System.out.println("beta = " + Functions.toDeg(beta));
		if (normal.dot(dir) > 0) {
			if (normal.dot(cr) > 0) return beta-alpha; else return 2*Math.PI - beta - alpha;
		}
		else { 
			if (normal.dot(cr) > 0) return 2*Math.PI - beta - alpha; else return beta - alpha;
		}
	}
	
	/** Draws the circle as 360 dots */
	public void toScene(J3DScene scene, double width, Color clr) {
		double step = Math.PI/180;
		double cosI;
		double sinI;
		double angle = -step;
		Vector u = normal.getOrthonormal().scaleToLength(radius);
		Vector nxu = normal.cross(u);
		for (int i = 0; i < 360; i++) {
			angle += step;
			cosI = Math.cos(angle);
			sinI = Math.sin(angle);
			new Point(u.x()*cosI + nxu.x()*sinI + center.x(),
					  u.y()*cosI + nxu.y()*sinI + center.y(),
					  u.z()*cosI + nxu.z()*sinI + center.z()).toScene(scene, width, clr);
		}
	}
	
	/** Draws an arc as series of dots 1 degree apart*/
	public void toSceneArc(J3DScene scene, double width, Color clr, int degree, Point start) {
		double step = Math.PI/180;
		double cosI;
		double sinI;
		double angle = -step;
		Vector u = new Vector(center, start).scaleToLength(radius);
		Vector nxu = normal.cross(u);
		for (int i = 0; i < degree; i++) {
			angle += step;
			cosI = Math.cos(angle);
			sinI = Math.sin(angle);
			new Point(u.x()*cosI + nxu.x()*sinI + center.x(),
					  u.y()*cosI + nxu.y()*sinI + center.y(),
					  u.z()*cosI + nxu.z()*sinI + center.z()).toScene(scene, width, clr);
		}
	}
	
	
	public Cylinder[] toScene(J3DScene scene, double width, int res, Color clr) {
		double step = 2*Math.PI/res;
		
		Vector a = normal.getOrthonormal();
		Vector b = a.cross(normal);
		double alpha = 0.0;
		Point p = new Point(center.x() + radius*a.x(), center.y() + radius*a.y(), center.z() + radius*a.z());
		Point q;
		LineSegment seg;
		Cylinder cyl[] = new Cylinder[res];
		for (int i = 0; i < res; i++) {
			alpha += step; 
			double cosAlpha = Math.cos(alpha);
			double sinAlpha = Math.sin(alpha);
			q = new Point(center.x() + radius*(a.x()*cosAlpha + b.x()*sinAlpha),
					      center.y() + radius*(a.y()*cosAlpha + b.y()*sinAlpha),
					      center.z() + radius*(a.z()*cosAlpha + b.z()*sinAlpha));
			seg = new LineSegment(p,q);
			cyl[i] = seg.toScene(scene, width, clr);
			p = q;
		}
		return cyl;
	}
	
	public Cylinder[] toScene(J3DScene scene, double width, int res) {
		return toScene(scene, width, res, Color.blue);
	}
	
	public void fromScene(J3DScene scene, Cylinder[] cyl) {
		for (int i = 0; i < cyl.length; i++) scene.removeShape(cyl[i]);
	}
	
//	private static void intersectionInPlane() {
//		Circle c1 = new Circle(new Point(0,0,0), 0.5, new Vector(0,0,1));
//		Circle c2 = new Circle(new Point(-0.1,0.3,0), 0.6, new Vector(0, 0, 1));
//		J3DScene scene = J3DScene.createJ3DSceneInFrame();
//		c1.toScene(scene, 0.005, Color.blue);
//		c2.toScene(scene, 0.005, Color.red);
//		new LineSegment(c1.center, c2.center).toScene(scene, 0.002, Color.black);
//		double cosg = (c1.center.distanceSquared(c2.center) + c1.radius*c1.radius - c2.radius*c2.radius)/(2*c1.radius*c1.center.distance(c2.center));
//		double alpha = Functions.toDeg(Math.acos(cosg));
//		Vector cc = new Vector(c1.center, c2.center).scaleToLength(c1.radius);
//		cc = c1.normal.rotateIn(cc, Math.acos(cosg));
//		Point q = c1.center.add(cc);
//		q.toScene(scene, 0.03, Color.green);
//		scene.addText("q", q.add(-0.075, 0,0));
//		new LineSegment(c1.center, q).toScene(scene, 0.002, Color.black);
//		new LineSegment(c2.center, q).toScene(scene, 0.002, Color.black);
//		cc = c1.normal.rotateIn(cc, -2*Math.acos(cosg));
//		q = c1.center.add(cc);
//		q.toScene(scene, 0.03, Color.green);
//		scene.addText("c1", c1.center.add(0,-0.075,0));
//		scene.addText("c2", c2.center.add(0,-0.075,0));
//	
//	}
	
	
	@SuppressWarnings("unused")
	private static void intersectionsInSpace() {
		Circle c1 = new Circle(new Point(1,1,0), 0.5, new Vector(0,0,1));
		Circle c2 = new Circle(new Point(-0.1,-0.3,-0.4), 0.6, new Vector(1,0-5,0.5));
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		scene.autoZoom();
		c1.toScene(scene, 0.005, Color.blue);
		c2.toScene(scene, 0.005, Color.red);
		Plane pl1 = c1.getPlane();
		Plane pl2 = c2.getPlane();
		pl1.toScene(scene, new Color(0,0,0,10), 2);
		pl2.toScene(scene, new Color(0,0,0,10), 2);
		Line ln = pl1.getIntersection(pl2);
		ln.toScene(scene, 0.005, Color.black);
		new LineSegment(c1.center, c2.center).toScene(scene, 0.002, Color.black);

	}
	
	public static LineSegment getFurthestDistance_centers(Circle c1, Circle c2) {
		Vector v = c1.center.vectorTo(c2.center);
		Vector y2 = c2.normal.cross(v);
		Vector x2 = c2.normal.cross(y2); x2.multiplyThis(c2.radius/x2.length());
		Point p21 = c2.center.add(x2);
		Point p22 = c2.center.subtract(x2);

		Vector y1 = c1.normal.cross(v).normalizeThis();
		Vector x1 = c1.normal.cross(y1);  x1.multiplyThis(c1.radius/x1.length());
		Point p11 = c1.center.add(x1);
		Point p12 = c1.center.subtract(x1);

		double d11_21 = p11.distanceSquared(p21);
		double d12_21 = p12.distanceSquared(p21);
		double d11_22 = p11.distanceSquared(p22);
		double d12_22 = p12.distanceSquared(p22);
		if(d11_21>d12_21 && d11_21>d11_22 && d11_21>d12_22) return new LineSegment(p11, p21);
		if(d12_21>d11_22 && d12_21>d12_22) return new LineSegment(p12, p21);
		if(d11_22>d12_22) return new LineSegment(p11, p22);
		return new LineSegment(p12, p22);
	}
	
	public static Point getFurthestPoint_bruteForce(Circle c, Point p){
		if(c.radius<Constants.EPSILON)
			return c.center;
		
		int divisions = 8;
		Vector x = new Vector(1.001,1.002,1.003).crossThis(c.normal).normalizeThis();
		Vector y = c.normal.cross(x);
		
		double range = 2*Math.PI;
		double maxT = 0; 
		double maxDist = 0;
		while(range>Math.PI/10){
			for(double t=maxT-range/2;t<maxT+range/2;t+=range/divisions){
				Point pc = c.center.add(  x.multiply(c.radius*Math.cos(t)).addThis(y.multiply(c.radius*Math.sin(t)))  );
				double dist = p.distanceSquared(pc);
				if(dist>maxDist){
					maxT = t;
					maxDist = dist;
				}
			}
			range/=8;
		}
		Point pc = c.center.add(  x.multiply(c.radius*Math.cos(maxT)).addThis(y.multiply(c.radius*Math.sin(maxT)))  );
		return pc;
	}

	public static LineSegment getFurthestDistance_bruteForce(Circle c1, Circle c2) {
		if(c1.radius<Constants.EPSILON) 
			return new LineSegment(c1.center, getFurthestPoint_bruteForce(c2, c1.center));
		
		int divisions = 8;
		Vector x1 = new Vector(1.001,1.002,1.003).crossThis(c1.normal).normalizeThis();
		Vector y1 = c1.normal.cross(x1);
		Vector x2 = new Vector(1.001,1.002,1.003).crossThis(c2.normal).normalizeThis();
		Vector y2 = c2.normal.cross(x2);
		
		double range = 2*Math.PI;
		double maxT1 = 0; 
		double maxT2 = 0;
		double maxDist = 0;
		while(range>Math.PI/10){
			for(double t1=maxT1-range/2;t1<maxT1+range/2;t1+=range/divisions){
				Point p1 = c1.center.add(  x1.multiply(c1.radius*Math.cos(t1)).addThis(y1.multiply(c1.radius*Math.sin(t1)))  );
				for(double t2=maxT2-range/2;t2<maxT2+range/2;t2+=range/divisions){
					Point p2 = c2.center.add(  x2.multiply(c2.radius*Math.cos(t2)).addThis(y2.multiply(c2.radius*Math.sin(t2)))  );
					double dist = p1.distanceSquared(p2);
					if(dist>maxDist){
						maxT1 = t1;
						maxT2 = t2;
						maxDist = dist;
					}
				}
			}
			range/=8;
		}
		Point p1 = c1.center.add(  x1.multiply(c1.radius*Math.cos(maxT1)).addThis(y1.multiply(c1.radius*Math.sin(maxT1)))  );
		Point p2 = c2.center.add(  x2.multiply(c2.radius*Math.cos(maxT2)).addThis(y2.multiply(c2.radius*Math.sin(maxT2)))  );
		
		return new LineSegment(p1, p2);
	}
	
	
	public static void main(String[] args) {
	
		Point p = new Point (4, 1, 4);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		p.toScene(scene,0.03,Color.black);
		Circle c = new Circle(new Point(0, 2, 0), 1, new Vector(1,3,1).normalize());
		c.toScene(scene, 0.02, Color.red);
		Circle cb = new Circle(new Point(3, 2, 1), 1, new Vector(2, 1, 3).normalize());
		cb.toScene(scene, 0.02, Color.blue);
		c.getClosestDistance(cb);

/*		Point q = c.getFarthestPoint(p);
		q.toScene(scene,0.03,Color.blue);
		LineSegment s = new LineSegment(p,q);
		s.toScene(scene, 0.01, Color.blue);
*/	}

}

