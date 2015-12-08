package ProGAL.geom3d;

import java.awt.Color;

import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.viewer.TextShape;
import ProGAL.geom3d.volumes.Cylinder;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Constants;
import ProGAL.math.Functions;

/**
 * A plane in (x,y,z)-space represented by a point and a normal. 
 * 
 * Assuming that <i>p</i> is a point on the plane and <i>n</i> is the normal vector, 
 * the half-space <i>{q|pq·n>0}</i> is called the 'upper halfspace' wrt. the plane 
 * and vice versa for the 'lower halfspace'.
 * If the unit normal vector is n = (nx,ny,nz) and the point on the plane is p = (px,py,pz), then the plane has 
 * equation  ax + bx + cx = d  where a = nx, b = ny, c = nz, and d = px*nx + py*ny + pz*nz  
 */
public class Plane implements Shape{
	/** Normal vector of the plane. */
	protected Vector normal; 
	/** Point in the plane. */
	protected Point point;  

	/** Constructs a plane with the normal vector n containing point p. */
	public Plane(Point p, Vector n) {
		this.normal = n.normalizeThis();
		this.point = p;
	}

	/** Constructs a plane with the normal vector n containing the point (0,0,0). */
	public Plane(Vector n) {
		this.normal = n.normalizeThis();
		this.point = new Point(0,0,0);
	}

	/** Constructs a plane with the normal vector n at distance d from the origin.
	 * If d > 0 then the origin is in the half-space determined by the direction of n  
	 */
	public Plane(Vector n, double d) {
		this.normal = n.normalizeThis();
		this.point = new Point(-d*n.x(), -d*n.y(), -d*n.z());
	}
	
	/** 
	 * Constructs a plane through three points using the first point as defining point. 
	 * The normal of the plane will be decided by the order of p, q and r such that if 
	 * the right hand follows the rotation from p to q to r the thumb points in the 
	 * normals direction. TODO: Test this
	 * 
	 * An error is thrown if the points are collinear.
	 */
	public Plane(Point p, Point q, Point r) {
		if(Point.collinear(p, q, r)) throw new Error("Cant construct plane: Points are collinear");
		normal = p.vectorTo(q).crossThis(p.vectorTo(r)).normalizeThis();
		this.point = p;
	}
	
	/** Constructs a plane bisecting two points */
	public Plane(Point p, Point q) {
		normal = new Vector(p, q).normalizeThis();
		point = Point.getMidpoint(p, q);
	}
	
	private double getD(){
		return -normal.x()*point.x() - normal.y()*point.y() - normal.z()*point.z();
	}
	
	/** Get the point defining this plane. */
	public Point getPoint(){	return point;	}
	
	/** Return the normal defining this plane. */
	public Vector getNormal(){ 	return normal; 	}

	/** Set the normal to n. */
	public void setNormal(Vector n) { 
		this.normal = n; 
	}

	/** Get the projection of p onto this plane. */
	public Point projectPoint(Point p) {
		//Daisy
		/*
		System.out.println("Plane, point = "+point.toString());
		System.out.println("Plane, p = "+p.getCoord(0)+", "+p.getCoord(1)+", "+p.getCoord(2));
		Vector v = point.subtract(p).toVector();
		System.out.println("Plane, v = "+v.toString());
		double dist = v.dot(normal);
		System.out.println("Plane, dist = "+dist);
		Point projPoint = point.subtract(normal.multiply(dist));
		return projPoint;*/
		//Rasmus
		
		double t = normal.x()*p.x() + normal.y()*p.y() + normal.z()*p.z() + getD();
		return new Point(p.x() - normal.x()*t, p.y() - normal.y()*t, p.z() - normal.z()*t);
		
	}

	/** Returns 1/0/-1 if point p is above/on/below this plane. */
	public int above(Point p) {
		double dotP = normal.dot(p.toVector());
		double d = getD();
		if (dotP > -d) return 1;
		if (dotP < -d) return -1; 
		return 0;
	}

	/** Returns 1/0/-1 if point p is below/on/above this plane */
	public int below(Point p) { return -above(p); }

	/** Get the distance of point p to this plane */
	public double getDistance(Point p) { return Math.abs(normal.dot(p.toVector()) + getD()); }

	/** Get the unsigned angle between this plane and p. */
	public double getUnsignedDihedralAngle(Plane p){
		return Math.acos(normal.dot(p.normal));
	}

	/** Get the intersection of a line with the plane. Returns null if the line is 
	 * parallel to plane. */
	public Point getIntersection(Line line) {
		double denom = normal.dot(line.getDir());
		if (denom==0) return null;
		else {
			Point a = line.getP();
			Vector pa = point.vectorTo(a);
			double u = normal.dot(pa)/denom;
			return new Point(a.x() - u*line.dir.x(), a.y() - u*line.dir.y(), a.z() - u*line.dir.z());
		}
	}
	
	/** Get the line-parameter of the intersection between a plane and a line. Returns infinity 
	 * if line is parallel to plane. TODO: Consider moving to Line3d 
	 * */
	public double getIntersectionParameter(Line line) {
		double denom = normal.dot(line.getDir());
		if (denom == 0) return Double.POSITIVE_INFINITY;
		else {
			Point a = line.getP();
			Vector pa = point.vectorTo(a);
			double u = normal.dot(pa)/denom;
			return u;
		}
	}

	/** Get the intersection of a segment with the plane. Returns null if line is 
	 * parallel to plane.*/
	public Point getIntersection(LineSegment sgm) {
		//Daisy
		double dist0 = normal.dot(sgm.getA().subtract(point));
		double dist1 = normal.dot(sgm.getB().subtract(point));
		if (dist0*dist1>0) return null;
		Vector x = (Vector)(sgm.getB().subtract(sgm.getA())).multiplyThis(1/sgm.getB().distance(sgm.getA())).toVector();
		double cos = normal.dot(x);
		if (Math.abs(cos)>=Constants.EPSILON) {
			return sgm.getB().subtract(x.multiply(dist1/cos));
		} else return null;
		//Rasmus
/*		Vector dir = sgm.getAToB();
		double denom = normal.dot(dir);
		if (denom == 0) return null;
		else {
			Vector pa = point.vectorTo(sgm.a);
			double u = normal.dot(pa)/denom;
			if ((u < 0) || (u > 1)) return null;
			else return new Point(sgm.a.x() + u*dir.x(), sgm.a.y() + u*dir.y(), sgm.a.z() + u*dir.z());
		}*/
	}

	public Double getIntersectionAngle(Circle circle, Point p, Vector dir, J3DScene scene) {
		Vector nC = circle.getNormal();
		if (nC.isParallel(normal)) return null;
		Plane circlePlane = new Plane(circle.getCenter(), nC);
		Line line = getIntersection(circlePlane);                       // line.toScene(scene, 0.01, Color.blue);
		double dist = line.getDistance(circle.getCenter());
		if (dist > circle.getRadius() - Constants.EPSILON) return null;
		return circle.getFirstIntersection(line, p, dir, scene);
	}
	
	/** Get intersection of a circle with the plane. Returns null if the plane does not intersect the circle
	 * or if the circle is in the plane. Returns the touch point if the plane is tangent to the circle. 
	 * Otherwise returns 2 points
	 */
	public Point[] getIntersection(Circle circle) {
		Vector u = circle.getNormal().getOrthonormal().multiply(circle.getRadius());
		return getIntersection(circle, u);
	}
	
	public Point[] getIntersection(Circle circle, Vector u) {
		Vector nC = circle.getNormal();
		if (nC.isParallel(normal)) return null;
		Plane circlePlane = new Plane(circle.getCenter(), nC);
		Line line = getIntersection(circlePlane);
		double dist = line.getDistance(circle.getCenter());
		if (dist > circle.getRadius() + Constants.EPSILON) return null;
		if (dist > circle.getRadius() - Constants.EPSILON) {
			Point intPoints[] = new Point[1];
			intPoints[0] = line.orthogonalProjection(circle.getCenter());
			return intPoints;
		}
		Vector v = u.clone();
		nC.rotateIn(v, Math.PI/2);
		Vector cp = new Vector(circle.getCenter(), point);
		double a = u.dot(normal);
		double b = v.dot(normal);
		double c = cp.dot(normal);
		double r = Math.sqrt(a*a + b*b);
		double x = Math.atan2(b/r, a/r);   // atan2(sin x, cos x)
		double alpha1 = Math.acos(c/r); 
		double alpha2 = 2*Math.PI - alpha1;
		double t1 = alpha1 + x;
		double t2 = alpha2 + x;
		Vector[] intVectors = new Vector[2];
		intVectors[0] = u.clone();
		nC.rotateIn(intVectors[0], t1);
		intVectors[1] = u.clone();
		nC.rotateIn(intVectors[1], t2);
		Point intPoints[] = new Point[2];
		intPoints[0] = circle.getCenter().clone().add(intVectors[0]);
		intPoints[1] = circle.getCenter().clone().add(intVectors[1]);
		return intPoints;
		
	}
	
	/** returns the intersection line with another plane*/
	public Line getIntersection(Plane pl) {
		Vector dir = normal.cross(pl.getNormal());
		if (dir.isZeroVector()) return null;
		double h1 = normal.dot(new Vector(point));
		double h2 = pl.getNormal().dot(new Vector(pl.getPoint()));
		double dd = normal.dot(pl.getNormal());
		double denom = 1 - dd*dd;
		double c1 = (h1 - h2*dd)/denom;
		double c2 = (h2 - h1*dd)/denom;
		Point q = new Point(c1*normal.x() + c2*pl.getNormal().x(),
							c1*normal.y() + c2*pl.getNormal().y(),
							c1*normal.z() + c2*pl.getNormal().z());
		return new Line(q, dir);
	}

	
	
	/** Get intersection of a sphere with the plane. Returns null if the plane does not intersect the sphere.
	 * Returns a circle with radius 0 if the plane is tangent to the sphere.
	 * Otherwise returns circle
	 */
	public Circle getIntersection(Sphere sphere) {
		double dist = this.getDistance(sphere.getCenter());
		double rad = sphere.getRadius();
		if (dist - rad > Constants.EPSILON) return null;
		Point center = projectPoint(sphere.getCenter());
		if (dist - rad > -Constants.EPSILON) return new Circle(center, 0, null);
		return new Circle(center, Math.sqrt(rad*rad - center.distanceSquared(sphere.getCenter())), normal); 
	}
	
	/** Returns the defining point for this plane. The center of a plane is not well-defined, so  
	 * to implement the shape interface the defining point is simply used. */
	public Point getCenter() {
		return point.clone();
	}

	public CTetrahedron toScene(J3DScene scene, Color clr, double size) {
		Line line = new Line(point, normal);
		Vector dir = normal.getOrthonormal();
		dir.multiplyThis(size);
		CVertex p1 = new CVertex(point.add(dir), 0); 
		CVertex p2 = new CVertex(line.rotate(p1, Math.PI/2), 0);
		CVertex p3 = new CVertex(line.rotate(p2, Math.PI/2), 0);
		CVertex p4 = new CVertex(line.rotate(p3, Math.PI/2), 0);
		p4.addThis(normal.multiply(0.01));
		CTetrahedron tetr = new CTetrahedron(p1, p2, p3, p4);
		scene.addShape(tetr, clr);	
		return tetr;
	}

	private static void testing31() {
		final Color transp1 = new Color(255,0,0,50);
		final Color transp2 = new Color(0,255,0,50);
		Point a = new Point(0.7, 0.-0.2,  0.3);
		Point b = new Point(0.4, 0.6, 0.8);
		Point c = new Point(0.4, 0.6, 0.2);
		Point d = new Point(-0.6, 0.1, -0.3);
		double alpha = 1;
		Line lineABC = new Line(a, b, c);
		Point projD = lineABC.orthogonalProjection(a);
		double dist = Math.sqrt(alpha*alpha - lineABC.getDistanceSquared(a));
		Vector dir = lineABC.dir.normalize().scaleToLength(dist);
		Point c1 = projD.add(dir);
		Point c2 = projD.subtract(dir);
		Sphere C1 = new Sphere(c1, alpha);
		Sphere C2 = new Sphere(c2, alpha);
		Circle dOrbit = new Circle(new Point(0, 0, d.z()), Math.sqrt(d.x()*d.x() + d.y()*d.y()), new Vector(0, 0, 1));
		Point[] p1 = C1.getIntersections(dOrbit);
		Point[] p2 = C2.getIntersections(dOrbit);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		Sphere sa = new Sphere(a, 0.01);
		scene.addShape(sa, Color.red);
		Sphere sb = new Sphere(b, 0.01);
		scene.addShape(sb, Color.red);
		Sphere sc = new Sphere(c, 0.01);
		scene.addShape(sc, Color.red);
		Sphere sd = new Sphere(d, 0.01);
		scene.addShape(sd, Color.green);
		lineABC.toScene(scene, 0.005, Color.black);
		C1.toScene(scene, transp1);
		C2.toScene(scene, transp2);
		dOrbit.toScene(scene, 0.001, 32);
		for (int i = 0; i < p1.length; i++) scene.addShape(new Sphere(p1[i],0.01), Color.magenta);
		for (int i = 0; i < p2.length; i++) scene.addShape(new Sphere(p2[i],0.01), Color.pink);
	}

	
	private static void testing22() {
		Sphere sphere;
		double alpha = 1.1;
		Point a = new Point(0.6, -0.1,  0.41);
		Point b = new Point(0.4, 0.3, -0.11);
		System.out.println("square radius of the red circle " + (alpha*alpha - a.distanceSquared(b)/4));
		if (a.distanceSquared(b) > alpha*alpha) {
			System.out.println("Points a and b too far apart. No sphere of radius " + alpha + " can contain them both.");
			System.exit(0);
		}

		Point c = new Point(0.4, -0.4, 0.2);  
		Point d = new Point(0.3, -0.2, -0.4);  
		System.out.println("square radius of the green circle " + (alpha*alpha - c.distanceSquared(d)/4));
		if (c.distanceSquared(d) > alpha*alpha) {
			System.out.println("Points c and d too far apart. No sphere of radius " + alpha + " can contain them both.");
			System.exit(0);
		}
		Vector z = new Vector(0 ,0, 1);
		Vector cd = new Vector(c ,d);
		Vector cdBis = new Vector(c, d);
		System.out.println("Angle check: " + Functions.toDeg(z.angle(cd)));
		System.out.println("Angle check bis: " + Functions.toDeg(z.angle(cdBis)));

		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		new Sphere(a, 0.05).toScene(scene, Color.red);
		new Sphere(b, 0.05).toScene(scene, Color.red);
		new LineSegment(a,b).toScene(scene, 0.01, Color.red);
		Plane planeAB = new Plane(a, b);
		Circle circleAB = new Circle(Point.getMidpoint(a, b), Math.sqrt(alpha*alpha - a.distanceSquared(b)/4), planeAB.normal);
		circleAB.toScene(scene,  0.005, 32, Color.red);

		
		new Sphere(c, 0.05).toScene(scene, Color.green);
		new Sphere(d, 0.05).toScene(scene, Color.green);
		Point midPointCD = Point.getMidpoint(c, d);
		new Sphere(midPointCD, 0.05).toScene(scene, Color.green);
		Shape shape = new LineSegment(c,d).toScene(scene, 0.01, Color.green);		
		Plane planeCD = new Plane(c, d);
		Circle circleCD = new Circle(midPointCD, Math.sqrt(alpha*alpha - c.distanceSquared(d)/4), planeCD.normal);

		Circle cOrbit = new Circle(new Point(0, 0, c.z()), Math.sqrt(c.x()*c.x() + c.y()*c.y()), new Vector(0, 0, 1));
		cOrbit.toScene(scene, 0.001, 32);
		Circle dOrbit = new Circle(new Point(0, 0, d.z()), Math.sqrt(d.x()*d.x() + d.y()*d.y()), new Vector(0, 0, 1));
		dOrbit.toScene(scene, 0.001, 32);
		Circle cdOrbit = new Circle(new Point(0, 0, (c.z()+d.z())/2), Math.sqrt(midPointCD.x()*midPointCD.x() + midPointCD.y()*midPointCD.y()), new Vector(0, 0, 1));
		cdOrbit.toScene(scene, 0.001, 32);		
		
		Line zaxis = new Line(new Point(0,0,0), new Vector(0,0,1));
		zaxis.toScene(scene, 0.005, Color.black);

		Cylinder cyl[] = circleCD.toScene(scene,  0.005, 32, Color.green);
		

		int steps = 200;
		double angle = 0.0;
		double delta = (2*Math.PI)/steps;

		for (int i = 0; i < steps; i++) {
			scene.removeShape(shape);
			c.rotationCW(zaxis.dir, delta);
			d.rotationCW(zaxis.dir, delta);
			cd = new Vector(c ,d);
			cdBis = new Vector(cdBis.x()*Math.cos(delta) - cdBis.y()*Math.sin(delta), cdBis.y()*Math.cos(delta) + cdBis.x()*Math.sin(delta), cdBis.z());
			System.out.println("Angle check: " + Functions.toDeg(z.angle(cd)));
			System.out.println("Angle check bis: " + Functions.toDeg(z.angle(cdBis)));
			midPointCD.rotationCW(zaxis.dir, delta);
			scene.repaint();
			new LineSegment(c,d);
			shape = new LineSegment(c, d).toScene(scene, 0.01, Color.green);
			
			planeCD = new Plane(c, d);
//			lineABCD = planeAB.getIntersection(planeCD);
//			lineABCD.toScene(scene, 0.01, Color.blue);
			
			circleCD = new Circle(midPointCD, Math.sqrt(alpha*alpha - c.distanceSquared(d)/4), planeCD.normal);
//			circleCD.fromScene(scene, cyl);
			cyl = circleCD.toScene(scene,  0.001, 32, new Color(i%256,100, 0));
			Point[] intersections = planeAB.getIntersection(circleCD);
			if (intersections != null) {
				new Sphere(intersections[0], 0.01).toScene(scene, Color.blue);
				System.out.println(intersections[0].distance(circleAB.getCenter()) - alpha);
				if (intersections.length == 2) {
					new Sphere(intersections[1], 0.01).toScene(scene, Color.blue);
					System.out.println("   " + (intersections[1].distance(circleAB.getCenter()) - alpha));
				}
			}
			sphere = new Sphere(a, b, c, d);
			angle += delta;
			new Point(angle, sphere.getRadius(), 0.0).toScene(scene, 0.02, Color.black, 32);
			scene.repaint();
		}
	}

	private static Vector m(Vector v, double beta) {
		return new Vector(v.x()*Math.cos(beta) - v.y()*Math.sin(beta), v.y()*Math.cos(beta) + v.x()*Math.sin(beta), v.z());
	}
	private static Point q(Point q0, double beta) {
		return new Point(q0.distance()*Math.cos(beta), q0.distance()*Math.sin(beta), q0.z());
	}

	private static void  testing() {
		double r1 = 2;
		double r2 = 3;
		double beta1 = Math.PI/10;
		double beta2 = Math.PI/13;
		double beta = Math.PI/9;
		double m1x = 0.6;
		double m2x =-3;
		double m1y = 7;
		double m2y = 12;
		double m1z = 6;
		double m2z = 3;
		double x = 12;
		double y = 6;
		double z = 2;
		double alpha = 1;
		
		double value = 
		 (r2*Math.cos(beta2+beta) - r1*Math.cos(beta1+beta))*(x - 0.5*r1*Math.cos(beta1+beta) -0.5*r2*Math.cos(beta2+beta)) +
		 (r2*Math.sin(beta2+beta) - r1*Math.sin(beta1+beta))*(y - 0.5*r1*Math.sin(beta1+beta) -0.5*r2*Math.sin(beta2+beta)) +
		 (m2z-m1z)*(z - 0.5*(m1z+m2z));
				
		double value1 =
		 x*r2*Math.cos(beta2+beta) - x*r1*Math.cos(beta1+beta) -
		 0.5*r1*r2*Math.cos(beta2+beta)*Math.cos(beta1+beta) - 0.5*r2*r2*Math.cos(beta2+beta)*Math.cos(beta2+beta) +
		 0.5*r1*r1*Math.cos(beta1+beta)*Math.cos(beta1+beta) + 0.5*r1*r2*Math.cos(beta1+beta)*Math.cos(beta2+beta) +
		 y*r2*Math.sin(beta2+beta) -y*r1*Math.sin(beta1+beta) -
		 0.5*r1*r2*Math.sin(beta2+beta)*Math.sin(beta1+beta) - 0.5*r2*r2*Math.sin(beta2+beta)*Math.sin(beta2+beta) +
		 0.5*r1*r1*Math.sin(beta1+beta)*Math.sin(beta1+beta) + 0.5*r1*r2*Math.sin(beta1+beta)*Math.sin(beta2+beta) +
		 z*(m2z - m1z) - 0.5*(m2z-m1z)*(m1z+m2z);
		double value2 =
		 x*r2*Math.cos(beta2+beta) - x*r1*Math.cos(beta1+beta) +
		 y*r2*Math.sin(beta2+beta) -y*r1*Math.sin(beta1+beta) +
		 z*(m2z - m1z) - 0.5*(m2z-m1z)*(m1z+m2z) + 0.5*(r1*r1-r2*r2);

		double mvalue =
		 (x - 0.5*r1*Math.cos(beta1+beta) + 0.5*r2*Math.cos(beta2+beta))*(x - 0.5*r1*Math.cos(beta1+beta) + 0.5*r2*Math.cos(beta2+beta)) +
		 (y - 0.5*r1*Math.sin(beta1+beta) + 0.5*r2*Math.sin(beta2+beta))*(y - 0.5*r1*Math.sin(beta1+beta) + 0.5*r2*Math.sin(beta2+beta)) +
		 (z - 0.5*(m1z+m2z))*(z - 0.5*(m1z+m2z)) - alpha*alpha + 0.25*((m2x-m1x)*(m2x-m1x)+(m2y-m1y)*(m2y-m1y)+(m2z-m1z)*(m2z-m1z));

		double mvalue1 =
		 x*x - x*r1*Math.cos(beta1+beta) + x*r2*Math.cos(beta2+beta) -
		 0.5*r1*r2*Math.cos(beta1+beta)*Math.cos(beta2+beta) + 0.25*r1*r1*Math.cos(beta1+beta)*Math.cos(beta1+beta) + 0.25*r2*r2*Math.cos(beta2+beta)*Math.cos(beta2+beta) +
		 y*y - y*r1*Math.sin(beta1+beta) + y*r2*Math.sin(beta2+beta) -
		 0.5*r1*r2*Math.sin(beta1+beta)*Math.sin(beta2+beta) + 0.25*r1*r1*Math.sin(beta1+beta)*Math.sin(beta1+beta) + 0.25*r2*r2*Math.sin(beta2+beta)*Math.sin(beta2+beta) +
		 z*z - z*(m1z+m2z) + 0.25*(m1z+m2z)*(m1z+m2z) - alpha*alpha + 0.25*((m2x-m1x)*(m2x-m1x)+(m2y-m1y)*(m2y-m1y)+(m2z-m1z)*(m2z-m1z));
	
		double mvalue2 =
		 x*x - x*r1*Math.cos(beta1+beta) + x*r2*Math.cos(beta2+beta) -
		 0.5*r1*r2*Math.cos(beta1+beta)*Math.cos(beta2+beta) +
		 y*y - y*r1*Math.sin(beta1+beta) + y*r2*Math.sin(beta2+beta) -
		 0.5*r1*r2*Math.sin(beta1+beta)*Math.sin(beta2+beta) +
		 z*z - z*(m1z+m2z) + 0.25*(m1z+m2z)*(m1z+m2z) + 0.25*(r1*r1+r2*r2) - alpha*alpha + 0.25*((m2x-m1x)*(m2x-m1x)+(m2y-m1y)*(m2y-m1y)+(m2z-m1z)*(m2z-m1z));
	
		double mvalue3 =
		 x*x - x*r1*Math.cos(beta1+beta) + x*r2*Math.cos(beta2+beta) +
		 y*y - y*r1*Math.sin(beta1+beta) + y*r2*Math.sin(beta2+beta) +
		 z*z - z*(m1z+m2z) +
		 0.25*(m1z+m2z)*(m1z+m2z) + 0.25*(r1*r1 + r2*r2) - 0.5*r1*r2*(Math.cos(beta1)*Math.cos(beta2) + Math.sin(beta1)*Math.sin(beta2)) - alpha*alpha + 0.25*((m2x-m1x)*(m2x-m1x)+(m2y-m1y)*(m2y-m1y)+(m2z-m1z)*(m2z-m1z));
				
		double difference =
		 x*x + y*y + z*z - 2*z*m2z +
		 0.5*(m2z*m2z - m1z*m1z) + 0.25*(m1z+m2z)*(m1z+m2z) + 0.75*r2*r2 - 0.25*r1*r1  - 
		 0.5*r1*r2*(Math.cos(beta1)*Math.cos(beta2) + Math.sin(beta1)*Math.sin(beta2)) - alpha*alpha + 0.25*((m2x-m1x)*(m2x-m1x)+(m2y-m1y)*(m2y-m1y)+(m2z-m1z)*(m2z-m1z));

		double differencem = mvalue3 - value2;
		
		alpha = 1;
		
		Line zaxis = new Line(new Point(0,0,0), new Vector(0,0,1));
		Line xaxis = new Line(new Point(0,0,0), new Vector(1,0,0));
		Line yaxis = new Line(new Point(0,0,0), new Vector(0,1,0));

//		Point a = new Point(0.6, -0.1,  0.21);
		Point a = new Point(0.3, -0.51,  0.341);
		Point b = new Point(0.4, 0.1, -0.11);
		System.out.println("square radius of the red circle " + (alpha*alpha - a.distanceSquared(b)/4));
		if (a.distanceSquared(b) > alpha*alpha) {
			System.out.println("Points a and b too far apart. No sphere of radius " + alpha + " can contain them both.");
			System.exit(0);
		}

		Point c = new Point(0.4, 0.0, 0.2);  
		Point d = new Point(0.1, 0.0, -0.2);  

//		Point c = new Point(0.4, 0.2, 0.2);  
//		Point d = new Point(-0.1, -0.2, -0.2);  
//		Point d = new Point(0.4, 0.2, -0.2);
		System.out.println("square radius of the green circle " + (alpha*alpha - c.distanceSquared(d)/4));
		if (c.distanceSquared(d) > alpha*alpha) {
			System.out.println("Points c and d too far apart. No sphere of radius " + alpha + " can contain them both.");
			System.exit(0);
		}
		Point midPointCD = Point.getMidpoint(c, d);

		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		zaxis.toScene(scene, 0.002,	Color.red);
		xaxis.toScene(scene, 0.002,	Color.black);
		yaxis.toScene(scene, 0.002,	Color.black);
		new Sphere(a, 0.05).toScene(scene, Color.red);
		new Sphere(b, 0.05).toScene(scene, Color.red);
		Plane planeAB = new Plane(a, b);
		Circle circleAB = new Circle(Point.getMidpoint(a, b), Math.sqrt(alpha*alpha - a.distanceSquared(b)/4), planeAB.normal);
		circleAB.toScene(scene,  0.005, 32, Color.red);

		Vector n  = new Vector(a, b);
		Point p = Point.getMidpoint(a,  b);
//		n.toScene(scene, p, Color.red, 0.005);
		
		Vector m0 = new Vector(c, d);
		Point q0 = Point.getMidpoint(c, d);
		
		
		int steps = 20;
		double delta = (2*Math.PI)/steps;

		if (Math.abs(n.z()) < Constants.EPSILON) {
			System.out.println("n.z too small ");
			System.exit(0);
		}

		double Ax = b.x() - a.x();
		double Ay = b.y() - a.y();
		double Az = b.z() - a.z();
		double A = (a.distanceSquared() - b.distanceSquared())/2;
		
		double Bx = -(a.x() + b.x());
		double By = -(a.y() + b.y());
		double Bz = -(a.z() + b.z());
		double B = (Bx*Bx + By*By + Bz*Bz)/4 -alpha*alpha + a.distanceSquared(b)/4;
		
		double Cz = -2*d.z();
		double cosb1 = c.x()/Math.sqrt(c.x()*c.x() + c.y()*c.y());
		double cosb2 = d.x()/Math.sqrt(d.x()*d.x() + d.y()*d.y());
		double sinb1 = c.y()/Math.sqrt(c.x()*c.x() + c.y()*c.y());
		double sinb2 = d.y()/Math.sqrt(d.x()*d.x() + d.y()*d.y());

		double C = (d.z()*d.z() - c.z()*c.z())/2 +
				   (c.z()+ d.z())*(c.z()+ d.z())/4 +
				   0.75*(d.x()*d.x() + d.y()*d.y()) -
				   0.25*(c.x()*c.x() + c.y()*c.y()) -
				   0.5*Math.sqrt((c.x()*c.x() + c.y()*c.y())*(d.x()*d.x() + d.y()*d.y()))*(cosb1*cosb2 + sinb1*sinb2) -
				   alpha*alpha + Point.getMidpoint(c, d).distanceSquared()/4;
		
		double ratio = m0.z()/n.z();
		Point onLine;
		for (int i = 0; i < steps; i++) {
			new Sphere(c, 0.03).toScene(scene, Color.green);
			new Sphere(d, 0.03).toScene(scene, Color.green);
			new Sphere(midPointCD, 0.01).toScene(scene, Color.green);
			Plane planeCD = new Plane(c, d);
			Circle circleCD = new Circle(midPointCD, Math.sqrt(alpha*alpha - c.distanceSquared(d)/4), planeCD.normal);
			circleCD.toScene(scene, 0.001, 32, Color.green);

			Point[] blackPoints = planeAB.getIntersection(circleCD);
			if (blackPoints != null)
				for (int k = 0; k < blackPoints.length; k++) {
					new Sphere(blackPoints[k], 0.01).toScene(scene, Color.black);
			}
			Point[] bluePoints = planeCD.getIntersection(circleAB);
			if (bluePoints != null)
				for (int k = 0; k < bluePoints.length; k++) {
					new Sphere(bluePoints[k], 0.01).toScene(scene, Color.blue);
			}
			
			Vector mBeta = m(m0, i*delta);   // direction vector between the rotating points C an D
			Point qBeta = q(q0, i*delta);    // midpoint of rotating points C and D
//			mBeta.toScene(scene, qBeta, Color.green, 0.005);

			Vector nxmBeta = n.cross(mBeta).normalize(); // direction vector of the intersection line of planes AB (fixed) and CD (rotating) 
			
			if (Math.abs(n.y()) > Constants.EPSILON) {
				z = (n.y()*qBeta.dot(mBeta) - mBeta.y()*p.dot(n))/(n.y()*mBeta.z() - mBeta.y()*n.z());
				onLine = new Point(0, (p.dot(n) - n.z()*z)/n.y(), z); // point on the intersection line
			}
			else {
				onLine = new Point(0, 0, p.dot(n)/n.z());
			}
			Line l = new Line(onLine, nxmBeta);     // intersection line
			l.toScene(scene, 0.0005, Color.magenta);
//			new Sphere(onLine, 0.01).toScene(scene, Color.magenta);
			mBeta.toScene(scene, qBeta, Color.green, 0.001);
			n.toScene(scene, qBeta, Color.red, 0.001);
			nxmBeta.toScene(scene, qBeta, Color.magenta, 0.001);
			c.rotationCW(zaxis.dir, delta);
			d.rotationCW(zaxis.dir, delta);
			midPointCD.rotationCW(zaxis.dir, delta);

		}
		

	}
	
	private static void testing21(Point a, Point b, Point c, Point d) {
		final Color transp1 = new Color(255,0,0,50);
		final Color transp2 = new Color(0,255,0,50);

		double alpha = 0.7;
		Plane planeAB = new Plane(a, b);
		Circle circleAB = new Circle(Point.getMidpoint(a, b), Math.sqrt(alpha*alpha - a.distanceSquared(b)/4), planeAB.normal);
		Line axis = new Line(new Point(0,0,0), new Vector(0,0,1));
		J3DScene scene = J3DScene.createJ3DSceneInFrame();

		Sphere sa = new Sphere(a, 0.01);
		scene.addShape(sa, Color.black);
		Sphere sb = new Sphere(b, 0.01);
		scene.addShape(sb, Color.blue);
		Sphere sc = new Sphere(c, 0.01);
		Sphere sd = new Sphere(d, 0.01);
		planeAB.toScene(scene, transp1, 1);
		circleAB.toScene(scene,  0.005, 32);

		Circle cOrbit = new Circle(new Point(0, 0, c.z()), Math.sqrt(c.x()*c.x() + c.y()*c.y()), new Vector(0, 0, 1));
		cOrbit.toScene(scene, 0.001, 32);
		Vector cOrbitNormal = cOrbit.getNormal();
		Point cOrbitCenter = cOrbit.getCenter();
		Circle dOrbit = new Circle(new Point(0, 0, d.z()), Math.sqrt(d.x()*d.x() + d.y()*d.y()), new Vector(0, 0, 1));
		dOrbit.toScene(scene, 0.001, 32);
		Vector dOrbitNormal = dOrbit.getNormal();
		Point dOrbitCenter = dOrbit.getCenter();

		Sphere sphereC = new Sphere(c, alpha);
		Circle intCircleC = null;
		Cylinder[] sIntCircleC = null;
		Sphere sphereD = new Sphere(d, alpha);
		Circle intCircleD = null;
		Cylinder[] sIntCircleD = null;
		int steps = 100;
		for (int i = 0; i < steps; i++) {
			intCircleC = planeAB.getIntersection(sphereC);
			scene.addShape(sc, Color.red);			
			sphereC.toScene(scene, transp2);
			if (intCircleC != null) sIntCircleC = intCircleC.toScene(scene, 0.005, 32, Color.red);
			intCircleD = planeAB.getIntersection(sphereD);
			scene.addShape(sd, Color.green);			
			sphereD.toScene(scene, transp1);
			if (intCircleD != null) sIntCircleD = intCircleD.toScene(scene, 0.005, 32, Color.green);
			scene.removeShape(sphereC); scene.removeShape(sphereD);
			if (intCircleC != null) intCircleC.fromScene(scene, sIntCircleC);
			scene.removeShape(sc);
			c.rotation(cOrbitNormal, (2*Math.PI)/steps, cOrbitCenter);
			if (intCircleD != null) intCircleD.fromScene(scene, sIntCircleD);
			scene.removeShape(sd);
			c.rotation(cOrbitNormal, (2*Math.PI)/steps, cOrbitCenter);
			d.rotation(dOrbitNormal, (2*Math.PI)/steps, dOrbitCenter);
		}
	}
	
	private static void test() {
		Point S1 = new Point(0.4, 0.2, 0.3);
		Point S2 = new Point(-0.2, 0.6, -0.3);
		Point R1 = new Point(0.4, 0.6, -0.4);
		Point R2 = new Point(0.8, 0.4, -0.5);
		
		Vector n1 = new Vector(0, 0, 1);
		Vector u1 = new Vector(R1.x(), R1.y(), 0).normalize();
		Vector v1 = n1.cross(u1);
		System.out.println("length: " + v1.length());
		Circle C1 = new Circle(new Point(0, 0, R1.z()), Math.sqrt(R1.x()*R1.x() + R1.y()*R1.y()), n1);

		Vector n2 = new Vector(0, 0, 1);
		Vector u2 = new Vector(R2.x(), R2.y(), 0).normalize();
		Vector v2 = n2.cross(u2);
		System.out.println("length: " + v2.length());
		Circle C2 = new Circle(new Point(0, 0, R2.z()), Math.sqrt(R2.x()*R2.x() + R2.y()*R2.y()), n2);

		double alpha = 3;
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		
		Plane planeS12 = new Plane(S1, S2);
		Circle circleS12 = new Circle(Point.getMidpoint(S1, S2), Math.sqrt(alpha*alpha - S1.distanceSquared(S2)/4), planeS12.normal);
		circleS12.toScene(scene,  0.005, 32, Color.red);

		S1.toScene(scene, 0.05, Color.red);
		S2.toScene(scene, 0.05, Color.red);
		Sphere S;
		
		int steps = 100;
		double delta = (2*Math.PI)/steps;
		double theta = 0.0;
		Point Q1 = R1;
		Point Q2 = R2;
		Sphere sph1 = null;
		Sphere sph2 = null;

		for (int i = 0; i < steps; i++) {
			Q1 = C1.getCenter().add(u1.multiply(C1.getRadius()*Math.cos(theta)).add(v1.multiply(C1.getRadius()*Math.sin(theta))));
			Q2 = C2.getCenter().add(u2.multiply(C2.getRadius()*Math.cos(theta)).add(v2.multiply(C2.getRadius()*Math.sin(theta))));
			scene.removeShape(sph1);
			scene.removeShape(sph2);
			sph1 = Q1.toScene(scene, 0.05, Color.blue);
			sph2 = Q2.toScene(scene, 0.05, Color.blue);
			S = new Sphere(S1, S2, Q1, Q2); 
			S.toScene(scene, new Color(255,0,0,2));
			S.getCenter().toScene(scene, 0.02, Color.black);
			System.out.println(S.getRadius());
			theta += delta;
		}
	}

	
	private static void testOnceAgain(Point A, Point B, Point C, Point D, Line L) {
		/* assumption 1: R1 and R2 rotate around the z-axis. If not, translation and rotation of points can achieve this 
		 * assumption 2: R1 is in the xy-plane. If not, translation along z-coordinate can be applied. It does not affect the rotation axis
		 * assumption 3: radius of R1 is 1. If not, appropriate scaling can be carried out. It does not affect previous assumptions.
		 * assumption 4: The initial position of R1 is (1, 0, 0). Otherwise appropriate rotation of all 4 points around the z-axis can be carried out. 
		 * 				 It does not affect the previous assumptions.  
		 */
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();

		Point S1 = new Point(A);
		Point S2 = new Point(B);	
		Point R1 = new Point(C);
		Point R2 = new Point(D);
		Line rotAxis = L.clone();
		Cylinder cylRotAxis = rotAxis.toScene(scene, 0.001, Color.red);		
		Sphere origo = rotAxis.p.toScene(scene, 0.01, Color.pink);
		Point projR1 = rotAxis.orthogonalProjection(R1);
		Circle cR1 = new Circle(projR1, R1.distance(projR1), rotAxis.dir);
		Cylinder[] cylCR1 = cR1.toScene(scene, 0.004, 32, Color.green);
		Point projR2 = rotAxis.orthogonalProjection(R2);
		Circle cR2 = new Circle(projR2, R2.distance(projR2), rotAxis.dir);
		Cylinder[] cylCR2 = cR2.toScene(scene, 0.004, 32, Color.red);
		
		Shape SS1 = S1.toScene(scene, 0.01, Color.black);
		Shape SS2 = S2.toScene(scene, 0.01, Color.blue);
		Shape SR1 = R1.toScene(scene, 0.01, Color.green);
		Shape SR2 = R2.toScene(scene, 0.01, Color.red);
		
		Vector exVector = new Vector(1,0,0);
		Vector eyVector = new Vector(0,1,0);
		Vector ezVector = new Vector(0,0,1);
		LSS xAxisShape = exVector.toScene(scene, Color.black, 0.002);
		LSS yAxisShape = eyVector.toScene(scene, Color.black, 0.002);
		LSS zAxisShape = ezVector.toScene(scene, Color.black, 0.002);
		TextShape xText = new TextShape("x", new Point(1,0,0), 0.05);
		TextShape yText = new TextShape("y", new Point(0,1,0), 0.05);
		TextShape zText = new TextShape("z", new Point(0,0,1), 0.05);
		scene.addShape(xText);
		scene.addShape(yText);
		scene.addShape(zText);
		
		// Translates the points so that the rotation axis goes through origo
		S1.subtractThis(rotAxis.p); scene.removeShape(SS1); 
			SS1 = S1.toScene(scene, 0.01, Color.black);
		S2.subtractThis(rotAxis.p); scene.removeShape(SS2); 
			SS2 = S2.toScene(scene, 0.01, Color.blue);

		R1.subtractThis(rotAxis.p); scene.removeShape(SR1); SR1 = R1.toScene(scene, 0.01, Color.green);
		cR1.fromScene(scene, cylCR1);
		cR1 = new Circle(cR1.getCenter().subtract(rotAxis.p), cR1.getRadius(), cR1.getNormal());
		cylCR1 = cR1.toScene(scene, 0.004, 32, Color.green);

		R2.subtractThis(rotAxis.p); scene.removeShape(SR2); SR2 = R2.toScene(scene, 0.01, Color.red);
		cR2.fromScene(scene, cylCR2);
		cR2 = new Circle(cR2.getCenter().subtract(rotAxis.p), cR2.getRadius(), cR2.getNormal());
		cylCR2 = cR2.toScene(scene, 0.004, 32, Color.red);

		rotAxis.p = new Point(0,0,0);
		scene.removeShape(cylRotAxis);
		cylRotAxis = rotAxis.toScene(scene, 0.001, Color.red);	
		scene.removeShape(origo);
		origo = rotAxis.p.toScene(scene, 0.01, Color.pink);

		
		// changes coordinate system so that the rotation axis becomes identical with the z-axis
		Vector rotV = rotAxis.dir.cross(new Vector(0,0,1));
		LSS lss = rotV.toScene(scene, Color.magenta, 0.01);
		double angle = rotAxis.dir.angle(new Vector(0,0,1));
		S1.rotationCW(rotV, angle);
		S2.rotationCW(rotV, angle);
		R1.rotationCW(rotV, angle);
		R2.rotationCW(rotV, angle);
		scene.removeShape(SS1); SS1 = S1.toScene(scene, 0.01, Color.black);
		scene.removeShape(SS2); SS2 = S2.toScene(scene, 0.01, Color.blue);
		scene.removeShape(SR1); SR1 = R1.toScene(scene, 0.01, Color.green);
		scene.removeShape(SR2); SR2 = R2.toScene(scene, 0.01, Color.red);
		scene.removeShape(cylRotAxis);
		scene.removeShape(origo);


		rotAxis.dir = new Vector(0,0,1);		
		cR1.fromScene(scene, cylCR1);
		projR1 = rotAxis.orthogonalProjection(R1);
		cR1 = new Circle(projR1, R1.distance(projR1), rotAxis.dir);
		cylCR1 = cR1.toScene(scene, 0.004, 32, Color.green);

		cR2.fromScene(scene, cylCR2);
		projR2 = rotAxis.orthogonalProjection(R2);
		cR2 = new Circle(projR2, R2.distance(projR2), rotAxis.dir);
		cylCR2 = cR2.toScene(scene, 0.004, 32, Color.red);
		scene.removeShape(lss);

		// translates the points so that R1 is on the xy-plane
		Vector transl = new Vector(0,0, R1.z());
		S1.subtractThis(transl); scene.removeShape(SS1); SS1 = S1.toScene(scene, 0.01, Color.black);
		S2.subtractThis(transl); scene.removeShape(SS2); SS2 = S2.toScene(scene, 0.01, Color.blue);
		R1.subtractThis(transl); scene.removeShape(SR1); SR1 = R1.toScene(scene, 0.01, Color.green);
		R2.subtractThis(transl); scene.removeShape(SR2); SR2 = R2.toScene(scene, 0.01, Color.red);

		cR1.fromScene(scene, cylCR1);
		cR1 = new Circle(cR1.getCenter().subtract(transl), cR1.getRadius(), cR1.getNormal());
		cylCR1 = cR1.toScene(scene, 0.004, 32, Color.green);
		cR2.fromScene(scene, cylCR2);
		cR2 = new Circle(cR2.getCenter().subtract(transl), cR2.getRadius(), cR2.getNormal());
		cylCR2 = cR2.toScene(scene, 0.004, 32, Color.red);

		// scales so the radius of R1 is 1
		double factor = 1/R1.distance();
		S1.scaleThis(factor); 
		scene.removeShape(SS1); SS1 = S1.toScene(scene, 0.01, Color.black);
		S2.scaleThis(factor); scene.removeShape(SS2); SS2 = S2.toScene(scene, 0.01, Color.blue);
		R1.scaleThis(factor); scene.removeShape(SR1); SR1 = R1.toScene(scene, 0.01, Color.green);
		R2.scaleThis(factor); scene.removeShape(SR2); SR2 = R2.toScene(scene, 0.01, Color.red);

		cR1.fromScene(scene, cylCR1);
		cR1 = new Circle(cR1.getCenter(), 1, cR1.getNormal());
		cylCR1 = cR1.toScene(scene, 0.004, 32, Color.green);
		cR2.fromScene(scene, cylCR2);
		cR2 = new Circle(new Point(0, 0, R2.z()), cR2.getRadius()*factor, cR2.getNormal());
		cylCR2 = cR2.toScene(scene, 0.004, 32, Color.red);
		
		
		// rotates so that R1 is on the positive x-axis
		angle = new Vector(1,0,0).angle(R1.toVector());
		S1.rotationCW(rotAxis.dir, angle);
		S2.rotationCW(rotAxis.dir, angle);
		R1.rotationCW(rotAxis.dir, angle);
		R2.rotationCW(rotAxis.dir, angle);
		scene.removeShape(SS1); SS1 = S1.toScene(scene, 0.01, Color.black);
		scene.removeShape(SS2); SS2 = S2.toScene(scene, 0.01, Color.blue);
		scene.removeShape(SR1); SR1 = R1.toScene(scene, 0.01, Color.green);
		scene.removeShape(SR2); SR2 = R2.toScene(scene, 0.01, Color.red);

		
		Vector n1 = new Vector(0,0,1);
		
		Vector n2 = new Vector(0, 0, 1);
		Vector u2 = new Vector(R2.x(), R2.y(), 0).normalize();
		Vector v2 = n2.cross(u2);
		double r2 = Math.sqrt(R2.x()*R2.x() + R2.y()*R2.y());
		Circle C2 = new Circle(new Point(0, 0, R2.z()), r2, n2);

		double alpha = 1.5;
		
		Plane planeS12 = new Plane(S1, S2);
		Circle circleS12 = new Circle(Point.getMidpoint(S1, S2), Math.sqrt(alpha*alpha - S1.distanceSquared(S2)/4), planeS12.normal);
		circleS12.toScene(scene,  0.005, 32, Color.red);
		 
		double a = -planeS12.normal.x()/planeS12.normal.z();
		double b = -planeS12.normal.y()/planeS12.normal.z();
		double c = planeS12.normal.dot(Point.getMidpoint(S1, S2))/planeS12.normal.z();

		S1.toScene(scene, 0.05, Color.black);
		S2.toScene(scene, 0.05, Color.blue);
		Sphere S = null;
		
		int steps = 100;
		double delta = (2*Math.PI)/steps;
		double theta = 0.0;
		Point R1t = R1;
		Point R2t = R2;
		Vector nt;
		Point Mt;
		Plane planeR12t;
		Line lt;
		double at, bt, ct;
		Sphere sph1 = null;
		Sphere sph2 = null;
		Sphere oldS = null;
		Point A1 = null;
		Point A2 = null;
		Point A3 = null; 
		Point A4 = null;
		LineSegment sgmR12;
		boolean increasing = false;
		boolean first = true;
		
		for (int i = 0; i < steps; i++) {
			R1t = new Point(Math.cos(theta), Math.sin(theta), 0);
			R2t = new Point(r2*u2.x()*Math.cos(theta) + r2*v2.x()*Math.sin(theta), r2*u2.y()*Math.cos(theta) + r2*v2.y()*Math.sin(theta), R2.z());
			nt = new Vector(R1t, R2t).normalize();
			Mt = Point.getMidpoint(R1t, R2t);
			planeR12t = new Plane(R2t, R1t);
			Circle circleR12 = new Circle(Point.getMidpoint(R1t, R2t), Math.sqrt(alpha*alpha - R1t.distanceSquared(R2t)/4), planeR12t.normal);
			circleR12.toScene(scene,  0.001, 32, Color.green);
			sgmR12 = new LineSegment(R1t, R2t);
			sgmR12.toScene(scene, 0.005, Color.pink);
			
			// verification of normal vector
			Vector ntCheck = new Vector((r2*u2.x()-1)*Math.cos(theta), (r2*v2.y()-1)*Math.sin(theta), R2.z()).normalize();
			System.out.println("Verification: nt = " + nt + " ntCheck = " + ntCheck); 
			
			lt = planeR12t.getIntersection(planeS12);
			
			// verification of the intersection line between two planes
			at = R2.z()*a + (r2*u2.x()-1)*Math.cos(theta);
			bt = R2.z()*b + (r2*v2.y()-1)*Math.sin(theta);
			ct = r2*r2*u2.x()*u2.x() +(R2.z()*R2.z() -1)/2;
			
			scene.removeShape(sph1); sph1 = R1t.toScene(scene, 0.05, Color.green);
			scene.removeShape(sph2); sph2 = R2t.toScene(scene, 0.05, Color.red);
			oldS = S;
			S = new Sphere(S1, S2, R1t, R2t); 
//			S.toScene(scene, new Color(255,0,0,15));
			S.getCenter().toScene(scene, 0.02, Color.black);
			System.out.println(S.getRadius());
			if ((oldS != null) && (S.getRadius() > oldS.getRadius())) increasing = true;
			else {
				if (oldS != null) {
					if (increasing) {
						if (first) {
							A1 = S.getCenter();
							A2 = oldS.getCenter();
							first = false;
							increasing = false;
						}
						else {
							A3 = S.getCenter();
							A4 = oldS.getCenter();
							increasing = false;
						}
					}
				}
			}
			theta += delta;
		}
		Line asympt1 = new Line(A1, A2);
		asympt1.toScene(scene, 0.01, Color.yellow);
		Line asympt2 = new Line(A3, A4);
		asympt2.toScene(scene, 0.01, Color.yellow);

	}
	
	
	private static void fourSpheres(Point S1, Point S2, Point R1, Point R2) {
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		double alpha = 1;
		Sphere sphereS1 = new Sphere(S1, alpha);
		scene.addShape(sphereS1, new Color(255, 0, 0, 50));
		Sphere sphereS2 = new Sphere(S2, alpha);
		scene.addShape(sphereS2, new Color(255, 0, 0, 50));
		Sphere sphereR1;
		Sphere sphereR2;
		int steps = 20;
		double delta = (2*Math.PI)/steps;
		double theta = 0.0;
		Vector n2 = new Vector(0, 0, 1);
		Vector u2 = new Vector(R2.x(), R2.y(), 0).normalize();
		Vector v2 = n2.cross(u2);
		double r2 = Math.sqrt(R2.x()*R2.x() + R2.y()*R2.y());
		Circle C2 = new Circle(new Point(0, 0, R2.z()), r2, n2);
		Point R1t = R1;
		Point R2t = R2;
		for (int i = 0; i < steps; i++) {
			R1t = new Point(Math.cos(theta), Math.sin(theta), 0);
			R2t = new Point(r2*u2.x()*Math.cos(theta) + r2*v2.x()*Math.sin(theta), r2*u2.y()*Math.cos(theta) + r2*v2.y()*Math.sin(theta), R2.z());
			sphereR1 = new Sphere(R1t, alpha);
			scene.addShape(sphereR1, new Color(0,0,255,20));
			sphereR2 = new Sphere(R2t, alpha);
			scene.addShape(sphereR2, new Color(0,0,255,20));
			theta += delta;
		}
		
	}
	public static void main(String[] args) {		
		Point S1 = new Point(0.4, 0.4, -0.1);
//		Point S2 = new Point(-0.2, 1, 0.4);       // hyperbola
//		Point S2 = new Point(0.2, 0.1, 0.4);      // ellipse
		Point S2 = new Point(0.4, 0.4, 0.4);       // hyperbola 2 intersection
		Point R1 = new Point(1.0, 0.0, 0.4);
		double angle = 0.2*Math.PI;
		double r2 = 1;
		Point R2 = new Point(r2*Math.cos(angle), r2*Math.sin(angle), 1);
		Line L = new Line(new Point(0,0,0), new Vector(0, 0, 1));
		Plane.testing22();
//		Plane.testing21(S1, S2, R1, R2);
//		Plane.testOnceAgain(S1, S2, R1, R2, L);
//		Plane.fourSpheres(S1,  S2, R1, R2);

	}

}
