package ProGAL.geom2d;

import java.awt.Color;
import java.math.BigDecimal;
import java.util.List;

import ProGAL.geom2d.Triangulation.TriangulationAlgorithm;
import ProGAL.geom2d.delaunay.Vertex;
import ProGAL.geom2d.viewer.J2DScene;
//import ProGAL.geom3d.Point;
//import ProGAL.geom3d.kineticDelaunay.Vertex;
import ProGAL.math.Constants;
import ProGAL.math.Functions;
import ProGAL.math.Trigonometry;

public class Circle implements Shape{
	protected Point center;
	protected double radius;

	/** Construct a circle with the given center and radius. */
	public Circle(Point center, double radius) {
		this.center = center;
		this.radius = radius;
	}
 
	/** Constructs a circle that is a copy of a given circle */
	public Circle(Circle c) {
		center = new Point(c.center.x(), c.center.y());
		radius = c.radius;
	}

	/** Constructs a circle through two given points and wth their midpoint as the center*/
	public Circle(Point p1, Point p2){
		this( Point.midPoint(p1, p2), p1.distance(p2)/2 );
	}

	/** Constructs a circle through 3 given points. */
	public Circle(Point a, Point b, Point c) {
		if (a.equals(b) || b.equals(c)) {
			center = Point.midPoint(a, c);
			radius = center.distance(a);
		}
		else 
			if (a.equals(c)) {
				center = Point.midPoint(a, b);
				radius = center.distance(a);
			}
			else {
				Line ab = Point.getBisector(a, b);
				Line bc = Point.getBisector(b, c);
				if (ab.isParallelWith(bc)) {
					Point p1 = a;
					Point p2 = b;
					if (a.distance(b) < a.distance(c)) p2 = c;
					if (b.distance(c) > p1.distance(p2)) p1 = b;
					center = Point.midPoint(p1, p2);
					radius = p1.distance(p2);
				}
				else {
					center = Line.getIntersection(ab, bc);
					if (center == null) throw new Error(a+" "+b+" "+c);
					radius = center.distance(b);
				}
			}
		}
	

	
//	/**
//	 * Creates circle through 2 given points and with center on given line
//	 */
//	public static Circle createEnclosingCircleOnLine(Line2d ab, Point b, Point c) {
//		Line2d bc = Point.getBisector(b, c);
//		center = Line2d.getIntersection(ab, bc);
//		radius = center.getDistance(c);
//	}
//  
	/** Creates smallest circle containing two given circles	 */
	public Circle(Circle c1, Circle c2){
		if(c1.contains(c2)){
			center = c1.center.clone();
			radius = c1.radius;
		}else if(c2.contains(c1)){
			center = c2.center.clone();
			radius = c2.radius;
		}else if(c1.center.equals(c2.center)){
			center = c1.center.clone();
			radius = Math.max(c1.radius, c2.radius);
		}else{
			Line l = new Line(c1.center, c2.center);
			center = l.getPoint(0.5*c1.center.distance(c2.center)-c1.radius/2+c2.radius/2);
			radius = center.distance(c1.center)+c1.radius;
		}
	}
	
	public Circle(Circle circle1, Circle circle2, Circle circle3) {
		Circle c1 = circle1;
		Circle c2 = circle2;
		Circle c3 = circle3;

		Circle tmp;
		if( (tmp=new Circle(c1,c2)).contains(c3)) { center = tmp.center; radius = tmp.radius; return;}
		if( (tmp=new Circle(c1,c3)).contains(c2)) { center = tmp.center; radius = tmp.radius; return;}
		if( (tmp=new Circle(c2,c3)).contains(c1)) { center = tmp.center; radius = tmp.radius; return;}

		Circle mec = ApolloniusSolver.solveApollonius(c1, c2, c3, 1, 1, 1);//68HOps
		this.center = mec.center;
		this.radius = mec.radius;
	}
	
	public Circle(Triangle tri){
		this(tri.getCorner(0), tri.getCorner(1), tri.getCorner(2));
	}


	public Point center() { return center; }
	public double getRadius() { return radius; }

	public void setCenter(Point p) { center = p; }
	public void setRadius(double r) { radius = r; }
	
	public void translate(Vector v) { center.addThis(v); }
	
	public Double enteringAngle(Point p, Circle C, boolean ccw) {
		double centerDist = center.distance(C.center);
		if (centerDist < C.radius + radius) {
			double alpha = Math.acos(((p.x()-center.x())*(C.center.x()-center.x()) + (p.y()-center.y())*(C.center.y()-center.y()))/(radius*centerDist)); 
			double beta = Math.acos((centerDist*centerDist + radius*radius - C.radius*C.radius)/(2*centerDist*radius)); 
			if (ccw) { 
				if (Point.rightTurn(C.center, center, p)) alpha = alpha + Math.PI - beta; else alpha = alpha - beta;
			}
			else {
				if (Point.leftTurn(C.center, center, p)) alpha = alpha + beta- 2*Math.PI; else alpha = beta - alpha; 
			}
			System.out.println("alpha = " + Functions.toDeg(alpha));
			return alpha;
		}
		return null;
	}
	
	public Double exitingAngle(Point p, Circle C, boolean ccw) {
		double centerDist = center.distance(C.center);
		Double alpha = Math.acos(((p.x()-center.x())*(C.center.x()-center.x()) + (p.y()-center.y())*(C.center.y()-center.y()))/(radius*centerDist));
		Double beta = Math.acos((centerDist*centerDist + radius*radius - C.radius*C.radius)/(2*centerDist*radius)); 
		if (ccw) {
			if (Point.rightTurn(C.center, center, p)) alpha = alpha - beta; else alpha = alpha + beta;
		}
		else {
			if (Point.leftTurn(C.center, center, p)) alpha = beta - alpha; else alpha = -alpha - beta;
		}
		System.out.println("alpha = " + Functions.toDeg(alpha));
		return alpha;
	}
	
	public Double enteringAngle(Point p, Line L, boolean ccw) {
		double dist = L.getDistance(p);
		if (dist < radius) {
			Point q = L.projectPoint(center);
			double alpha = Math.acos(((p.x()-center.x())*(q.x()-center.x()) + (p.y()-center.y())*(q.y()-center.y()))/(radius*dist)); 
			double beta = Math.acos(dist/radius);
			if (ccw) { 
				if (Point.rightTurn(q, center, p)) alpha = 2*Math.PI - alpha;
			}
			else {
				alpha = -alpha;
				beta = -beta;
				if (Point.leftTurn(q, center, p)) alpha = -2*Math.PI - alpha;
			}
			System.out.println("alpha = " + Functions.toDeg(alpha));
			System.out.println("beta  = " + Functions.toDeg(beta));
			return alpha - beta;
		}
		return null;
	}
	
	
	/** returns 0-2 intersections of this circle with another circle  */
	public Point[] intersections(Circle c){
		Point[] ret = new Point[2];
		double centerDist = center.distance(c.center); 
		if (centerDist < c.radius+radius) {
			// nested circles, no intersection
			if ((centerDist + radius < c.radius) || (centerDist + c.radius < radius)) return null;
			// nested circles, touch point from inside: this cricle is inside c
			if (centerDist + radius == c.radius) {
				ret[0] = c.center.add(c.center.vectorTo(center).scaleToLength(c.radius));
				ret[1] = null;
				return ret;
			}
			// nested circles, touch point from inside: c is inside this circle
			if (centerDist + c.radius == radius) {
				ret[0] = center.add(center.vectorTo(c.center).scaleToLength(radius));
				ret[1] = null;
				return ret;
			}
			// two intersections
			double dSq = c.center.distanceSquared(center); 
			double d = Math.sqrt(dSq);
			double RSq = radius*radius;
			double rSq = c.radius*c.radius;
			double tmp = dSq-rSq+RSq;

			double d1 = tmp/(2*d);
			double a = Math.sqrt( 4*dSq*RSq - tmp*tmp )/d; 
			
			Vector x = center.vectorTo(c.center).normalizeThis();
			Vector y = x.rotate90();
			x.multiplyThis(d1);
			y.multiplyThis(a/2);
			ret[0] = center.add(x).addThis(y);
			ret[1] = center.add(x).subtractThis(y);
			return ret;
		}
		// touch point from outside
		if (centerDist == c.radius+radius){
			ret[0] = center.add(center.vectorTo(c.center).scaleToLength(radius));
			ret[1] = null;
			return ret;
		}
		// circles are disjoint
		return null;
	}
	
	/** returns 0-2 intersections of this circle with a line, 
	 * source: http://mathworld.wolfram.com/Circle-LineIntersection.html
	 * some odd errors corrected 01-09-2013 - PW */
	public Point[] intersections(Line l) {
		double dx = l.n.y();
		double dy = -l.n.x();
		double x1 = l.p.x() - center.x();
		double y1 = l.p.y() - center.y();
		double x2 = x1 - dx;
		double y2 = y1 - dy;
		double dr2 = dx*dx + dy*dy;
		double D = x1*y2 - x2*y1;
		double delta = radius*radius*dr2 - D*D;
		if (delta < -Constants.EPSILON) return null;
		if (delta <  Constants.EPSILON) {
			Point[] ret = new Point[1];
			ret[0] = new Point(center.x() + D*dy/dr2, center.y() - D*dx/dr2);
			return ret;
		}
		Point[] ret = new Point[2];
		double deltaRoot = Math.sqrt(delta);
		if (dy < 0.0) {
			ret[0] = new Point(center.x() + (D*dy - dx*deltaRoot)/dr2, center.y() - (D*dx + dy*deltaRoot)/dr2);
			ret[1] = new Point(center.x() + (D*dy + dx*deltaRoot)/dr2, center.y() - (D*dx - dy*deltaRoot)/dr2);
		}
		else {
			ret[0] = new Point(center.x() + (D*dy + dx*deltaRoot)/dr2, center.y() - (D*dx - dy*deltaRoot)/dr2);
			ret[1] = new Point(center.x() + (D*dy - dx*deltaRoot)/dr2, center.y() - (D*dx + dy*deltaRoot)/dr2);
		}
		return ret;
	}

	
	public static Circle minimumEnclosingCircle_Welzl(List<Point> points){
		Point[] b = new Point[3];
		return findMEC(points.size(), points,0,b);
	}
	
	private static Circle findMEC(int n, List<Point> points, int m, Point[] b){
		Circle mec = new Circle(new Point(0,0),0);
		
		// Compute the Smallest Enclosing Circle defined by B
		if(m == 1)			mec = new Circle(b[0],0);
		else if(m == 2)		mec = new Circle(b[0], b[1]);
		else if(m == 3)		return new Circle( b[0], b[1], b[2]);
	
		// Check if all the points in p are enclosed
		for(int i=0; i<n; i++)	{
			if(!mec.contains(points.get(i))){
				// Compute B <--- B union P[i].
				b[m] = points.get(i);	
				mec = findMEC(i,points, m+1, b);// Recurse
			}
		}
		
		return mec;
	}
	
	double getPowerDistance(Point p) { return center.distanceSquared(p) - radius*radius; }
	
	double getPowerDistance(Circle c) { return getPowerDistance(c.center) -c.radius*c.radius; }
	
	/** returns TRUE if the interior of the circle (for a given eps reduction of the radius) is empty */
	public boolean isEmpty(List<Point> points, double eps) {
		for (Point p : points) if (contains(p, eps)) return false;
		return true;
	}

	
	boolean isOrthogonal(Circle c) { return getPowerDistance(c) == 0.0; }

	@Deprecated
	public static Circle minimumEnclosingCircle_bruteforce(List<Point> points){
		System.err.println("Warning: Please only use Circle.minimumEnclosingCircle_bruteforce for testing purposes!!");
		//System.out.println("MEC .. ");
		//for(Point p: points) System.out.println(p);

		//System.out.println("Checking pairs");
		double minRad = Double.POSITIVE_INFINITY;
		Circle minCircle = null;
		for(int i=0;i<points.size();i++){
			for(int j=i;j<points.size();j++){
				if(i==j) continue;
				Circle tmp = new Circle(points.get(i), points.get(j));
				boolean containsAll = true;
				for(Point p: points) if(!tmp.contains(p)) {
					containsAll=false;break;
				}
				//System.out.println(tmp+" .. contains all: "+containsAll);
				if(containsAll && tmp.radius<minRad){
					minRad = tmp.radius;
					minCircle = tmp;
				}
			}
		}

		for(int i=0;i<points.size();i++){
			for(int j=i;j<points.size();j++){
				for(int k=j;k<points.size();k++){
					if(i==j || i==k || j==k) continue;
					Circle tmp = new Circle(points.get(i), points.get(j), points.get(k));
					boolean containsAll = true;
					for(Point p: points) if(!tmp.contains(p)) {
						//System.out.println(" ! "+tmp.getCenter().getDistance(p)+" from "+p);
						containsAll=false;break;
					}
					if(containsAll && tmp.radius<minRad){
						minRad = tmp.radius;
						minCircle = tmp;
					}
					/*if(i==1 && j==2 && k==3){
						J3DScene scene = J3DScene.createJ3DSceneInFrame();
						scene.addShape(new Sphere3d(new Point3d(points.get(i).x, points.get(i).y, 0), 0.3), Color.BLACK);
						scene.addShape(new Sphere3d(new Point3d(points.get(j).x, points.get(j).y, 0), 0.3), Color.BLACK);
						scene.addShape(new Sphere3d(new Point3d(points.get(k).x, points.get(k).y, 0), 0.3), Color.BLACK);
						scene.addShape(new Sphere3d(new Point3d(points.get(0).x, points.get(0).y, 0), 0.3), Color.GRAY);
						scene.addShape(new Cylinder3d(new Point3d(tmp.center.x, tmp.center.y, 0), new Point3d(tmp.center.x, tmp.center.y, 0.1), tmp.radius), new Color(0,0,240, 100));
						scene.addShape(new Cylinder3d(new Point3d(minCircle.center.x, minCircle.center.y, 0), new Point3d(minCircle.center.x, minCircle.center.y, 0.1), minCircle.radius), new Color(0,0,240, 100));
					}*/
				}
			}

		}
		
		if(minCircle==null) throw new Error("minCircle not set .. "+points.size());
		return minCircle;
	}

	public boolean contains(Point p) {
//		return center.distance(p)<=(radius+0.0001);
		return center.distanceSquared(p) < radius*radius;
	}
	public boolean contains(Point p, double eps) {
		return center.distanceSquared(p) < radius*radius - eps;
	}
	
	public boolean contains(Circle c){
		return radius>=center.distance(c.center)+c.radius-0.000001;
	}

	public boolean onCircle(Point p) {
		return Math.abs(center.distanceSquared(p) - radius* radius) < Constants.EPSILON;
	}
	public String toString() { return toString(2); }
	
	public String toString(int dec) {
		return String.format("Circle[%s,%"+dec+"f]",center,radius);
	}

	public void toConsole() { toConsole(2); } 
	public void toConsole(int dec) { System.out.println(toString(dec)); }

	public void toScene(J2DScene scene) { scene.addShape(this, Color.black); }
	public void toScene(J2DScene scene, Color clr) { scene.addShape(this, clr); }


	public Point getCenter() {
		return center;
	} 
	
	public static void zerooneMove(Point A, Point B, Point C, Point D) {
		if (Point.leftTurn(A, B, C)) { oneoneMove(B, A, C, D); return; }
		double beta  = D.polarAngle(); double cosBeta  = Math.cos(beta);  double sinBeta  = Math.sin(beta);
		System.out.println("beta  = " + Functions.toDeg(beta));

		double aa  = A.getSquaredDistance(); 
		double bb  = B.getSquaredDistance(); 
		double cc  = C.getSquaredDistance();
		double dd  = D.getSquaredDistance(); double d = Math.sqrt(dd);
		
		double m41 = A.y()*(bb-cc) + B.y()*(cc-aa) + C.y()*(aa-bb);
		double m42 = A.x()*(bb-cc) + B.x()*(cc-aa) + C.x()*(aa-bb);
		double m43 = A.x()*(B.y()-C.y()) + B.x()*(C.y()-A.y()) + C.x()*(A.y()-B.y());
		double m44 = A.x()*(B.y()*cc -C.y()*bb) + B.x()*(C.y()*aa - A.y()*cc) + C.x()*(A.y()*bb - B.y()*aa);

		double coefSin = d*(sinBeta*m41 + cosBeta*m42);
		double coefCos = d*(sinBeta*m42 - cosBeta*m41);
		double coef    = m44 - dd*m43;
		Double[] roots = new Double[2];
		roots = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);

		if (roots != null) {
			J2DScene scene = J2DScene.createJ2DSceneInFrame();
			Point newD1 = D.clone(); newD1.rotation(roots[0]); newD1.toScene(scene, 0.05, Color.red);
			Point newD2 = D.clone(); newD2.rotation(roots[1]); newD2.toScene(scene, 0.05, Color.red);

			Circle circumABC = new Circle(A, C, B);
			Circle circumABD = new Circle(A, D, B);
			Circle circumACD = new Circle(A, C, D);
			Circle circumBCD = new Circle(B, C, D);
			LineSegment AC = new LineSegment(A, C);
			LineSegment BC = new LineSegment(B, C);
			LineSegment AB = new LineSegment(A, B);
			LineSegment AD = new LineSegment(A, D);
			LineSegment BD = new LineSegment(B, D);
			LineSegment CD = new LineSegment(C, D);
		
			circumABC.toScene(scene);
//			circumABD.toScene(scene);
//			circumACD.toScene(scene);
//			circumBCD.toScene(scene);
			circumABC.center.toScene(scene, 0.03, Color.blue);
			circumABD.center.toScene(scene, 0.03, Color.red);
			circumACD.center.toScene(scene, 0.03, Color.pink);
			circumBCD.center.toScene(scene, 0.03, Color.magenta);
			AC.toScene(scene);
			BC.toScene(scene);
			AD.toScene(scene);
			BD.toScene(scene);
			CD.toScene(scene);
			AB.toScene(scene);
			A.toScene(scene, 0.08, Color.black);
			B.toScene(scene, 0.08, Color.green);
			C.toScene(scene, 0.08, Color.blue);
			D.toScene(scene, 0.08, Color.red);
			Circle cir;
			double cos = Math.cos(0.005);
			double sin = Math.sin(0.005);
			for (double angle = 0; angle < 2*Math.PI; angle = angle + 0.005) {
				D.rotation(cos, sin);
				cir = new Circle(A, D, B);
				circumABD.center.set(cir.center);
				cir.center.toScene(scene, 0.03, Color.pink);			
				circumABD.radius = cir.radius;
				cir = new Circle(A, C, D);
				circumACD.center.set(cir.center);
				cir.center.toScene(scene, 0.03, Color.magenta);
				circumACD.radius = cir.radius;
				cir = new Circle(B, C, D);
				circumBCD.center.set(cir.center);
				cir.center.toScene(scene, 0.03, Color.red);
				circumBCD.radius = cir.radius;		
				try {
					Thread.sleep(50);
				} catch (InterruptedException e) {}
				scene.repaint();
			}
		}
	}

	
	
	public static void oneoneMove(Point A, Point B, Point C, Point D) {
		if (Point.leftTurn(A, B, C)) { oneoneMove(B, A, C, D); return; }
		double alpha = C.polarAngle(); double cosAlpha = Math.cos(alpha); double sinAlpha = Math.sin(alpha);
		double beta  = D.polarAngle(); double cosBeta  = Math.cos(beta);  double sinBeta  = Math.sin(beta);
		System.out.println("alpha = " + Functions.toDeg(alpha));
		System.out.println("beta  = " + Functions.toDeg(beta));

		double aa  = A.getSquaredDistance(); 
		double bb  = B.getSquaredDistance(); 
		double cc  = C.getSquaredDistance(); double c = Math.sqrt(cc);
		double dd  = D.getSquaredDistance(); double d = Math.sqrt(dd);
		
		double m11 = aa - bb;
		double m12 = A.y() - B.y();
		double m13 = A.y()*bb - B.y()*aa;
		double m22 = A.x() - B.x();
		double m23 = A.x()*bb - B.x()*aa;
		double m33 = A.x()*B.y() - A.y()*B.x();

		double coefSin = d*cosBeta*(m23 - cc*m22) + d*sinBeta*(m13 - cc*m12) + c*cosAlpha*(dd*m22 - m23) + c*sinAlpha*(dd*m12 - m13);
		double coefCos = d*cosBeta*(cc*m12 - m13) + d*sinBeta*(m23 - cc*m22) + c*cosAlpha*(m13 - dd*m12) + c*sinAlpha*(dd*m22 - m23);
		double coef    = d*c*m11*Math.sin(beta-alpha) + (cc - dd)*m33;
		
		Double[] roots = new Double[2];
		roots = Trigonometry.solveAsinXPlusBcosXplusC(coefSin, coefCos, coef);

		if (roots != null) {
			J2DScene scene = J2DScene.createJ2DSceneInFrame();
			Point newC1 = C.clone(); newC1.rotation(roots[0]); newC1.toScene(scene, 0.05, Color.blue);
			Point newC2 = C.clone(); newC2.rotation(roots[1]); newC2.toScene(scene, 0.05, Color.blue);
			Point newD1 = D.clone(); newD1.rotation(roots[0]); newD1.toScene(scene, 0.05, Color.red);
			Point newD2 = D.clone(); newD2.rotation(roots[1]); newD2.toScene(scene, 0.05, Color.red);

			Point origo = new Point(0,0);
			Circle circumABC = new Circle(A, C, B);
			Circle circumABD = new Circle(A, D, B);
			Circle circumACD = new Circle(A, C, D);
			Circle circumBCD = new Circle(B, C, D);
			LineSegment AC = new LineSegment(A, C);
			LineSegment BC = new LineSegment(B, C);
			LineSegment AB = new LineSegment(A, B);
			LineSegment AD = new LineSegment(A, D);
			LineSegment BD = new LineSegment(B, D);
			LineSegment CD = new LineSegment(C, D);
		
//		CO.toScene(scene, Color.blue);
//		DO.toScene(scene, Color.red);
			circumABC.toScene(scene);
			circumABD.toScene(scene);
			circumACD.toScene(scene);
			circumBCD.toScene(scene);
			circumABC.center.toScene(scene, 0.03, Color.blue);
			circumABD.center.toScene(scene, 0.03, Color.red);
			circumACD.center.toScene(scene, 0.03, Color.pink);
			circumBCD.center.toScene(scene, 0.03, Color.magenta);
			AC.toScene(scene);
			BC.toScene(scene);
			AD.toScene(scene);
			BD.toScene(scene);
			CD.toScene(scene);
			AB.toScene(scene);
			A.toScene(scene, 0.08, Color.black);
			B.toScene(scene, 0.08, Color.green);
			C.toScene(scene, 0.08, Color.blue);
			D.toScene(scene, 0.08, Color.red);
			Circle cir;
			double cos = Math.cos(0.005);
			double sin = Math.sin(0.005);
			for (double angle = 0; angle < 2*Math.PI; angle = angle + 0.005) {
				D.rotation(cos, sin);
				cir = new Circle(A, C, B);
				circumABC.center.set(cir.center);
				cir.center.toScene(scene, 0.03, Color.pink);			
				circumABC.radius = cir.radius;
				cir = new Circle(A, D, B);
				circumABD.center.set(cir.center);
				cir.center.toScene(scene, 0.03, Color.pink);			
				circumABD.radius = cir.radius;
				cir = new Circle(A, C, D);
				circumACD.center.set(cir.center);
				cir.center.toScene(scene, 0.03, Color.magenta);
				circumACD.radius = cir.radius;
				cir = new Circle(B, C, D);
				circumBCD.center.set(cir.center);
				cir.center.toScene(scene, 0.03, Color.red);
				circumBCD.radius = cir.radius;		
				try {
					Thread.sleep(50);
				} catch (InterruptedException e) {}
				scene.repaint();
			}
		}
	}
	
	public static void twozeroMove(Point A, Point B, Point C, Point D) {
		oneoneMove(A, D, B, C);
	}
	
/*
		
		Point origo = new Point(0,0);
		Circle circum = new Circle(s, r1, r2);
		LineSegment sr1 = new LineSegment(s, r1);
		LineSegment r1r2 = new LineSegment(r1, r2);
		LineSegment r2s = new LineSegment(r2, s);
		Circle r1c = new Circle(origo, origo.distance(r1));
		Circle r2c = new Circle(origo, origo.distance(r2));
		
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		r1c.toScene(scene);
		r2c.toScene(scene);
		circum.toScene(scene);
		sr1.toScene(scene);
		r1r2.toScene(scene);
		r2s.toScene(scene);
		s.toScene(scene, 0.05, Color.black);
		r1.toScene(scene, 0.05, Color.blue);
		r2.toScene(scene, 0.05, Color.blue);
		p.toScene(scene, 0.05, Color.red);
		Circle cir;
		circum.center.toScene(scene, 0.03, Color.pink);

		double cos = Math.cos(0.005);
		double sin = Math.sin(0.005);
		for (double angle = 0; angle < 2*Math.PI; angle = angle + 0.002) {
			r1.rotation(cos, sin);
			r2.rotation(cos, sin);
			cir = new Circle(s, r1, r2);
			circum.center = cir.center;
			circum.center.toScene(scene, 0.03, Color.pink);
			circum.radius = cir.radius;
			try {
				Thread.sleep(300);
			} catch (InterruptedException e) {}
			scene.repaint();
		}
	}
*/
	public static void main(String[] args) {
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		Point c = new Point(0.5,0.2);
		Circle cir = new Circle(c, 1);
		cir.toScene(scene);
		Point p = new Point(1,1.5);
		Line l = new Line(c, p);
		LineSegment seg = new LineSegment(c, p);
		seg.toScene(scene);
		Point[] inter =  cir.intersections(l);
		inter[0].toScene(scene, 0.02, Color.red);
		inter[1].toScene(scene, 0.02, Color.blue);
		System.out.println(inter);
//		zerooneMove(new Point(-8, 1), new Point(-6, 8), new Point(-1, 3), new Point(0.5, -5) );
//		oneoneMove(new Point(-8, 1), new Point(-6, 8), new Point(-1, 3), new Point(0.5, -5) );
//		oneoneMove(new Point(3, 1), new Point(2, 6), new Point(1, 6), new Point(2, 2) );
//		twozeroMove(new Point(1.35, 2.2), new Point(-0.1, 4.8), new Point(1.5, 1.8), new Point(0.7, 0.6) );
	}


}
