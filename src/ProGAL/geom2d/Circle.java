package ProGAL.geom2d;

import java.util.List;

public class Circle {
	protected Point center;
	protected double radius;

	/** Construct a circle with the given center and radius. */
	public Circle(Point center, double radius) {
		this.center = center;
		this.radius = radius;
	}


	
	public Circle(Point p1, Point p2){
		this( Point.getMidPoint(p1, p2), p1.getDistance(p2)/2 );
	}

	/** Creates circle through 3 given points. */
	public Circle(Point a, Point b, Point c) {
		Line ab = Point.getBisector(a, b);
		Line bc = Point.getBisector(b, c);
		if(a.equals(b) || b.equals(c)){
			center = Point.getMidPoint(a, c);
			radius = center.getDistance(a);
		}else if(a.equals(c)){
			center = Point.getMidPoint(a, b);
			radius = center.getDistance(a);
		}else if(ab.isParallelWith(bc)){
			Point p1 = a;
			Point p2 = b;
			if(a.getDistance(b)<a.getDistance(c)) p2 = c;
			if(b.getDistance(c)>p1.getDistance(p2)) p1 = b;
			center = Point.getMidPoint(p1, p2);
			radius = p1.getDistance(p2);
		}else{
			center = Line.getIntersection(ab, bc);
			if(center==null) throw new Error(a+" "+b+" "+c);
			radius = center.getDistance(c);
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
//	/**
//	 * Creates smallest circle containing two given circles
//	 */
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
			center = l.getPoint(0.5*c1.center.getDistance(c2.center)-c1.radius/2+c2.radius/2);
			radius = center.getDistance(c1.center)+c1.radius;
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


	public Point getCenter() { return center; }
	public double getRadius() { return radius; }


//	/**
//	 * returns true if none of the given points is in the circle
//	 */
//	public boolean isEmpty(PointSet2d points) {
//		double rr = radius*radius;
//		int size = points.getSize(); 
//		int i = 0;
//		while (i < size) if (((Point)points.get(i++)).getSquaredDistance(center) < rr - 0.000000001) return false;
//		return true;
//	}
//
//	/**
//	 * returns true if none of the given points except points a, b and c is in the circle
//	 */
//	public boolean isEmpty(PointSet2d points, Point a, Point b, Point c) {
//		double rr = radius*radius;
//		Point p;
//		int size = points.getSize(); 
//		int i = 0;
//		while (i < size) {
//			p = (Point)points.get(i++);
//			if ((p != a) && (p != b) && (p != c) && (p.getSquaredDistance(center) < rr - 0.000000001)) return false;
//		}
//		return true;
//	}
//
//
//	/**
//	 * returns the secant of the circle on a given line
//	 * @param line
//	 * @return
//	 */
//	public Segment2d getIntersection(Line2d line) {
//		Point p1 = line.p;
//		Point p2 = line.getPoint(1.0);
////		Point p3 = center;
//		double dx = p2.x - p1.x;
//		double dy = p2.y - p1.y;
//		double ex = p1.x - center.x;
//		double ey = p1.y - center.y;
//		double a = dx*dx + dy*dy;
//		double b = 2*(dx*ex + dy*ey);
//		double c = center.x*center.x + center.y*center.y + p1.x*p1.x + p1.y*p1.y - 2*(center.x*p1.x + center.y*p1.y) - radius*radius;
//		double delta = b*b - 4*a*c; 
//		if (delta < 0) return null;
//		double u1, u2;
//		if (delta == 0) u1 = u2 = -b/(2*a);
//		else {
//			double sqr = Math.sqrt(delta)/(2*a);
//			u1 = -b + sqr;
//			u2 = -b - sqr;
//		}
//		return new Segment2d(new Point(p1.x + u1*dx, p1.y + u1*dy),
//							 new Point(p1.x + u2*dx, p1.y + u2*dy));
//	}
	
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

	private boolean contains(Point p) {
		return center.getDistance(p)<=(radius+0.0001);
	}
	public boolean contains(Circle c){
		return radius>=center.getDistance(c.center)+c.radius-0.000001;
	}

	public String toString() { return toString(2); }
	
	public String toString(int dec) {
		return String.format("Circle[%s,%"+dec+"f]",center,radius);
	}

	public void toConsole() { toConsole(2); } 
	public void toConsole(int dec) { System.out.println(toString(dec)); } 

}
