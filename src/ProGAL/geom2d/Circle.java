package ProGAL.geom2d;

import java.util.List;

public class Circle implements Shape{
	protected Point center;
	protected double radius;

	/** Construct a circle with the given center and radius. */
	public Circle(Point center, double radius) {
		this.center = center;
		this.radius = radius;
	}


	
	public Circle(Point p1, Point p2){
		this( Point.midPoint(p1, p2), p1.distance(p2)/2 );
	}

	/** Creates circle through 3 given points. */
	public Circle(Point a, Point b, Point c) {
		Line ab = Point.getBisector(a, b);
		Line bc = Point.getBisector(b, c);
		if(a.equals(b) || b.equals(c)){
			center = Point.midPoint(a, c);
			radius = center.distance(a);
		}else if(a.equals(c)){
			center = Point.midPoint(a, b);
			radius = center.distance(a);
		}else if(ab.isParallelWith(bc)){
			Point p1 = a;
			Point p2 = b;
			if(a.distance(b)<a.distance(c)) p2 = c;
			if(b.distance(c)>p1.distance(p2)) p1 = b;
			center = Point.midPoint(p1, p2);
			radius = p1.distance(p2);
		}else{
			center = Line.getIntersection(ab, bc);
			if(center==null) throw new Error(a+" "+b+" "+c);
			radius = center.distance(c);
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


	public Point center() { return center; }
	public double getRadius() { return radius; }


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

	public boolean contains(Point p) {
		return center.distance(p)<=(radius+0.0001);
	}
	public boolean contains(Circle c){
		return radius>=center.distance(c.center)+c.radius-0.000001;
	}

	public String toString() { return toString(2); }
	
	public String toString(int dec) {
		return String.format("Circle[%s,%"+dec+"f]",center,radius);
	}

	public void toConsole() { toConsole(2); } 
	public void toConsole(int dec) { System.out.println(toString(dec)); }



	@Override
	public Point getCenter() {
		return center;
	} 

}
