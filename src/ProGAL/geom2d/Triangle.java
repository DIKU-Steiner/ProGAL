package ProGAL.geom2d;


public class Triangle implements Shape{
	protected final Point[] points;

	public Triangle(Point p1, Point p2, Point p3) { points = new Point[]{p1,p2,p3}; }

	public Point getCorner(int i) { return points[i]; }

	public void setCorner(Point p, int i){ points[i] = p; }

	public double getAltitude(int i) {
		Line line  = new Line(points[(i+1)%3], points[(i+2)%3]);
		return line.getDistance(points[i]);
	}
	
	/** Calculate the area from the three side-lengths of a triangle */
	public static double calculateArea(double a, double b, double c){
		//http://www.mathsisfun.com/geometry/herons-formula.html
		double s = (a+b+c)*0.5;
		return Math.sqrt(s*(s-a)*(s-b)*(s-c));
	}
	
	/** 
	 * Calculate the height from the the endpoint of a and b down to the c-edge based on the 
	 * lengths of the edges.
	 */
	public static double calculateHeight(double a, double b, double c) {
		double A = calculateArea(a,b,c);
		return 2*A/c;
	}

	public Point getCenter() { return points[0].clone(); }

	/** returns the cosinus of the angle opposite to the vertex with vertex i */
	public double getCos(int i) {
		double a01 = points[0].distance(points[1]);
		double a12 = points[1].distance(points[2]);
		double a20 = points[2].distance(points[0]);
		if (i == 0) return (a01*a01 + a20*a20 - a12*a12)/(2*a01*a20);
		if (i == 1) return (a01*a01 + a12*a12 - a20*a20)/(2*a01*a12); 
		return (a12*a12 + a20*a20 - a01*a01)/(2*a12*a20);
	}

	public Circle getCircumCircle(){
		Line bisectorAB = Point.getBisector(points[0],points[1]);
		Line bisectorAC = Point.getBisector(points[0],points[2]);
		Point d = Line.getIntersection(bisectorAB, bisectorAC);
		return new Circle(d, d.distance(points[0]));
		//		return new Circle(points[0], points[1], points[2]);
	}

	public Point getSteinerPoint(){
		Vector v = points[0].vectorTo(points[1]);
		if(Point.leftTurn(points[0], points[1], points[2])){
			v.rotateThis(-Math.PI/3.0);
			Point eqPoint = points[0].add(v);
			Triangle eqTri = new Triangle(eqPoint, points[1], points[0]);
			if(!Point.rightTurn(eqPoint, points[0], points[2])) return points[0];
			if(!Point.leftTurn(eqPoint, points[1], points[2])) return points[1];
			Circle ccirc = eqTri.getCircumCircle();
			Line l = new Line(eqPoint, eqPoint.vectorTo(points[2]));
			Point[] intersections = ccirc.intersections(l);
			return intersections[1];
		}else{
			v.rotateThis(Math.PI/3.0);
			Point eqPoint = points[0].add(v);
			Triangle eqTri = new Triangle(eqPoint, points[1], points[0]);
			if(!Point.leftTurn(eqPoint, points[0], points[2])) return points[0];
			if(!Point.rightTurn(eqPoint, points[1], points[2])) return points[1];
			Circle ccirc = eqTri.getCircumCircle();
			Line l = new Line(eqPoint, eqPoint.vectorTo(points[2]));
			Point[] intersections = ccirc.intersections(l);
			return intersections[1];
		}
	}

	public boolean inCircumCircle(Point s) {
		Point p = points[0];
		Point q = points[1];
		Point r = points[2];
		double pD2 = p.x()*p.x() + p.y()*p.y();
		double qD2 = q.x()*q.x() + q.y()*q.y();
		double rD2 = r.x()*r.x() + r.y()*r.y();
		double sD2 = s.x()*s.x() + s.y()*s.y();
		return -s.x()*(p.y()*(qD2-rD2) + q.y()*(rD2-pD2) + r.y()*(pD2-qD2)) + 
				s.y()*(p.x()*(qD2-rD2) + q.x()*(rD2-pD2) + r.x()*(pD2-qD2)) - 
				sD2*(p.x()*(q.y()-r.y()) + q.x()*(r.y()-p.y()) + r.x()*(p.y()-q.y())) + 
				p.x()*(q.y()*rD2-r.y()*qD2) + p.y()*(r.x()*qD2-q.x()*rD2) + pD2*(q.x()*r.y()-r.x()*q.y())  < 0.0;
	}

	public static void main(String[] args) {
		Triangle tr = new Triangle(new Point(-3.2,-1.5), new Point(-1,3), new Point(3.1,0));
		if (tr.inCircumCircle(new Point(0,-4))) System.out.println("TRUE"); else System.out.println("FALSE");
	}

	@Override
	public boolean contains(Point p) {
		boolean l1 = Point.leftTurn(points[0], points[1], p); 
		boolean l2 = Point.leftTurn(points[1], points[2], p);
		boolean l3 = Point.leftTurn(points[2], points[0], p);
		return l1==l2 && l1==l3;
	}


}
