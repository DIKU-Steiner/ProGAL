package ProGAL.geom2d;

public class LineSegment implements Shape{
	protected Point a,b;
	
	public LineSegment(Point a, Point b){
		this.a = a;
		this.b = b;
	}

	public Point getA(){ return a; }
	public Point getB(){ return b; }
	
	
	
	public double getLength(){
		return a.distance(b);
	}
	
	public LineSegment clone(){
		return new LineSegment(a.clone(), b.clone());
	}

	public Point getCenter() {
		return Point.midPoint(a, b);
	}
	
	public double distance(Point p){
		Vector v = a.vectorTo(b);
		Vector vP = a.vectorTo(p);
		double t = v.dot(vP)/v.getSquaredLength();
		return a.add(v.multiplyThis(Math.min(Math.max(0,t), 1))).distance(p);
	}
}
