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

	@Override
	public Point getCenter() {
		return Point.midPoint(a, b);
	}
}
