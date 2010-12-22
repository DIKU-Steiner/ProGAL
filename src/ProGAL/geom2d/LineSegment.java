package ProGAL.geom2d;

public class LineSegment {
	protected Point a,b;
	
	public LineSegment(Point a, Point b){
		this.a = a;
		this.b = b;
	}

	public Point getA(){ return a; }
	public Point getB(){ return b; }
	
	public double getLength(){
		return a.getDistance(b);
	}
	
	public LineSegment clone(){
		return new LineSegment(a.clone(), b.clone());
	}
}
