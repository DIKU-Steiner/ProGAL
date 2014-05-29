package ProGAL.geom2d;

/** 
 * A line-swept circle.
 */
public class LSC implements Shape{
	protected LineSegment segment;
	protected double radius;
	
	public LSC(LineSegment segment, double radius){
		this.segment = segment;
		this.radius = radius;
	}
	
	public LSC(Point p1, Point p2, double radius){
		this(new LineSegment(p1,p2),radius);
	}
	
	public void setRadius(double rad){	this.radius = rad;	}
	
	public LineSegment getSegment(){ return segment; }
	public double getRadius(){ return radius; }
	
	public double getArea(){
		return (Math.PI*radius+2*segment.getLength())*radius;
	}

	public Point getCenter() {
		return segment.getCenter();
	}
	
//	/**
//	 * Create an LSC given a canonical representation of three points and the direction and radius of the LSC. 
//	 * The canonical representation of three points (x1,y1), (x2,y2), (x3,y3) means that they are translated 
//	 * and scaled such that y1=y2=0, y3=1 and -x1=x2. It is assumed that the direction is 'near (1,0,0)' in the 
//	 * sense that P1 and P2 should be the points that span the hemicircles of the LSC
//	 */
//	public static LSC createLSCFromDirAndRad(double x2, double x3, Vector dir, double rad){
//		return null;
//	}

	public LSC clone(){
		return new LSC(segment.clone(), radius);
	}

	@Override
	public boolean contains(Point p) {
		double dist = segment.distance(p);
		return dist<=radius;
	}
	
}
