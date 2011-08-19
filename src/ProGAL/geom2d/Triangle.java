package ProGAL.geom2d;


public class Triangle {
	protected final Point[] points;
	
	public Triangle(Point p1, Point p2, Point p3){
		points = new Point[]{p1,p2,p3};
	}
	
	public Point getCorner(int i){
		return points[i];
	}
	public void setCorner(Point p, int i){
		points[i] = p;
	}
}
