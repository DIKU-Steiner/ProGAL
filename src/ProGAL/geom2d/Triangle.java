package ProGAL.geom2d;


public class Triangle implements Shape{
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

	public Point getCenter() {
		return points[0].clone();
	}
	
	public Circle getCircumCircle(){
		Line bisectorAB = Point.getBisector(points[0],points[1]);
		Line bisectorAC = Point.getBisector(points[0],points[2]);
		Point d = Line.getIntersection(bisectorAB, bisectorAC);
		return new Circle(d, d.distance(points[0]));
//		return new Circle(points[0], points[1], points[2]);
	}
}
