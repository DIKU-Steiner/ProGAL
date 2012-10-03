package ProGAL.geom2d;

import java.util.ArrayList;
import java.util.List;

public class Polygon implements Shape {
	private final List<Point> corners;
	
	public Polygon(List<Point> corners){
		this.corners = new ArrayList<Point>(corners);
	}
	
	public Polygon(Point[] corners){
		this.corners = new ArrayList<Point>();
		for(Point p:corners) this.corners.add(p);
	}

	public List<Point> getCorners(){
		return corners;
	}
	
	public boolean isConvex(){
		if(corners.size()<4) return true;
		
		Point p0=corners.get(0);
		Point p1=corners.get(1);
		Point p2=corners.get(2);
		boolean ccw = Point.leftTurn(p0,p1,p2);
		
		for(int i=1;i<corners.size();i++){
			p0=p1;
			p1=p2;
			p2=corners.get((i+2)%corners.size());
			if(ccw!=Point.leftTurn(p0, p1, p2)) return false;
		}
		return true;
	}
	
	@Override
	public Point getCenter() {
		Vector v = new Vector(0,0);
		for(Point p: corners) v.addThis(p);
		return new Point(v.x()/corners.size(), v.y()/corners.size());
	}

}
