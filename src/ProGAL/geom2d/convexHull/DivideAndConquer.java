package ProGAL.geom2d.convexHull;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom2d.Point;
import ProGAL.geom2d.Polygon;

public class DivideAndConquer {

	
	public static Polygon getConvexHull(List<Point> points) {
		return new Polygon(recursiveCH(points));
	}
	
	private static List<Point> recursiveCH(List<Point> points) {
		if(points.size()==3) {
			if(Point.rightTurn(points.get(0),points.get(1),points.get(2)))
				points.add(points.remove(1));
			return points;
		}else if(points.size()<3){
			return points;
		}

		List<Point> half = new LinkedList<Point>();
		int l = points.size();
		while(points.size()>l/2)
			half.add(points.remove(points.size()-1));

		List<Point> ch1 = recursiveCH(points);
		List<Point> ch2 = recursiveCH(half);
		
		
		return null;
	}
	
	private static List<Point> merge(List<Point> ch1, List<Point> ch2){
		List<Point> ret = new ArrayList<Point>();
		
		int p1 = lowest(ch1);
		int p2 = -1;
//		for(int)
		
		
		return ret;
	}
	
	private static double angle(Point p0, Point p1){
		return Math.atan2(p1.y()-p0.y(), p1.x()-p0.x());
	}
	
	private static int lowest(List<Point> ch){
		Point min = null;
		for(Point p: ch) 
			if(min==null || p.y()<min.y())
				min = p;
		
		return ch.indexOf(min);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println(180/Math.PI*angle(new Point(1,1), new Point(2,1)));
		System.out.println(180/Math.PI*angle(new Point(1,1), new Point(2,2)));
		System.out.println(180/Math.PI*angle(new Point(1,1), new Point(1,2)));
		System.out.println(180/Math.PI*angle(new Point(1,1), new Point(1,-1)));
		System.out.println(180/Math.PI*angle(new Point(1,1), new Point(0,1)));
	}

}
