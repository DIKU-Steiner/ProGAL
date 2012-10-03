package ProGAL.geom2d.convexHull;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Stack;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.Polygon;
import ProGAL.geom2d.viewer.J2DScene;
import ProGAL.math.Randomization;

public class GrahamScan {

	public static void main(String[] args) {
		ArrayList<Point> points = new ArrayList<Point>();
		for(int i=0;i<10;i++){
			points.add(new Point(Randomization.randBetween(0.0, 1.0),Randomization.randBetween(0.0, 1.0)));
		}
		Polygon pol = getConvexHull(points);

		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		for(Point p: points) scene.addShape(new Circle(p,0.01));
		scene.addShape(pol);

	}


	public static Polygon getConvexHull(List<Point> points) {
		int length = points.size();
		final Point pivot = lowestPoint(points);

		ArrayList<Point> sorted = new ArrayList<Point>(points);
		Collections.sort(sorted, new PolarAngleComparator(pivot));
		sorted.add(0, pivot);

		Stack<Point> stack = new Stack<Point>();
		stack.push(sorted.get(length - 1));
		stack.push(pivot);
		int i = 1;
		while (i < length) {
			Point pt1 = stack.pop();
			Point pt2 = stack.peek();
			stack.push(pt1);
			if ( Point.leftTurn(pt1, pt2, sorted.get(i)) ) {
				stack.push(sorted.get(i));
				i++;
			} else {
				stack.pop();
			}
		}
		ArrayList<Point> convex = new ArrayList<Point>(stack);
		convex.remove(convex.size() - 1);
		return new Polygon(convex);
	}

	private static double distance(Point a, Point b) {
		double dx = a.x() - b.x(), dy = a.y() - b.y();
		return dx * dx + dy * dy;
	}


	private static Point lowestPoint(List<Point> points){
		double minY = Double.POSITIVE_INFINITY;
		Point ret = null;
		for(Point p: points) 
			if(p.y()<minY){
				ret = p;
				minY = p.y();
			}
		return ret;
	}
	private static class PolarAngleComparator implements Comparator<Point>{
		private final Point pivot;
		PolarAngleComparator(Point pivot){
			this.pivot = pivot;
		}
		public int compare(Point o1, Point o2) {
			if (o1.equals(o2)) {
				return 0;
			}
			if (angle_cmp(pivot, o1, o2)) {
				return 1;
			} else {
				return -1;
			}	
		}
		private double area(Point a, Point b, Point c) {
			return a.x() * b.y() - a.y() * b.x() + b.x() * c.y() - b.y() * c.x() + c.x() * a.y() - c.y() * a.y();
		} 

		private boolean angle_cmp(Point pivot, Point a, Point b) {
			if (area(pivot, a, b) == 0) {
				return distance(pivot, a) < distance(pivot, b);
			}
			double d1x = a.x() - pivot.x(), d1y = a.y() - pivot.y();
			double d2x = b.x() - pivot.x(), d2y = b.y() - pivot.y();
			return (Math.atan2((double) d1y, (double) d1x) - Math.atan2((double) d2y, (double) d2x)) < 0;
		}
	}
}
