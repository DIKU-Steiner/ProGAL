package ProGAL.dataStructures;

import ProGAL.geom2d.Point;
import ProGAL.geom2d.*;

public class SortToolPoint2dAroundPoint implements SortTool{
	Point p;
	double x;
	double y;
	public SortToolPoint2dAroundPoint(Point p) {
		super();
		this.p = p;
		x = p.x();
		y = p.y();
	}
	public int compare(Object p1, Object p2) {
		double x1 = ((Point)p1).x() - x; 
		double x2 = ((Point)p2).x() - x;
		double y1 = ((Point)p1).y() - y;
		double y2 = ((Point)p2).y() - y;

		if ((p1 instanceof LineSegment) && (p2 instanceof LineSegment)) {
			if (y1 > y) {
				if (y2 < y) return COMP_LESS; 
				if (y2 > y) {
					if (Point.leftTurn(p, (Point)p1, (Point)p2)) return COMP_LESS;
					if ((Point.collinear(p, (Point)p1, (Point)p2)) &&
							((x1-x)*(x1-x) + (y1-y)*(y1-y) < (x2-x)*(x2-x) + (y2-y)*(y2-y))) return COMP_LESS; else return COMP_GRTR;
				}
				if (x2 >= x) return COMP_GRTR; else return COMP_LESS;
			}
			if (y1 < y) {
				if (y2 >= y) return COMP_GRTR;
				if (Point.leftTurn(p, (Point)p1, (Point)p2)) return COMP_LESS;
				if ((Point.collinear(p, (Point)p1, (Point)p2)) &&
						((x1-x)*(x1-x) + (y1-y)*(y1-y) < (x2-x)*(x2-x) + (y2-y)*(y2-y))) return COMP_LESS; else return COMP_GRTR;
			}
			if (x1 >= x) {
				if ((y2 == y) && ((x1 < x2) || (x2 < x))) return COMP_LESS; else return COMP_GRTR;
			}
			if (y2 > y) return COMP_GRTR;						
			if (y2 < y) return COMP_LESS;
			if (x2 > x1)  return COMP_GRTR; else return COMP_LESS;		
		}
		else throw SortTool.err1;
	}
}
