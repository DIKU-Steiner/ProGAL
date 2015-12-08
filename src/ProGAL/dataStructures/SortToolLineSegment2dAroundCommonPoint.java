package ProGAL.dataStructures;

import ProGAL.geom2d.Point;
import ProGAL.geom2d.LineSegment;

public class SortToolLineSegment2dAroundCommonPoint implements SortTool {
	Point p;
	double x;
	double y;
	public SortToolLineSegment2dAroundCommonPoint(Point p) {
		super();
		this.p = p;
		x = p.x();
		y = p.y();
	}
	public int compare(Object s1, Object s2) {
		if ((s1 instanceof LineSegment) && (s2 instanceof LineSegment)) {
			Point p1 = ((LineSegment)s1).getB();
			Point p2 = ((LineSegment)s2).getB();
			double x1 = p1.x(); 
			double x2 = p2.x();
			double y1 = p1.y();
			double y2 = p2.y();
			if (y1 > y) {
				if (y2 < y) return COMP_LESS; 
				if (y2 > y) {
					if (Point.leftTurn(p, p1, p2)) return COMP_LESS;
					if ((Point.collinear(p, p1, p2)) &&
							((x1-x)*(x1-x) + (y1-y)*(y1-y) < (x2-x)*(x2-x) + (y2-y)*(y2-y))) return COMP_LESS; else return COMP_GRTR;
				}
				if (x2 >= x) return COMP_GRTR; else return COMP_LESS;
			}
			if (y1 < y) {
				if (y2 >= y) return COMP_GRTR;
				if (Point.leftTurn(p, p1, p2)) return COMP_LESS;
				if ((Point.collinear(p, p1, p2)) &&
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