package ProGAL.dataStructures;

import ProGAL.geom2d.Point;


public class SortToolPoint2dAroundOrigo implements SortTool{
	public int compare(Object p1, Object p2) {
		if ((p1 instanceof Point) && (p2 instanceof Point)) {
			double x1 = ((Point)p1).x(); 
			double x2 = ((Point)p2).x();
			double y1 = ((Point)p1).y();
			double y2 = ((Point)p2).y();
			if (y1 > 0.0) {
				if (y2 < 0) return COMP_LESS; 
				if (y2 > 0) {
					if (Point.leftTurn(Point.origo, (Point)p1, (Point)p2)) return COMP_LESS;
					if ((Point.collinear(Point.origo, (Point)p1, (Point)p2)) &&
							(x1*x1 + y1*y1 < x2*x2 + y2*y2)) return COMP_LESS; else return COMP_GRTR;
				}
				if (x2 >= 0) return COMP_GRTR; else return COMP_LESS;
			}
			if (y1 < 0.0) {
				if (y2 >= 0.0) return COMP_GRTR;
				if (Point.leftTurn(Point.origo, (Point)p1, (Point)p2)) return COMP_LESS;
				if ((Point.collinear(Point.origo, (Point)p1, (Point)p2)) &&
						(x1*x1 + y1*y1 < x2*x2 + y2*y2)) return COMP_LESS; else return COMP_GRTR;
			}
			if (x1 >= 0.0) {
				if ((y2 == 0.0) && ((x1 < x2) || (x2 < 0))) return COMP_LESS; else return COMP_GRTR;
			}
			if (y2 > 0.0) return COMP_GRTR;						
			if (y2 < 0.0) return COMP_LESS;
			if (x2 > x1)  return COMP_GRTR; else return COMP_LESS;		
		}
		else throw SortTool.err1;
	}

}

