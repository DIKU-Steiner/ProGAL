package ProGAL.dataStructures;

import ProGAL.geom2d.Point;

public class SortToolPoint2dXY implements SortTool {
	public int compare(Object p1, Object p2) {
		if ((p1 instanceof Point) && (p2 instanceof Point)) {
			double x1 = ((Point)p1).x();
			double x2 = ((Point)p2).x();
			if (x1 < x2) return COMP_LESS;
			else { 
				if (x1 > x2) return COMP_GRTR; 
				else {
					double y1 = ((Point)p1).y();
					double y2 = ((Point)p2).y();
					if (y1 < y2) return COMP_LESS;
					else { 
						if (y1 > y2) return COMP_GRTR; else return COMP_EQUAL; 
					}
				}
			}
		}
		else throw SortTool.err1;
	}
}

