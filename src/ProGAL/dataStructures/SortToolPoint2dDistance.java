package ProGAL.dataStructures;

import ProGAL.geom2d.Point;

public class SortToolPoint2dDistance implements SortTool {
	
	public int compare(Object p1, Object p2) {
		if ((p1 instanceof Point) && (p2 instanceof Point)) {
			double d1 = ((Point)p1).getSquaredDistance();
			double d2 = ((Point)p2).getSquaredDistance();
			if (d1 < d2) return COMP_LESS; else {
				if (d1 > d2) return COMP_GRTR; else return COMP_EQUAL;
			}
		}
		else throw SortTool.err1;
	}
}
