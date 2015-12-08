package ProGAL.dataStructures;

import ProGAL.geom2d.LineSegment;

public class SortToolLineSegment2dByLength implements SortTool {
	public int compare(Object p1, Object p2) {
		if ((p1 instanceof LineSegment) && (p2 instanceof LineSegment)) {
			double w1 = ((LineSegment)p1).getSquaredLength();
			double w2 = ((LineSegment)p2).getSquaredLength();
			if (w1 < w2) return COMP_LESS;
			if (w1 > w2) return COMP_GRTR; 
			return COMP_EQUAL; 
		}
		else throw SortTool.err1;
	}
}
