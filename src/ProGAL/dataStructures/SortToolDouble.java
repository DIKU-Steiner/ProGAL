package ProGAL.dataStructures;

public class SortToolDouble implements SortTool {
	public int compare(Object x1, Object x2) {
		if ((x1 instanceof Double) && (x2 instanceof Double)) {
			double d1 = ((Double)x1).doubleValue();
			double d2 = ((Double)x2).doubleValue();
			if (d1 < d2) return COMP_LESS;
			else { if (d1 > d2) return COMP_GRTR; else return COMP_EQUAL; }
		}
		else throw SortTool.err1;
	}
}

