package ProGAL.dataStructures;

public class SortToolString implements SortTool {
	public int compare(Object x1, Object x2) {
		if ((x1 instanceof String) && (x2 instanceof String)) {
			int c = ((String)x1).compareTo((String)x2);
			if (c < 0) return COMP_LESS; else if (c > 0) return COMP_GRTR; else return COMP_EQUAL;
		}
		else throw SortTool.err1;
	}
}

