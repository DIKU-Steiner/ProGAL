package ProGAL.dataStructures.rangeSearching;

import java.util.List;

import ProGAL.geomNd.Point;

public interface RangeSearchDataStructure {
//	public void addPoint(Point p);
	public List<Point> query(Point low, Point high);
}
