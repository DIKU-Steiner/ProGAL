package ProGAL.dataStructures.rangeSearching.rangeTree;

import ProGAL.geomNd.Point;

public class RangeTreeAssocValue {
	
	public Point point;
	public int leftindex;
	public int rightindex;
	
	/**
	 * Instanciates an value object in a associated structure.
	 */
	public RangeTreeAssocValue(Point point, int leftindex, int rightindex) {
		this.point = point;
		this.leftindex = leftindex;
		this.rightindex = rightindex;
	}
	
}
