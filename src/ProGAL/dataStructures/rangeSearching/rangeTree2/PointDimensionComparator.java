package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.Comparator;

import ProGAL.geomNd.Point;

/**
 * Class used to compare two points at a given dimension.
 * If two points are equal at the given dimension, the one point
 * with a value at a dimension (running 0..d) smaller than the other
 * is the smaller point.
 * 
 * @author Søren Lynnerup 11.01.2012
 */
public class PointDimensionComparator implements Comparator<Point> {
	int dimension;
	
	/**
	 * Initialize the comparator at a given dimension.
	 */
	public PointDimensionComparator(int dimension) {
		this.dimension = dimension;
	}

	/**
	 * Compares two points.
	 * Note that if p1[d] == p2[d] then we run through all of the dimensions,
	 * from 0..d, and the one point with a value smaller than the other will
	 * be smaller than the other point.
	 */
	public int compare(Point p1, Point p2) {
		if(p1.getCoord(this.dimension) < p2.getCoord(this.dimension)) {
			// p1 < p2
			return -1;
		} else if(p1.getCoord(this.dimension) > p2.getCoord(this.dimension)) {
			// p1 > p2
			return 1;
		} else {
			// p1 == p2
			int dimensions = p1.getDimensions();
			
			// Check the other coordinates from 0..dimensions
			for(int i = 0; i < dimensions; i++) {
				if(p1.getCoord(i) < p2.getCoord(i)) {
					// p1 < p2
					return -1;
				} else if(p1.getCoord(i) > p2.getCoord(i)) {
					// p1 < p2
					return 1;
				}
			}
			
			return 0;
		}
	}
}
