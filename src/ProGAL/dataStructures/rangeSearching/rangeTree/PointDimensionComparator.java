package ProGAL.dataStructures.rangeSearching.rangeTree;

import java.util.Comparator;

import ProGAL.geomNd.Point;

public class PointDimensionComparator implements Comparator<Point> {
	int dimension;
	
	public PointDimensionComparator(int dimension) {
		this.dimension = dimension;
	}

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
			
			// Check the other coordinates from 0..this.dimension
			for(int i = 0; i < this.dimension; i++) {
				if(p1.getCoord(i) < p2.getCoord(i)) {
					// p1 < p2
					return -1;
				} else if(p1.getCoord(i) > p2.getCoord(i)) {
					// p1 < p2
					return 1;
				}
			}
			
			// Check the other coordinates from this.dimension..dimensions
			for(int i = this.dimension; i < dimensions; i++) {
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
