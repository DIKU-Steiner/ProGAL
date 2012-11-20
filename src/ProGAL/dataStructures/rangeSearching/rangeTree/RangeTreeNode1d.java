package ProGAL.dataStructures.rangeSearching.rangeTree;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geomNd.Point;

public class RangeTreeNode1d extends RangeTreeNode {
	
	private List<Point> points;
	
	public RangeTreeNode1d(List< List<Point> > sortedpoints) {
		this.points = sortedpoints.get(0);
	}

	public List<Point> query(double[] low, double[] high) {
		int min = 0;
		int max = this.points.size();
		int medianindex;

		// Perform a binary search to find the smallest number larger or equal to 'low'.
		while(min != max) {
			medianindex = min + (max - min) / 2;
			Point point = this.points.get(medianindex);
			
			// Branch left or right.
			if(low[0] <= point.get(0)) {
				max = medianindex;
			} else {
				// Because of rounding, we check if we are at a stall.
				if(min == medianindex) {
					min = max;
				} else {
					min = medianindex;
				}
			}
		}

		List<Point> reportedpoints = new ArrayList<Point>();
		
		// Run through the array until the coordinates are too high, and report the points
		while(min < this.points.size() && points.get(min).get(0) <= high[0]) {
			reportedpoints.add(points.get(min));
			min++;
		}
		
		return reportedpoints;
	}
}
