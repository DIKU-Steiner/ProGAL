package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import ProGAL.geomNd.Point;

/**
 * This class models a simple list which implements the Query-method in O(n) time.
 * @author Søren Lynnerup 11.01.2012
 */
public class NaivRangeArray {
	
	private List<Point> points;
	private int dimensions;
	
	/**
	 * Initialize the NaivRangeArray.
	 */
	public NaivRangeArray(List<Point> points) {
		this.dimensions = points.get(0).getDimensions();
		this.points = points;
	}
	
	/**
	 * Return the points inside the given range.
	 */
	public List<Point> query(double[] low, double[] high) {
		// Make sure low <= high for every dimension.
		for(int dimension = 0; dimension < this.dimensions; dimension++) {
			if(low[dimension] > high[dimension]) {
				double temp = low[dimension];
				low[dimension] = high[dimension];
				high[dimension] = temp;
			}
		}
		
		// Initialize variables.
		Iterator<Point> iter = this.points.iterator();
		List<Point> reported = new ArrayList<Point>();
		
		// Run through every point and check if it is inside the interval.
		while(iter.hasNext()) {
			Point p = iter.next();
			boolean report = true;
			
			// Check that each dimension is inside the interval.
			for(int d = 0; d < this.dimensions; d++) {
				report = report && (p.get(d) >= low[d] && p.get(d) <= high[d]);
			}
			
			// Check if the point should be reported.
			if(report) {
				reported.add(p);
			}
		}
		
		return reported;
	}

}
