package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import ProGAL.geomNd.Point;

/**
 * Class that models a Range Tree of Points.
 * The dimension of the Range Tree is the dimension of the points.
 * @author Søren Lynnerup 11.01.2012
 */
public class RangeTree {
	
	RangeTreeNode root;
	int dimensions;
	
	/**
	 * Initializes the Range Tree.
	 * This either creates a 1-, 2- og N-dimensional Range Tree.
	 */
	public RangeTree(List<Point> points, boolean use_fc) {
		// Get the number of dimensions from the points.
		this.dimensions = points.get(0).getDimensions();
		
		// Initialize the root of the 0th dimension.
		if(this.dimensions == 1) {
			this.root = buildRangeTree1d(points, 0);
		} else if(this.dimensions == 2) {
			if(use_fc) {
				this.root = buildRangeTree2d(points, 0);
			} else {
				this.root = buildRangeTreeNd(points, 0, use_fc);
			}
		} else {
			this.root = buildRangeTreeNd(points, 0, use_fc);
		}
	}
	
	/**
	 * Returns the root of a 1-dimensional Range Tree.
	 * This method presorts the points before initializing the root.
	 */
	public static RangeTreeNode buildRangeTree1d(List<Point> points, int dimension) {
		// Create a new instance of the list.
		points = new ArrayList<Point>(points);
		Collections.sort(points, new PointDimensionComparator(dimension));
		return new RangeTreeNode1d(points, dimension);
	}
	
	/**
	 * Returns the root of a 2-dimensional Range Tree.
	 * This method presorts the points according to both of these dimensions 
	 * before initializing the root.
	 */
	public static RangeTreeNode buildRangeTree2d(List<Point> points, int dimension) {
		List< List<Point> > sortedpoints = new ArrayList< List<Point> >();
		
		// Presort the points according to each dimension.
		for(int d = 0; d < 2; d++) {
			sortedpoints.add(d, new ArrayList<Point>(points));
			Collections.sort(sortedpoints.get(d), new PointDimensionComparator(dimension + d));
		}
		
		return new RangeTreeNode2d(sortedpoints, dimension);
	}
	
	/**
	 * Returns the root of a N-dimensional Range Tree.
	 * This method presorts the points according to the dimension of the
	 * total tree this Range Tree represents.
	 */
	public static RangeTreeNode buildRangeTreeNd(List<Point> points, int dimension, boolean use_fc) {
		// Create a new instance of the list.
		points = new ArrayList<Point>(points);
		Collections.sort(points, new PointDimensionComparator(dimension));
		return new RangeTreeNodeNd(points, dimension, use_fc);
	}
	
	/**
	 * Query the Range Tree given a low and a high value for each dimension.
	 * The values in high does not necessary have to be higher than the values in low.
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
		
		return this.root.query(low, high);
	}
	
	/**
	 * Query the Range Tree given with the interval the two points span.
	 */
	public List<Point> query(Point p1, Point p2) {
		return this.query(p1.getCoords(), p2.getCoords());
	}

}
