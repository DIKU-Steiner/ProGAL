package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import ProGAL.geomNd.Point;

/**
 * Class that models a 2-dimension Range Tree using Fractional Cascading.
 * @author Søren Lynnerup 11.01.2012
 */
public class RangeTreeNode2d extends RangeTreeNode {
	
//	protected RangeTreeNode2d leftchild;
//	protected RangeTreeNode2d rightchild;
	protected List<RangeTreeAssocValue> assocstruct = new ArrayList<RangeTreeAssocValue>();
	
	/**
	 * Initializes the RangeTreeNode for 2 dimensions.
	 * We expect that sortedpoints.size() == 2.
	 */
	public RangeTreeNode2d(List< List<Point> > sortedpoints, int dimension) {
		// The number of points we are working with.
		int numberofpoints = sortedpoints.get(0).size();
		this.dimension = dimension;
		
		// Find the median point.
		int medianindex = numberofpoints / 2;
		this.median = sortedpoints.get(0).get(medianindex);
		
		// Take the first half of the x-values and put in the left child.
		List< List<Point> > leftpoints = new ArrayList< List<Point> >();
		leftpoints.add(0, sortedpoints.get(0).subList(0, medianindex));
		leftpoints.add(1, new ArrayList<Point>());
		
		// Take the second half of the x-values and put in the right child.
		List< List<Point> > rightpoints = new ArrayList< List<Point> >();
		rightpoints.add(0, sortedpoints.get(0).subList(medianindex, numberofpoints));
		rightpoints.add(1, new ArrayList<Point>());
		
		// Initialize the comparator for this first dimension (x-values).
		PointDimensionComparator comparator = new PointDimensionComparator(this.dimension);
		
		// We have reached the bottom of the tree.
		if(numberofpoints == 1) {
			this.isleave = true;
			this.assocstruct.add(new RangeTreeAssocValue(this.median, 0, 0));
			return;
		}
		
		// Initialize some iterators.
		Iterator<Point> iterY = sortedpoints.get(1).iterator();
		List<Point> leftY = leftpoints.get(1);
		List<Point> rightY = rightpoints.get(1);
		
		// Split the y-values onto the two children.
		while(iterY.hasNext()) {
			Point p = iterY.next();
			
			int leftassocindex = leftY.size();
			int rightassocindex = rightY.size();
			
			// Add to the correct child.
			if(comparator.compare(p, this.median) < 0) {
				leftY.add(p);
			} else {
				rightY.add(p);
			}
			
			// Add to the associated structure.
			this.assocstruct.add(new RangeTreeAssocValue(p, leftassocindex, rightassocindex));
		}
		
		// Create the child nodes.
		this.leftchild = new RangeTreeNode2d(leftpoints, this.dimension);
		this.rightchild = new RangeTreeNode2d(rightpoints, this.dimension);
		
		super.leftchild = this.leftchild;
		super.rightchild = this.rightchild;
	}
	
	
	/**
	 * Return the points inside the given range.
	 */
	public List<Point> query(double[] low, double[] high) {
		List<Point> reportedpoints = new ArrayList<Point>();
		return query(low,high,reportedpoints);
		
	}
	
	
	protected List<Point> query(double[] low, double[] high, List<Point> reportedpoints){
//	public List<Point> query(double[] low, double[] high) {
		// Find the splitting node on x-value.
		RangeTreeNode2d splitnode = (RangeTreeNode2d) findSplitNode(low[this.dimension], high[this.dimension], this.dimension);
		
//		List<Point> reportedpoints = new ArrayList<Point>();
		
		// If splitnode is a leaf, check if it should be reported.
		if(splitnode.isleave) {
			splitnode.reportPoint(reportedpoints, low, high);
			
			return reportedpoints;
		}
		
		RangeTreeAssocValue value;
		
		// **: The way the associated structure is built, is that if there
		// is no value greater or equal to in a substructure, then we
		// have an out of bounds exception => no feasible y-values down there.
		try {
			value = splitnode.findValueInAssocStruct(low[this.dimension + 1]);
		} catch(IndexOutOfBoundsException e) {
			return new ArrayList<Point>();
		}
		
		// Instantiate a couple of variables.
		RangeTreeNode2d n = (RangeTreeNode2d) splitnode.leftchild;
		RangeTreeAssocValue v;
		int newindex = value.leftindex;
		
		// Run down left child until we reach low[dimension] reporting right children on the way.
		while(!n.isleave) {
			
			// As commented above marked with **.
			try {
				v = n.assocstruct.get(newindex);
			} catch(IndexOutOfBoundsException e) {
				break;
			}
			
			if(low[this.dimension] < n.median.get(this.dimension)) {
				// Report the right child
				((RangeTreeNode2d) n.rightchild).reportPoints(reportedpoints, v.rightindex, high[this.dimension + 1]);
				n = (RangeTreeNode2d) n.leftchild;
				newindex = v.leftindex;
			} else {
				n = (RangeTreeNode2d) n.rightchild;
				newindex = v.rightindex;
			}
		}
		
		// Check if the leaf should be reported.
		n.reportPoint(reportedpoints, low, high);
		
		// Re-instantiate a could of variables.
		n = (RangeTreeNode2d) splitnode.rightchild;
		newindex = value.rightindex;
		
		// Run down right child until we reach high[dimension] reporting left children on the way.
		while(!n.isleave) {
			try {
				v = n.assocstruct.get(newindex);
			} catch(IndexOutOfBoundsException e) {
				break;
			}
			
			if(high[this.dimension] >= n.median.get(this.dimension)) {
				// Report the left child
				((RangeTreeNode2d) n.leftchild).reportPoints(reportedpoints, v.leftindex, high[this.dimension + 1]);
				n = (RangeTreeNode2d) n.rightchild;
				newindex = v.rightindex;
			} else {
				n = (RangeTreeNode2d) n.leftchild;
				newindex = v.leftindex;
			}
		}
		
		// Check if the leaf should be reported.
		n.reportPoint(reportedpoints, low, high);
		
		return reportedpoints;
	}
	
	/**
	 * Performs a binary search in the associated structure to find
	 * the smallest value larger or equal to the parameter.
	 */
	private RangeTreeAssocValue findValueInAssocStruct(double value) {
		int min = 0;
		int max = this.assocstruct.size();
		int medianindex;

		// Perform a binary search.
		while(min != max) {
			medianindex = min + (max - min) / 2;
			RangeTreeAssocValue obj = this.assocstruct.get(medianindex);
			
			// Branch left or right.
			if(value <= obj.point.get(this.dimension + 1)) {
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
		
		return this.assocstruct.get(min);
	}
	
	/**
	 * Reports the points from the associated list inside the y-value interval.
	 */
	private void reportPoints(List<Point> reportedpoints, int associndex, double high) {
		Iterator<RangeTreeAssocValue> iter = this.assocstruct.listIterator(associndex);
		
		// Run through the associated list and report the points with y-value <= high.
		while(iter.hasNext()) {
			Point point = iter.next().point; 
			if(point.get(this.dimension + 1) <= high) {
				reportedpoints.add(point);
			}
		}
	}
	
	/**
	 * Checks if the point should be reported.
	 */
	private void reportPoint(List<Point> reportedpoints, double[] low, double[] high) {
		if(this.median.get(this.dimension) >= low[this.dimension] && this.median.get(this.dimension) <= high[this.dimension] && 
		this.median.get(this.dimension + 1) >= low[this.dimension + 1] && this.median.get(this.dimension + 1) <= high[this.dimension + 1]) {
			reportedpoints.add(this.median);
		}
	}
	
}
