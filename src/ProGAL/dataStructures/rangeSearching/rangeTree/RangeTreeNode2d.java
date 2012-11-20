package ProGAL.dataStructures.rangeSearching.rangeTree;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import ProGAL.geomNd.Point;

public class RangeTreeNode2d extends RangeTreeNode {
	
	protected RangeTreeNode2d leftchild;
	protected RangeTreeNode2d rightchild;
	protected List<RangeTreeAssocValue> assocstruct = new ArrayList<RangeTreeAssocValue>();
	
	/**
	 * Initialize the RangeTreeNode for 2 dimensions.
	 * We expect that sortedpoints.size() == 2.
	 */
	public RangeTreeNode2d(List< List<Point> > sortedpoints, int dimension) {
		// The number of points we are working with.
		int numberofpoints = sortedpoints.get(0).size();
		this.dimension = dimension;
		
		// Find the median point.
		int medianindex = numberofpoints / 2;
		this.median = sortedpoints.get(0).get(medianindex);
		
		//System.out.println(this.dimension + "d " + (numberofpoints == 1 ? "Leaf" : "Median") + ": " + this.median);
		
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
	 * We know that this is a 2D tree so we know the absolute dimension
	 */
	public List<Point> query(double[] low, double[] high) {
		// Find the splitting node on x-value.
		RangeTreeNode2d splitnode = (RangeTreeNode2d) findSplitNode(low[this.dimension], high[this.dimension], this.dimension);
		
		// If splitnode is a leaf, check if it should be reported.
		if(splitnode.isleave) {
			return splitnode.reportPoint(low, high);
		}
		
		// Find the associated value based on y-value.
		RangeTreeAssocValue value = splitnode.findValueInAssocStruct(low[this.dimension + 1]);
		
		// Instantiate a couple of variables.
		List<Point> reportedpoints = new ArrayList<Point>();
		RangeTreeNode2d n = splitnode.leftchild;
		RangeTreeAssocValue v;
		int newindex = value.leftindex;
		
		// Run down left child.
		while(!n.isleave) {
			
			// The way the associated structure is built, is that if there
			// is no value greater or equal to in a substructure, then we
			// have an out of bounds exception => no feasible y-values down there.
			try {
				v = n.assocstruct.get(newindex);
			} catch(IndexOutOfBoundsException e) {
				break;
			}
			
			if(low[this.dimension] < n.median.get(this.dimension)) {
				// Report the right child
				reportedpoints.addAll(n.rightchild.reportPoints(v.rightindex, high[this.dimension + 1]));
				n = n.leftchild;
				newindex = v.leftindex;
			} else {
				n = n.rightchild;
				newindex = v.rightindex;
			}
		}
		
		// Check if the leaf should be reported.
		reportedpoints.addAll(n.reportPoint(low, high));
		
		// Re-instantiate a could of variables.
		n = splitnode.rightchild;
		newindex = value.rightindex;
		
		// Run down right child.
		while(!n.isleave) {
			try {
				v = n.assocstruct.get(newindex);
			} catch(IndexOutOfBoundsException e) {
				break;
			}
			
			if(high[this.dimension] >= n.median.get(this.dimension)) {
				// Report the left child
				reportedpoints.addAll(n.leftchild.reportPoints(v.leftindex, high[this.dimension + 1]));
				n = n.rightchild;
				newindex = v.rightindex;
			} else {
				n = n.leftchild;
				newindex = v.leftindex;
			}
		}
		
		// Check if the leaf should be reported.
		reportedpoints.addAll(n.reportPoint(low, high));
		
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
	private List<Point> reportPoints(int associndex, double high) {
		List<Point> points = new ArrayList<Point>();
		
		Iterator<RangeTreeAssocValue> iter = this.assocstruct.listIterator(associndex);
		
		// Run through the associated list and report the points with y-value <= high.
		while(iter.hasNext()) {
			Point point = iter.next().point; 
			if(point.get(this.dimension + 1) <= high) {
				points.add(point);
			}
		}
		
		return points;
	}
	
	/**
	 * Checks if the point should be reported.
	 */
	private List<Point> reportPoint(double[] low, double[] high) {
		List<Point> points = new ArrayList<Point>();
		
		if(this.median.get(this.dimension) >= low[this.dimension] && this.median.get(this.dimension) <= high[this.dimension] && 
		this.median.get(this.dimension + 1) >= low[this.dimension + 1] && this.median.get(this.dimension + 1) <= high[this.dimension + 1]) {
			points.add(this.median);
		}
		
		return points;
	}
	
}
