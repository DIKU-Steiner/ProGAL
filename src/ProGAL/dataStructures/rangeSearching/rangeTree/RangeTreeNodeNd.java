package ProGAL.dataStructures.rangeSearching.rangeTree;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import ProGAL.geomNd.Point;

public class RangeTreeNodeNd extends RangeTreeNode {
	
	private int dimensions;
	
	protected RangeTreeNodeNd leftchild;
	protected RangeTreeNodeNd rightchild;
	protected RangeTreeNode assocstruct;
	
	/**
	 * Initialize the RangeTreeNode for N dimensions.
	 * When the associated structure reaches 2 dimensions, we use 
	 * a RangeTreeNode2d instead.
	 */
	public RangeTreeNodeNd(List< List<Point> > sortedpoints, int dimension) {
		// Define dimensions and the number of points we are working with.
		this.dimensions = sortedpoints.size();
		this.dimension = dimension;
		int numberofpoints = sortedpoints.get(0).size();
		
		// Find the median point.
		int medianindex = numberofpoints / 2;
		this.median = sortedpoints.get(0).get(medianindex);
		
		//System.out.println(this.dimension + "d " + (numberofpoints == 1 ? "Leaf" : "Median") + ": " + this.median);
		
		// Create the associated structure. Check if we have reach the last two dimensions.
		if(this.dimensions - (this.dimension + 1) == 2) {
			this.assocstruct = new RangeTreeNode2d(sortedpoints.subList(1, this.dimensions), this.dimension + 1);
		} else {
			this.assocstruct = new RangeTreeNodeNd(sortedpoints.subList(1, this.dimensions), this.dimension + 1);
		}
		
		// We have reached the bottom of the tree.
		if(numberofpoints == 1) {
			this.isleave = true;
			return;
		}
		
		// Take the first half of the 0d-values and put in the left child.
		List< List<Point> > leftpoints = new ArrayList< List<Point> >();
		leftpoints.add(0, sortedpoints.get(0).subList(0, medianindex));
		
		// Take the second half of the 0d-values and put in the right child.
		List< List<Point> > rightpoints = new ArrayList< List<Point> >();
		rightpoints.add(0, sortedpoints.get(0).subList(medianindex, numberofpoints));
		
		// Initialize the comparator for this first dimension.
		PointDimensionComparator comparator = new PointDimensionComparator(this.dimension);
		
		int dimensionid = 1;
		Iterator< List<Point> > iterD = sortedpoints.listIterator(dimensionid);
		
		// Split the rest of the dimensions in sortedpoints into the two children.
		while(iterD.hasNext()) {
			List<Point> list = iterD.next();
			
			leftpoints.add(dimensionid, new ArrayList<Point>());
			rightpoints.add(dimensionid, new ArrayList<Point>());
			
			List<Point> left = leftpoints.get(dimensionid);
			List<Point> right = rightpoints.get(dimensionid);
			
			Iterator<Point> iter = list.iterator();
			
			// If several points are equal to the median, some of them might have
			// be in one child and the other ones might be in the other child.
			// Remember the index where to put them, if they have to be inserted.
			List<Point> medianequals = new ArrayList<Point>();
			int leftmedianindex = 0;
			int rightmedianindex = 0;
			
			// Loop the points in some dimension.
			while(iter.hasNext()) {
				Point p = iter.next();
				
				// Add to the correct child.
				if(comparator.compare(p, this.median) < 0) {
					left.add(p);
				} else if(comparator.compare(p, this.median) > 0) {
					right.add(p);
				} else {
					// p is equal to the median. Remember where to put the medianequals.
					leftmedianindex = left.size();
					rightmedianindex = right.size();
					medianequals.add(p);
				}
			}
			
			// Put the right amount of medianequals into the children.
			int splitindex = medianindex - left.size();
			left.addAll(leftmedianindex, medianequals.subList(0, splitindex));
			right.addAll(rightmedianindex, medianequals.subList(splitindex, medianequals.size()));
			
			dimensionid++;
		}
		
		// Create the child nodes.
		this.leftchild = new RangeTreeNodeNd(leftpoints, this.dimension);
		this.rightchild = new RangeTreeNodeNd(rightpoints, this.dimension);
		
		super.leftchild = this.leftchild;
		super.rightchild = this.rightchild;
	}
	
	/**
	 * Return the points inside the given range.
	 */
	public List<Point> query(double[] low, double[] high) {
		// Find the splitting node on x-value.
		RangeTreeNodeNd splitnode = (RangeTreeNodeNd) findSplitNode(low[this.dimension], high[this.dimension], this.dimension);
		
		// If splitnode is a leaf, check if it should be reported.
		if(splitnode.isleave) {
			return splitnode.reportPoint(low, high);
		}
		
		List<Point> reportedpoints = new ArrayList<Point>();
		RangeTreeNodeNd n = splitnode.leftchild;
		
		// Run down left child.
		while(!n.isleave) {
			if(low[this.dimension] < n.median.get(this.dimension)) {
				// Report the right child
				reportedpoints.addAll(n.rightchild.query(low, high));
				n = n.leftchild;
			} else {
				n = n.rightchild;
			}
		}
		
		// Check if the leaf should be reported.
		reportedpoints.addAll(n.reportPoint(low, high));
		
		n = splitnode.rightchild;
		
		// Run down right child.
		while(!n.isleave) {
			if(high[this.dimension] >= n.median.get(this.dimension)) {
				// Report the left child
				reportedpoints.addAll(n.leftchild.query(low, high));
				n = n.rightchild;
			} else {
				n = n.leftchild;
			}
		}
		
		// Check if the leaf should be reported.
		reportedpoints.addAll(n.reportPoint(low, high));
		
		return reportedpoints;
	}
	
	/**
	 * Checks if the point should be reported.
	 */
	private List<Point> reportPoint(double[] low, double[] high) {
		List<Point> points = new ArrayList<Point>();
		
		if(this.median.get(this.dimension) >= low[this.dimension] && this.median.get(this.dimension) <= high[this.dimension]) {
			points.addAll(this.assocstruct.query(low, high));
		}
		
		return points;
	}

}
