package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geomNd.Point;

/**
 * Class that models a 1-dimensional Range Tree.
 * @author Søren Lynnerup 11.01.2012
 */
public class RangeTreeNode1d extends RangeTreeNode {
	
	protected RangeTreeNode1d leftchild;
	protected RangeTreeNode1d rightchild;
	protected RangeTreeNode assocstruct;
	
	/**
	 * Initialize the Range Tree given a list of presorted points.
	 */
	public RangeTreeNode1d(List<Point> sortedpoints, int dimension) {
		// Define dimensions and the number of points we are working with.
		this.dimension = dimension;
		int numberofpoints = sortedpoints.size();
		
		// Find the median point.
		int medianindex = numberofpoints / 2;
		this.median = sortedpoints.get(medianindex);
		
		// We have reached the bottom of the tree.
		if(numberofpoints == 1) {
			this.isleave = true;
			return;
		}
		
		// Create the child nodes.
		this.leftchild = new RangeTreeNode1d(sortedpoints.subList(0, medianindex), this.dimension);
		this.rightchild = new RangeTreeNode1d(sortedpoints.subList(medianindex, numberofpoints), this.dimension);
		
		super.leftchild = this.leftchild;
		super.rightchild = this.rightchild;
	}
	
	/**
	 * Queries the Range Tree given a interval.
	 * 
	 * This method runs in O(lgn + k) time.
	 */
//	public List<Point> query(double[] low, double[] high) {
		
	public List<Point> query(double[] low, double[] high) {
			List<Point> reportedpoints = new ArrayList<Point>();
			return query(low,high,reportedpoints);
			
		}
		
	protected List<Point> query(double[] low, double[] high, List<Point> reportedpoints){
		// Find the splitting node on x-value.
		RangeTreeNode1d splitnode = (RangeTreeNode1d) findSplitNode(low[this.dimension], high[this.dimension], this.dimension);
		
//		List<Point> reportedpoints = new ArrayList<Point>();
		
		// If splitnode is a leaf, check if it should be reported.
		if(splitnode.isleave) {
			if(splitnode.median.get(this.dimension) >= low[this.dimension] && splitnode.median.get(this.dimension) <= high[this.dimension]) {
				reportedpoints.add(splitnode.median);
			}
			
			return reportedpoints;
		}
		
		RangeTreeNode1d n = splitnode.leftchild;
		
		// Run down left child.
		while(!n.isleave) {
			if(low[this.dimension] < n.median.get(this.dimension)) {
				// Report the right child
				n.rightchild.reportSubtree(reportedpoints);
				n = n.leftchild;
			} else {
				n = n.rightchild;
			}
		}
		
		// Check if the leaf should be reported.
		if(n.median.get(this.dimension) >= low[this.dimension] && n.median.get(this.dimension) <= high[this.dimension]) {
			reportedpoints.add(n.median);
		}
		
		n = splitnode.rightchild;
		
		// Run down right child.
		while(!n.isleave) {
			if(high[this.dimension] >= n.median.get(this.dimension)) {
				// Report the left child
				n.leftchild.reportSubtree(reportedpoints);
				n = n.rightchild;
			} else {
				n = n.leftchild;
			}
		}
		
		// Check if the leaf should be reported.
		if(n.median.get(this.dimension) >= low[this.dimension] && n.median.get(this.dimension) <= high[this.dimension]) {
			reportedpoints.add(n.median);
		}
		
		return reportedpoints;
	}
	
	private void reportSubtree(List<Point> reportedpoints) {
		if(this.isleave) {
			reportedpoints.add(this.median);
		} else {
			this.leftchild.reportSubtree(reportedpoints);
			this.rightchild.reportSubtree(reportedpoints);
		}
	}
}
