package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.LinkedList;
import java.util.List;

import ProGAL.geomNd.Point;

/**
 * Class that models a N-dimension Range Tree.
 * When the tree reaches the last two dimensions, a 2DRangeTree
 * is used for modeling these dimensions.
 * @author Soren Lynnerup 11.01.2012
 */
public class RangeTreeNodeNd extends RangeTreeNode {
	
//	protected RangeTreeNodeNd leftchild;
//	protected RangeTreeNodeNd rightchild;
	protected RangeTreeNode assocstruct;
	
	/**
	 * Initializes the RangeTreeNode for N dimensions.
	 */
	public RangeTreeNodeNd(List<Point> sortedpoints, int dimension, boolean use_fc) {
		// Define dimensions and the number of points we are working with.
		int dimensionsleft = sortedpoints.get(0).getDimensions() - dimension - 1;
		this.dimension = dimension;
		int numberofpoints = sortedpoints.size();
		
		// Find the median point.
		int medianindex = numberofpoints / 2;
		this.median = sortedpoints.get(medianindex);
		
		// Create the associated structure. Check if we have reach the last two dimensions.
		if(use_fc) {
			if(dimensionsleft == 2) {
				this.assocstruct = RangeTree.buildRangeTree2d(sortedpoints, this.dimension + 1);
			} else {
				this.assocstruct = RangeTree.buildRangeTreeNd(sortedpoints, this.dimension + 1, use_fc);
			}
		} else {
			if(dimensionsleft == 1) {
				this.assocstruct = RangeTree.buildRangeTree1d(sortedpoints, this.dimension + 1);
			} else {
				this.assocstruct = RangeTree.buildRangeTreeNd(sortedpoints, this.dimension + 1, use_fc);
			}
		}
		
		// We have reached the bottom of the tree.
		if(numberofpoints == 1) {
			this.isleave = true;
			return;
		}
		
		// Create the child nodes.
		this.leftchild = new RangeTreeNodeNd(sortedpoints.subList(0, medianindex), this.dimension, use_fc);
		this.rightchild = new RangeTreeNodeNd(sortedpoints.subList(medianindex, numberofpoints), this.dimension, use_fc);
		
		super.leftchild = this.leftchild;
		super.rightchild = this.rightchild;
	}
	
	/**
	 * Return the points inside the given range.
	 */
	public List<Point> query(double[] low, double[] high) {
		List<Point> reportedpoints = new LinkedList<Point>();
		return query(low,high,reportedpoints);
		
	}
	
	protected List<Point> query(double[] low, double[] high, List<Point> reportedpoints){
		// Find the splitting node on x-value.
		RangeTreeNodeNd splitnode = (RangeTreeNodeNd) findSplitNode(low[this.dimension], high[this.dimension], this.dimension);
		
//		List<Point> reportedpoints = new ArrayList<Point>();
		
		// If splitnode is a leaf, check if it should be reported.
		if(splitnode.isleave) {
			splitnode.reportPoint(reportedpoints, low, high);
			return reportedpoints;
		}
		
		RangeTreeNodeNd n = (RangeTreeNodeNd) splitnode.leftchild;
		
		// Run down left child.
		while(!n.isleave) {
			if(low[this.dimension] < n.median.get(this.dimension)) {
				// Report the right child
//				reportedpoints.addAll(n.rightchild.query(low, high));
				n.rightchild.query(low, high, reportedpoints);
				n = (RangeTreeNodeNd) n.leftchild;
			} else {
				n = (RangeTreeNodeNd) n.rightchild;
			}
		}
		
		// Check if the leaf should be reported.
		n.reportPoint(reportedpoints, low, high);
		
		n = (RangeTreeNodeNd) splitnode.rightchild;
		
		// Run down right child.
		while(!n.isleave) {
			if(high[this.dimension] >= n.median.get(this.dimension)) {
				// Report the left child
//				reportedpoints.addAll(n.leftchild.query(low, high));
				n.leftchild.query(low, high, reportedpoints);
				n = (RangeTreeNodeNd) n.rightchild;
			} else {
				n = (RangeTreeNodeNd) n.leftchild;
			}
		}
		
		// Check if the leaf should be reported.
		n.reportPoint(reportedpoints, low, high);
		
		return reportedpoints;
	}
	
	/**
	 * Checks if the point should be reported.
	 */
	private void reportPoint(List<Point> reportedpoints, double[] low, double[] high) {
		if(this.median.get(this.dimension) >= low[this.dimension] && this.median.get(this.dimension) <= high[this.dimension]) {
//			reportedpoints.addAll(this.assocstruct.query(low, high));
			this.assocstruct.query(low, high, reportedpoints);
		}
	}

}
