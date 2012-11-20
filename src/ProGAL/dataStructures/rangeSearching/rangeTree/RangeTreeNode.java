package ProGAL.dataStructures.rangeSearching.rangeTree;

import java.util.List;

import ProGAL.geomNd.Point;

public abstract class RangeTreeNode {
	
	protected Point median;
	protected int dimension;
	
	protected RangeTreeNode leftchild;
	protected RangeTreeNode rightchild;
	
	protected boolean isleave;
	
	/**
	 * Search the tree to find the node where v1 and v2 are in different child trees.
	 */
	protected RangeTreeNode findSplitNode(double v1, double v2, int dimension) {
		RangeTreeNode n = this;
		
		// Run down the tree to find the split node.
		while(!n.isleave && (v2 < n.median.get(dimension) || v1 >= n.median.get(dimension))) {
			if(v2 < n.median.get(dimension)) {
				n = n.leftchild;
			} else {
				n = n.rightchild;
			}
		}
		
		return n;
	}
	
	/**
	 * Return the points inside the given range.
	 */
	public abstract List<Point> query(double[] low, double[] high);
	
}
