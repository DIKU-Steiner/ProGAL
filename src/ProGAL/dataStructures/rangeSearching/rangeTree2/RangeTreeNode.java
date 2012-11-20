package ProGAL.dataStructures.rangeSearching.rangeTree2;

import java.util.List;

import ProGAL.geomNd.Point;

/**
 * Class that models a abstract node in a Range Tree.
 * This includes some instance variables as well at the 
 * findSplitNode method.
 * @author Søren 11.01.2012
 */
public abstract class RangeTreeNode {
	
	protected Point median;
	protected int dimension;
	
	protected RangeTreeNode leftchild;
	protected RangeTreeNode rightchild;
	
	protected boolean isleave;
	
	/**
	 * Searches the Range Tree to find the node where v1 and v2 are in different child trees.
	 */
	protected RangeTreeNode findSplitNode(double v1, double v2, int dimension) {
		RangeTreeNode n = this;
		
		// Run down the tree to find the split node.
		while(!n.isleave && (v2 < n.median.get(dimension) || v1 >= n.median.get(dimension))) {
			if(v2 < n.median.get(dimension)) {
				// The splitnode is down the left subtree.
				n = n.leftchild;
			} else {
				// The splitnode is down the right subtree.
				n = n.rightchild;
			}
		}
		
		return n;
	}

	
	/**
	 * Returns the points inside the given range.
	 */
	public abstract List<Point> query(double[] low, double[] high);
	protected abstract List<Point> query(double[] low, double[] high, List<Point> reportedpoints);
	
}
