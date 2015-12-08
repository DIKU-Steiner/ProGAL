package ProGAL.dataStructures.rangeSearching.rangeTree2;

import ProGAL.geomNd.Point;

/**
 * Class that models a value in the associated structure of a 
 * 2-dimensional Range Tree using Fractional Cascading.
 * 
 * This structure stores a point and two indices. The indices
 * point to the smallest point larger than or equal to this point,
 * in the associated structures in the two children of the node
 * that this associated structure is associated to.
 * 
 * @author Søren Lynnerup 11.01.2012
 */
public class RangeTreeAssocValue {
	
	public Point point;
	public int leftindex;
	public int rightindex;
	
	/**
	 * Instanciates a value object in a associated structure to a 2DRangeTreeNode.
	 */
	public RangeTreeAssocValue(Point point, int leftindex, int rightindex) {
		this.point = point;
		this.leftindex = leftindex;
		this.rightindex = rightindex;
	}
	
}
