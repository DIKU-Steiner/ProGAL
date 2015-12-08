package ProGAL.steiner.bnb;

import java.util.Arrays;

/** 
 * A node in the BnB-tree. 
 * @author RFonseca
 */
public class Node {
	final int depth;
	final Node parent;
	final int edgeSplit;
	public final int siteInserted;
	public double lowerBound;
	
	public Node(){
		this.depth = 0;
		this.parent = null;
		this.edgeSplit = -1;
		this.siteInserted = 2;
	}
	
	public Node(Node parent, int splitEdge, int siteInserted){
		this.parent = parent;
		this.depth = parent.depth+1;
		this.edgeSplit = splitEdge;
		this.siteInserted = siteInserted;
	}
	
	public String toString(){
		Node n = this;
		int[] kVec = new int[this.depth];
		while(n.parent!=null){ // && n!=prevNode ){
			kVec[n.depth-1] = n.edgeSplit;
			n = n.parent;
		}
		return String.format("Node[lb=%.3f,depth=%d,kVec=%s]",lowerBound, depth, Arrays.toString(kVec));
	}
	
}
