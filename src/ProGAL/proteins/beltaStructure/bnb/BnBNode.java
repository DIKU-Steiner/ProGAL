package ProGAL.proteins.beltaStructure.bnb;

/**
 * A node in a branch and bound tree corresponding to a partial protein structure.
 * The protein structure consists of parts where a part can be either a beta-sheet
 * or a helix/coil segment. The <code>part</code> field indicates the index of the 
 * part corresponding to this node. The <code>structure</code> field indicates 
 * which structure the specified part has. To recreate the partial structure 
 * corresponding to this node, simply iterate to the root of the tree using the 
 * <code>parent</code> pointers and for each node set the indicated part to have 
 * the indicated structure.
 * 
 * This class is primarily a storage class and doesnt have any functionality itself.
 * 
 * @author R.Fonseca
 */
public class BnBNode implements Comparable<BnBNode>{
	
	protected final BnBNode parent;
	protected final int structure, part;
	protected double lowerBound;
	
	protected BnBNode(BnBNode parent, int structure){
		this.parent = parent;
		this.part = parent==null?-1:parent.part+1;
		this.structure = structure;
		this.lowerBound = Double.NEGATIVE_INFINITY;
	}

	@Override
	public int compareTo(BnBNode n) {
		return Double.compare(lowerBound, n.lowerBound);
	}
	
	public String toString(){
		return String.format("BnBNode[part %d, structure %d, lower bound %.2f]",part,structure,lowerBound);
	}
}
