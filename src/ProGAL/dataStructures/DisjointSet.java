package ProGAL.dataStructures;

public class DisjointSet {

	public DisjointSet() {}
	
	public DSNode makeSet(Object object) { return new DSNode(object); }
	
	public DSNode find(DSNode node) {
		DSNode nd = node;
		while (nd != nd.parent) nd = nd.parent;
		compress(node, nd);
		return nd;
	}
	
	public DSNode union(DSNode nd1, DSNode nd2) {
		if (nd1.rank > nd2.rank) {
			nd2.parent = nd1;
			return nd1;
		}
		else {
			if (nd1.rank < nd2.rank) {
				nd1.parent = nd2;
				return nd2;
			}
			else {
				nd2.parent = nd1;
				nd1.rank = nd1.rank + 1;
				return nd1;
			}
		}
	}
	
	private void compress(DSNode node, DSNode root) {
		DSNode nd = node; 
		DSNode ndNext = nd.parent;
		while (nd != ndNext) {
			nd.parent = root;
			nd = ndNext;
			ndNext = nd.parent;
		}
	}
	
	public static class DSNode {
		protected Object object;
		protected DSNode parent;
		protected int rank; 
		
		public DSNode(Object object) {
			this.object = object;
			parent = this;
			rank = 0;
		}
	}
}
