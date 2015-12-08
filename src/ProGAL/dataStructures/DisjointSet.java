package ProGAL.dataStructures;

import java.util.HashMap;
import java.util.Map;

public class DisjointSet {
	private final Map<Object,DSNode> nodeMap = new HashMap<Object,DSNode>();
	private int size;
	
	public DisjointSet() {
		this.size = 0;
	}
	
	public DSNode makeSet(Object object) {
		nodeMap.put(object, new DSNode(object));
		this.size++;
		return nodeMap.get(object); 
	}
	
	public DSNode find(Object o) {
		if(nodeMap.containsKey(o)) return find(nodeMap.get(o)); 
		return null;
	}
	
	public DSNode find(DSNode node) {
		DSNode nd = node;
		while (nd != nd.parent) nd = nd.parent;
		compress(node, nd);
		return nd;
	}
	
	public DSNode union(DSNode nd1, DSNode nd2) {
		this.size--;
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
	
	public int nrSets() {
		return size;
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
