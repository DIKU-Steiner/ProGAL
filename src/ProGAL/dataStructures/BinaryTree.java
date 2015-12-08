package ProGAL.dataStructures;

public class BinaryTree {

	protected Node root;

	private class Node {
		private Object object;
		private Node left;
		private Node right;
		private Node father;
			
		protected Node(Object object) {
			this.object = object;
			left = right = father = null;
		}
		
		/** returns the object stored in the node */
		protected Object getObject() { return object; }
		
		/** returns TRUE if this node is a child of node p */
		protected boolean isChildOf(Node p) { return (p.left == this) || (p.right == this); }
		
		/** returns TRUE if this node is a leaf */
		protected boolean isLeaf() {return (left == null) && (right == null); }
		
		/** returns TRUE if this node is the root of binary tree */
		protected boolean isRoot() { return father == null; }
	}
	
	public BinaryTree() {}
		
	public void insert(Object object) {
		Node nd = new Node(object);
	}
	
}
