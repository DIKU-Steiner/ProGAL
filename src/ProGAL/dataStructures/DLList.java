package ProGAL.dataStructures;

public class DLList {
	protected DLNode first;
	protected DLNode last;
	protected Integer size;

	public DLList() {
		first = null;
		last = null;
		size = null;
	}
	
	/*
	 * returns true if the list is empty
	 */
	public boolean isEmpty() { return first == null; }
	
	public DLNode getFirst() { return first; }
	
	/*
	 * finds the element containing given object
	 */
	public DLNode findNode(Object obj) {
		DLNode elem = first;
		while ((elem !=null) && (elem.obj != obj)) elem = elem.next;
		return elem;
	}
	
	/*
	 * adds an object at the end of the list
	 */
	public void add(Object obj) {
		DLNode nd = new DLNode(obj,last,null);
		if (last != null)  last.next = nd;
		if (first == null) first = nd;
		last = nd;
		if (size != null) size++;
	}
	
	/* 
	 * adds an object at the front of the list
	 */
	public void push(Object obj) {
		DLNode nd = new DLNode(obj,null,first);
		if (first != null) first.prev = nd;
		if (last == null) last = nd;
		first = nd;
		if (size != null) size++;
	}
	
	/*
	 * deletes a node from the list and returns its object
	 */
	public Object delete(DLNode nd) {
		if (first == nd) first = nd.next;
		if (last  == nd) last = nd.prev;
		if (nd.next != null) nd.next.prev = nd.prev; 
		if (nd.prev != null) nd.prev.next = nd.next;
		if (size != null) size--;
		return nd.clear();
	}
	
	public Integer getSize() {
		if (size != null) return size;
		else {
			DLListIterator iter = new DLListIterator(this);
			size = 0;
			while (iter.hasNext()) {
				size++;
				iter.next();
			}
			return size;
		}
	}

	public static void main(String[] args) {
		System.out.println("here");
		DLList L = new DLList();
		Integer i1 = new Integer(1);
		Integer i2 = new Integer(2);
		L.add(i1);
		L.add(i2);
		System.out.println(L.getSize());

	}

	
	class DLNode {
		protected Object obj;
		protected DLNode prev;
		protected DLNode next;
		
		public DLNode(Object obj) { 
			this.obj = obj;
			prev = null;
			next = null; 
		}
		public DLNode(Object obj, DLNode prev, DLNode next) {
			this.obj = obj;
			this.prev = prev;
			this.next = next;
		}
		
		/** Returns the object represented by the node */
		public Object getObject() { return obj; }
		
		/** Returns previous node of the doubly linked list  */
		public DLNode getPrev() { return prev; }
		
		/** Return next node of the doubly linked list */
		public DLNode getNext() { return next; }
		
		/** Clears the node and returns its object */
		public Object clear() {
			prev = null;
			next = null;
			return obj;
		}
	}
	
	class DLListIterator implements java.util.Iterator<Object> {

		private DLList lst;
		private DLNode current;
		
		public DLListIterator(DLList lst) {
			this.lst = lst;
			current = null;
		}
		
		public boolean hasNext() { 
			if (lst.isEmpty()) return false; else return current != lst.last; 
		}
		
		public Object next() {
			if (hasNext()) {
				if (current == null) current = lst.first; else current = current.next;
				return current.obj;
			}
			else return null; 
		}
		
		public void remove() {}
	}
}
