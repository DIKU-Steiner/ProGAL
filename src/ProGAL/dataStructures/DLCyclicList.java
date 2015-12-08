package ProGAL.dataStructures;

import java.util.Iterator;

public class DLCyclicList<T> implements Iterable<T>{
	protected DLNode<T> entry;
	protected int size;

	public DLCyclicList() {
		entry = null;
		size = 0;
	}
	
	/** Returns true if the list is empty */
	public boolean isEmpty() { return entry == null; }
	
	public DLNode<T> getFirst() { return entry; }
	
	/** Finds the element containing given object. Worst-case O(n)-time. */
	public DLNode<T> findNode(Object obj) {
		DLNode<T> elem = entry;
		while ((elem !=null) && (elem.obj != obj)) elem = elem.next;
		return elem;
	}
	
	public DLNode<T> getEntry(){ return entry; }
	public void setEntry(DLNode<T> n){	entry = n;	}
	
	/** Adds the object before the entry of the list. Worst-case O(1)-time. */
	public void pushBefore(T obj) {
		pushBefore(obj, entry);
	}
		
	
	/** Adds an object before the specified node in the list. Worst-case O(1)-time. */
	public void pushBefore(T obj, DLNode<T> n) {
		DLNode<T> nd;
		if(n==null) nd = new DLNode<T>(obj);
		else		nd = new DLNode<T>(obj,n.prev,n);
		if(entry==null) entry = nd;
		size++;
	}
	
	/** Adds an object after the specified node in the list. Worst-case O(1)-time. */
	public void pushAfter(T obj, DLNode<T> n) {
		DLNode<T> nd;
		if(n==null) nd = new DLNode<T>(obj);
		else		nd = new DLNode<T>(obj,n,n.next);
		if(entry==null) entry = nd;
		size++;
	}
	
	/** Deletes a node from the list and returns its object. Worst-case O(1)-time. */
	public T delete(DLNode<T> nd) {
		if(entry==null) throw new RuntimeException("Cannot delete from empty list");
		if(entry == nd) entry = nd.next;
		nd.next.prev = nd.prev;
		nd.prev.next = nd.next; 
		size--;
		return nd.obj;
	}
	
	/** Returns the number of elements in this cyclic list. */
	public int getSize() {
		return size;
	}
	
	public Iterator<T> iterator() {
		return new DLListIterator<T>(this);
	}

	public static void main(String[] args) {
		DLCyclicList<Integer> L = new DLCyclicList<Integer>();
		L.pushBefore(0, L.getEntry());
		L.pushBefore(1, L.getEntry());
		L.pushBefore(2, L.getEntry());
		L.pushBefore(3, L.getEntry());
		L.pushBefore(4, L.getEntry());
		L.pushBefore(5, L.getEntry());
		for(Integer i: L){
			System.out.println(i);
		}
	}

	
	
	public static class DLNode<T> {
		protected T obj;
		protected DLNode<T> prev;
		protected DLNode<T> next;
		
		public DLNode(T obj) { 
			this.obj = obj;
			prev = this;
			next = this; 
		}
		public DLNode(T obj, DLNode<T> prev, DLNode<T> next) {
			this.obj = obj;
			this.prev = prev;
			this.next = next;
			prev.next = this;
			next.prev = this;
		}
		
		/** Returns the object represented by the node */
		public T getObject() { return obj; }
		
		/** Returns previous node of the doubly linked list  */
		public DLNode<T> getPrev() { return prev; }
		
		/** Return next node of the doubly linked list */
		public DLNode<T> getNext() { return next; }
		
		/** Clears the node and returns its object */
		public T clear() {
			prev = null;
			next = null;
			return obj;
		}
	}
	
	private static class DLListIterator<T> implements java.util.Iterator<T> {

		private DLCyclicList<T> lst;
		private DLNode<T> current;
		
		public DLListIterator(DLCyclicList<T> lst) {
			this.lst = lst;
			current = null;
		}
		
		public boolean hasNext() { 
			if (lst.isEmpty()) return false; else return current != lst.entry; 
		}
		
		public T next() {
			if (hasNext()) {
				T ret = null;
				if (current == null) {
					ret = lst.entry.obj;
					current = lst.entry.next;
				}else{
					ret = current.obj;
					current = current.next;
				}
				return ret;
			}
			else return null; 
		}
		
		public void remove() {}
	}

}
