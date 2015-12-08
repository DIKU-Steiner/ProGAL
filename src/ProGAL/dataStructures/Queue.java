package ProGAL.dataStructures;


public class Queue {
	protected QueueNode first, last;
	protected int size;
	
	private class QueueNode {
		protected Object obj = null;
		protected QueueNode next = null;
		
		public QueueNode(Object obj, QueueNode prev) { this.obj = obj; }
	}
	
	public Queue() { first = null; size = 0; }
	
	public void clear() { while (!isEmpty()) pop(); }
	
	public boolean isEmpty() {return first == null; }
	
	public void push(Object obj) {
		QueueNode nd = new QueueNode(obj, last);
		if (isEmpty()) first = nd; else last.next = nd;
		last = nd;
		size++;
	}
	
	public Object pop() {
		if (isEmpty()) return null;
		Object obj = first.obj;
		first = first.next;
		size--;
		if (isEmpty()) last = null;
		return obj;
	}
	
	public Object peek() { if (isEmpty()) return null; else return first.obj; }
	
	public int getSize() { return size; }
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	}
	

}
