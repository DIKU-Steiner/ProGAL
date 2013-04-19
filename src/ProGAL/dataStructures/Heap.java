package ProGAL.dataStructures;

	/**
	 * A heap-based priority queue, without any concurrency control
	 * (i.e., no blocking on empty/full states).
	 * This class provides the data structure mechanics for BoundedPriorityQueue.
	 * <p>
	 * The class currently uses a standard array-based heap, as described
	 * in, for example, Sedgewick's Algorithms text. All methods
	 * are fully synchronized. In the future,
	 * it may instead use structures permitting finer-grained locking.
	 * <p>[<a href="http://gee.cs.oswego.edu/dl/classes/EDU/oswego/cs/dl/util/concurrent/intro.html"> Introduction to this package. </a>]
	 **/

	public class Heap  {
	  protected Object[] nodes;  // tree nodes, packed into an array
	  protected int count = 0;   // number of used slots
	  protected SortTool tool;   // for ordering

	  /* ****************************************************************************
	   * Creates a heap with the given initial capacity and with the given SortTool
	   * @exception IllegalArgumentException if capacity less or equal to zero.
	   ******************************************************************************/
	  public Heap(int capacity, SortTool tool) 
	  throws IllegalArgumentException {
		  if (capacity < 0) capacity = 0;
		  nodes = new Object[capacity];
		  this.tool = tool;
	  }

	  /*
	   * Creates a heap with the given capacity, and relying on natural ordering.
	   */
	  public Heap(int capacity) { this(capacity, null); }
	  
	  /*
	   * Creates a heap of the given set in O(n) time, with the given SortTool
	   */
	  public Heap(Set<Object> set, SortTool tool) {
		  count = set.getSize();
		  nodes = new Object[count];
		  this.tool = tool;
		  System.arraycopy(set.getElements(), 0, nodes, 0, count);		  
		  for (int i = parent(count-1); i >= 0; i--) siftDown(i);
	  }

	  // indexes of heap parents and children
	  protected final int parent(int k) { return (k - 1) / 2;  }
	  protected final int left(int k)   { return 2 * k + 1; }
	  protected final int right(int k)  { return 2 * (k + 1); }

	  public Object[] getObjects() { return nodes; }
	  
	  /**
	   * returns true if the heap is empty.
	   * @return
	   */
	  public boolean isEmpty() {  return count == 0; }
	  public int getSize() { return count; }
	  
	  /**
	   * returns the i-th object in the binary heap. This is not a standard operation
	   * @param i
	   * @return
	   */
	  public Object getItem(int i) { return nodes[i]; }
	  
	  public void setItem(int i, Object x) {
		  nodes[i] = x;
	  }
	  
	  public void siftUp(int k) {
		  int q = k;
		  int p;
		  Object temp;
		  while ((q > 0) && (tool.compare(nodes[q], nodes[parent(q)]) < 0)) { 
			  p = parent(q);
			  temp = nodes[q];
			  nodes[q] = nodes[p]; 
			  nodes[p] = temp;
			  q = p; 
		  }
	  }
	  
	  public void siftDown(int k) {
		  int q = k;
		  int l = left(q);
		  int r;
		  int min;
		  while (l < count) {
			  r = right(q);
			  if ((r == count) || (tool.compare(nodes[l], nodes[r]) < 0)) min = l; else min = r; 
			  if (tool.compare(nodes[min], nodes[q]) < 0) {
				  Object temp = nodes[q];
				  nodes[q] = nodes[min];
				  nodes[min] = temp;
				  q = min;
				  l = left(q);
			  }
			  else break;
		  }
	  }
	  
	  /*
	   * insert an element, resize if necessary
	   */
	  public synchronized void insert(Object x) {
		  if (count >= nodes.length) {
			  int newcap =  3 * nodes.length / 2 + 1;
			  Object[] newnodes = new Object[newcap];
			  System.arraycopy(nodes, 0, newnodes, 0, nodes.length);
			  nodes = newnodes;
		  }
		  nodes[count++] = x;
		  siftUp(count-1);
	  }
	    
	  

	  /**
	   * Return and remove least element, or null if empty
	   **/

	  public synchronized Object extract() {
	    if (count < 1) return null;
	    Object x = nodes[0];
	    --count;
	    nodes[0] = nodes[count];
	    nodes[count] = null;
	    siftDown(0);
	    return x;
	  }

	  /** Return least element without removing it, or null if empty **/
	  public synchronized Object peek() {
	    if (count > 0) return nodes[0]; else return null;
	  }

	  /** Return number of elements **/
	  public synchronized int size() { return count; }
	  
	  /** remove all elements **/
	  public synchronized void clear() {
	    for (int i = 0; i < count; ++i) nodes[i] = null;
	    count = 0;
	  }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	}

}

