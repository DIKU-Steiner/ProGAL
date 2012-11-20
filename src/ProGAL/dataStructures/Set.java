package ProGAL.dataStructures;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Random;

import ProGAL.geom2d.Point;

public class Set<T> implements Iterable<T>{
	protected Object[] elements;
	protected int n;
	
	/*
	 * creates an empty set with capacity 0
	 */
	public Set() { this(0); }
	
	/*
	 * creates an empty set with the specified capacity
	 */
	public Set(int capacity)  {
		if (capacity < 0) throw new IllegalArgumentException();
		elements = new Object[capacity];
		n = 0;
	}
	
	/*
	 * creates a set containing the specified element and with capacity 1
	 */
	public Set(Object element) {
		n = 1;
		elements = new Object[n];
		elements[0] = element;
	}
	
	/*
	 * creates a set with the specified elements 
	 * and with the capacity equal to the number of elements.
	 */
	public Set(Object[] elements) {
		n = elements.length;
		this.elements = new Object[n];
		System.arraycopy(elements, 0, this.elements, 0, n);
	}
	
	/*
	 * creates a set that is a copy of the specified set
	 */
	public Set(Set<T> set) {
		n = set.getSize();
		elements = new Object[n];
		for (int i = 0; i < n; i++) elements[i] = set.get(i);
	}
	
	/*
	 * makes this set empty
	 */ 
	public void clear() {
		for (int i = 0; i < n; i++) elements[i] = null;
		n = 0;
	}
	
	public boolean isEmpty() { return n == 0; }
	
	@SuppressWarnings("unchecked")
	public T get(int i) {
		if ((i < 0) || (i >= n)) throw new IllegalArgumentException("index out of range");
		else return (T)elements[i];
	}
	
	@SuppressWarnings("unchecked")
	public T getFirst() {
		if (n == 0) throw new IllegalArgumentException("the set is empty");
		else return (T)elements[0];

	}
	
	@SuppressWarnings("unchecked")
	public T getLast() {
		if (n == 0) throw new IllegalArgumentException("the set is empty");
		else return (T)elements[n-1];
	}
	
	public void set(int i, Object element) { elements[i] = element; }
	
	public Object[] getElements() { return elements; }
	
	public int getSize() { return n; }
	
	public void randomPermutation() {
		Random random = new Random();
		int j;
		for (int i = 0; i < n; i++) {
			j = random.nextInt(n-i);
			swap(i,i+j);
		}
	}
	
	/*
	 * returns index of specified object
	 */
	public int findIndex(Object object) {
		for (int k = 0; k < n; k++) { if (elements[k] == object) return k; }
		return -1;
	}
	
	/*
	 * returns TRUE if the specified object is in the set.
	 */
	public boolean contains(Object object) { return findIndex(object) != -1; }
		
	public boolean isMember(Object object) { return findIndex(object) != -1; }
	
	/*
	 * inserts new object
	 */
	public void insert(T object) {
		if (n >= elements.length) {
			int newCapacity = 3*elements.length/2 + 1;
			Object[] newElements = new Object[newCapacity];
			System.arraycopy(elements, 0, newElements, 0, elements.length);
			elements = newElements;
		}
		elements[n++] = object;
	}
	
	public void append(Set<T> setToAppend) {
		int size = setToAppend.getSize();
		for (int i = 0; i < size; i++) insert(setToAppend.get(i));
	}
	
	public void append(T[] array){ for(T t: array) insert(t); }
	
	
	public void reverse() {
		Object[] tempElements = new Object[n];
		System.arraycopy(elements, 0, tempElements, 0,n);
		for (int i = 0; i < n; i++ ) elements[i] = tempElements[n-i-1];
	}
	
	/*
	 * deletes specified object. Note!!! Order is not preserved !!!
	 */
	public void delete(Object object) { deleteIndex(findIndex(object)); }

	/*
	 * deletes object with specified index. Note !!! Order is not preserved !!!
	 */
	public Object deleteIndex(int k) {
		if ((k < 0) || (k >= n)) throw new IllegalArgumentException("object not in the set");
		Object object = elements[k];
		elements[k] = elements[n-1];
		elements[--n] = null;
		if (n <= elements.length/4) {
			int newCapacity = elements.length/2 + 1;
			Object[] newElements = new Object[newCapacity];
			System.arraycopy(elements, 0, newElements, 0, n);
			elements = newElements;
		}
		return object;
	}

	/*
	 * deletes the first object (with lowest index) in the set.
	 */
	public Object deleteFirst() { return deleteIndex(0); }
	
	/*
	 * deletes the last object (with highest index) in the set.
	 */
	public Object deleteLast() { return deleteIndex(n-1); }
	
	/*
	 * swaps two elements of the array
	 */
	public void swap(int i, int j) {
		if ((i < 0) || (i >= n) || (j < 0) || (j >= n)) 
			throw new IllegalArgumentException("object not in the set");
		Object temp = elements[i];
		elements[i] = elements[j];
		elements[j] = temp;
	}

	/*
	 * shifts the elements so that the i-th element becomes first element.
	 */
	public void shift(int i) { 
		Object[] elementsCopy = new Object[i];
		System.arraycopy(elements,0,elementsCopy,0, i);
		for (int j = i; j < n; j++) elements[j-i] = elements[j];
		for (int j = 0; j < i; j++) elements[getSize()-i+j] = elementsCopy[j];
	}
	
	public boolean isEqual(Set<T> set2) {
		for (int i = 0; i < n; i++) if (!set2.contains(get(i))) return false;
		return true;
	}
	
	public Object binarySearch(SortTool tool, Object object) { return binarySearch(tool, object, 0, getSize()-1); }
	
	private Object binarySearch(SortTool tool, Object object, int left, int right) {
		if (right == left) return null;
		else {
			int mid = left + (right-left)/2;
			if (tool.compare(object, get(mid)) < 0) return binarySearch(tool, object, left, mid-1);
			else {
				if (tool.compare(object, get(mid)) > 0) return binarySearch(tool, object, mid+1, right);
				else return get(mid);
			}
		}
	}
	
	public int partition(SortTool tool, int left, int right) {
		int i = left;
		Object pivot = get(right);
		for (int k = left; k < right; k++) {
			if (tool.compare(get(k), pivot) < 0) swap(i++,k);
		}
		swap(i,right);
		return i;
	}
	
	public void sort() {
		if (!isEmpty()) {
			SorterQuick sorter = new SorterQuick();
			Object obj = get(0);
			if (obj instanceof Integer) sorter.Sort(this, new SortToolInteger());
			else {
				if (obj instanceof Double) sorter.Sort(this, new SortToolDouble());
				else {
					if (obj instanceof String) sorter.Sort(this, new SortToolString());
					else {
						if (obj instanceof Point) sorter.Sort(this,  new SortToolPoint2dDistance());
					}
				}
			}
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Set<Integer> s = new Set<Integer>(100);
		Random randGen = new Random();
		for (int i = 0; i < 100; i++) s.insert(new Integer(randGen.nextInt(1000)));
		for (int i = 0; i < 100; ++i) System.out.print(s.get(i) + " ");
		System.out.println();
		s.sort();
		for (int i = 0; i < 100; ++i) System.out.print(s.get(i) + " ");
		System.out.println();
	}

	public Iterator<T> iterator() {
		return new Iterator<T>(){
			int cur = 0;
			public boolean hasNext() {
				return cur<Set.this.getSize();
			}

			@SuppressWarnings("unchecked")
			public T next() {
				if(hasNext())
					return (T)elements[cur++];
				else
					throw new NoSuchElementException();
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
			
		};
	}

}
