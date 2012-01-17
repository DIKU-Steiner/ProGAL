package ProGAL.dataStructures;


public class SetIterator implements java.util.Iterator<Object> {

	private Set set;
	private int current;
	
	public SetIterator(Set set) {
		this.set = set;
		current = 0;
	}
	
	public boolean hasNext() { 
		return current != set.getSize(); 
	}
	
	public Object next() {
		if (hasNext()) return set.get(current++); else return null; 
	}
	
	public void remove() {}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

