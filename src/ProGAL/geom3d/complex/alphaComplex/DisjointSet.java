package ProGAL.geom3d.complex.alphaComplex;

import java.util.LinkedList;
import java.util.List;

/** 
 * A class for maintaining disjoint sets using union-by-rank and path-compression. This gives an amortized 
 * runningtime of each union and find-operation of O(\alpha(n)). Except that a set is represented as an internal 
 * class that is associated with objects using a hash-map, this implementation follows the description in 
 * Cormen et. al 2001 and on http://en.wikipedia.org/wiki/Disjoint-set_data_structure (Jan. 2011).

 * @author R. Fonseca
 */
public class DisjointSet {
	private List<Set> sets = new LinkedList<Set>();
	
	/** TODO: Improve so this is not O(n) */
	private Set getSet(Object o){
		for(Set s: sets) if(s.representative==o) return s;
		return null;
	}
	

	/** Merge two sets into a single set. Warning: To associate objects with sets the set has to be looked up in 
	 * a list which takes O(n) time. */
	public void union(Object x, Object y){
		union(getSet(x),getSet(y));
	}
	
	/** Merge two sets into a single set */
	public void union(Set x, Set y){
		Set xRoot = find(x);
		Set yRoot = find(y);
		if(xRoot.rank>yRoot.rank){
			yRoot.parent = xRoot;
		}else if(xRoot!=yRoot){
			xRoot.parent = yRoot;
			if(xRoot.rank==yRoot.rank) 
				yRoot.rank++;
		}
	}
	
	/** 
	 * Determine which set a particular element is in. If one wishes to determine if two objects belong to 
	 * the same set, the references returned by <code>find</code> can be compared by equality (i.e. "==" ). 
	 */
	public Set find(Set s){
		if(s.parent == s)
			return s;
		else{
			s.parent = find(s.parent);
			return s.parent;
		}
	}
	
	/** 
	 * Determine which set a particular element is in. If one wishes to determine if two objects belong to 
	 * the same set, the references returned by <code>find</code> can be compared by equality (i.e. "==" ). 
	 * Warning: To associate objects with sets the set has to be looked up in a list which takes O(n) time.
	 */
	public Set find(Object o){
		Set s = getSet(o);
		if(s==null) return null;
		return find(s);
	}
	
	
	/** Create a new set containing the element o */
	public Set makeSet(Object o){
		Set s = getSet(o);
		if(s==null) {
			s = new Set();
			s.parent = s;
			s.representative = o;
			sets.add(s);
		}
		return s;
	}
	
	/** A class representing a set. Internally this class is used to store a parent and a rank. A reference to 
	 * a Set can also be compared to other Set-references. */
	public class Set{
		private Set parent;
		private int rank = 0;
		private Object representative;
	}
}
