//UnionFind.java, created Wed Feb 23 13:06:57 2000 by cananian
//Copyright (C) 1999 C. Scott Ananian <cananian@alumni.princeton.edu>
//Licensed under the Modified BSD Licence; see COPYING for details.
//Note: Taken from jutil, with the graceful premission of Scott
//C. Ananian.  (in jutil and earlier jpaul versions, it was called
//DisjointSet). Modified heavily by Alex Salcianu.
package ProGAL.dataStructures;
//package jpaul.DataStructs;

import java.io.Serializable;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.LinkedHashSet;

/**
 * <code>UnionFind</code> is a datastructure for performing
 * unification and lookup operations.  It uses the path compression
 * and union-by-rank heuristics to achieve O(<code>m * alpha(m,
 * n)</code>) runtime, where <code>m</code> is the total number of
 * operations, <code>n</code> is the total number of elements in the
 * set, and <code>alpha</code> denotes the *extremely* slowly-growing
 * inverse Ackermann function.
 *
 * The abstract state of the data structure at each moment is
 * determined by the previously executed unifications ({@link
 * #union}).  Each element of the universe of discourse is part of
 * exactly one equivalence class.  Initially, each elements sits in
 * its own equivalence class; each {@link #union} operation merges the
 * equivalence classes for its arguments.  Each look-up operation
 * ({@link #find}) finds the representative of the equivalence class
 * of its argument, i.e., one of the elements from that equivalence
 * class.
 * 
 * @author  C. Scott Ananian - cananian@alumni.princeton.edu
 * @author  Alexandru Salcianu - salcianu@alum.mit.edu
 */
public class UnionFind<E>  implements Serializable {
	private static final long serialVersionUID = -9141049716473285037L;
	private final Map<E,Node<E>> elmap = new LinkedHashMap<E,Node<E>>();

	/** Creates a <code>UnionFind</code>. */
	public UnionFind() { /*nothing to do*/ }


	/** Unifies the elements <code>e1</code> and <code>e2</code> and
	returns the representative of the resulting equivalence class. */
	public E union(E e1, E e2) {
		Node<E> node1 = _get_or_create_root_set(e1);
		Node<E> node2 = _get_or_create_root_set(e2);
		return _link(node1, node2).element;
	}

	/** Returns the representative of the equivalence class of
     <code>e</code>. */
	public E find(E e) {
		Node<E> node = elmap.get(e);
		if(node == null) return e;
		return _find_root(node).element;
	}


	/** Checks whether the elements <code>e1</code> and
     <code>e2</code> are unified in this union-find structure. */
	public boolean areUnified(E e1, E e2) {
		return find(e1).equals(find(e2));
	}

	/** Returns true if the element <code>e</code> has not been
     unified yet with any DIFFERENT element.  Complexity: O(1).  */
	public boolean unUnified(E e) {
		Node<E> node = elmap.get(e);
		if(node == null) {
			// never been the argument of a call to union -> true
			return true;
		}
		// even if e was sent to union in the past, maybe it was
		// unified only with equal elements; the following tests
		// check this quickly
		if(node.parent != node) {
			// already have a parent, so cannot be unUnified -> false
			return false;
		}
		// node is a root! last chance is that it's only a root of itself!
		return node.nbElems == 1;
	}


	/** If e is already in a non-trivial equivalence class (i.e., with
     more than 2 elements), then return the Node<E> corresponding
     to the representative element.  Otherwise, i.e., if e sits in
     an equivalence class by himself, then create a Node<E>(e), put
     it into elmap and return it. */
	private Node<E> _get_or_create_root_set(E e) {
		Node<E> node = elmap.get(e);
		if(node == null) {
			// no node yet; create one now
			node = new Node<E>(e);
			elmap.put(e, node);
			return node;
		}
		return _find_root(node);
	}

	/** Given a Node<E> node, walk the parent field as far as
     possible, until reaching the root, i.e., the Node<E>
     corresponding to the current representative of this
     equivalence class.  To achieve the low complexity, we also
     compress the path, by making each node a direct son of the
     root. */
	private Node<E> _find_root(Node<E> node) {
		if (node.parent != node) {
			node.parent = _find_root(node.parent);
		}
		return node.parent;
	}


	/** Unifies the tree rooted in node1 with the tree rooted in
     node2, by making one of them the subtree of the other one;
     returns the root of the new tree. */
	private Node<E> _link(Node<E> node1, Node<E> node2) {
		// maybe there is no real unification here
		if(node1 == node2) return node2;
		// from now on, we know that we unify really disjoint trees
		if(node1.rank > node2.rank) {
			node2.parent = node1;
			node1.nbElems += node2.nbElems; 
			return node1;
		}
		else {
			node1.parent = node2;
			if (node1.rank == node2.rank) {
				node2.rank++;
			}
			node2.nbElems += node1.nbElems;
			return node2;
		}
	}


	/** Returns an unmodifiable collection containing all equivalence
     classes with more than one element.  Each equivalence class is
     represented as an unmodifiable set that contains all elements
     from that class.  The reason we return only equivalence
     classes with more than one element is that this structure does
     not know about the elements of your universe of discourse that
     have not been unified with anyone yet: it knows only about the
     elements that have been unified already.

	<p> This method is quite slow: its execution time is at least
	linear in the size of the underlying forest.  Hence, it has a
	long, inconvenient name!

	<p><strong>Note:</strong> The returned collection is in fact a
	set, but maintaining it as a (hash-)set, with the associated
	.equals on sets would be too expensive.  */
	// public Collection<Set<E>> allNonTrivialEquivalenceClasses() {
	//	Relation<E,E> relRepr2Elems = new MapSetRelation<E,E>();
	//	for (E elem : elmap.keySet()) {
	//	    E repr = find(elem);
	//	    relRepr2Elems.add(repr, elem);
	//	}
	//	List<Set<E>> result = new LinkedList<Set<E>>();
	//	for(E repr : relRepr2Elems.keys()) {
	//	    Set<E> equivClass = relRepr2Elems.getValues(repr);
	//	    // skip the trivial, 1-element equivalence classes
	//	    if(equivClass.size() > 1) {
	//		result.add(relRepr2Elems.getValues(repr));
	//	    }
	//	}
	//	return Collections.unmodifiableList(result);
	// }


	/** Returns the equivalence class that element <code>e</code> is
     part of, as an unmodifiable set containing all elements from
     that class (including <code>e</code>). 

	<p> This method is quite slow: its execution time is at least
	linear in the size of the underlying forest! */
	public Set<E> equivalenceClass(E e) {
		// Case 1. an element by itself
		if(unUnified(e)) {
			return Collections.singleton(e);
		}
		// Case 2. an element from a larger equivalence class
		Set<E> equivClass = new LinkedHashSet<E>();
		// as e is not unUnified, elmap.get(e) != null
		Node<E> reprRoot = _find_root(elmap.get(e));
		for(E elem : elmap.keySet()) {
			Node<E> elemRoot = _find_root(elmap.get(elem));
			if(elemRoot == reprRoot) {
				equivClass.add(elem);
			}
		}
		return Collections.unmodifiableSet(equivClass);
	}


	/** Returns an unmodifiable set containing all elements that this
     <code>UnionFind</code> structure has knowledge about.  It is
     important to understand that this structure knows only the
     elements that it was asked to unify via calls to {@link #union
     union}; it doesn't know the other elements of the universe of
     discourse. */
	public Set<E> allKnownElements() {
		return Collections.unmodifiableSet(elmap.keySet());
	}


	// /** Returns a human-readable representation of the UnionFind. */
	// public String toString() {
	//	Relation<E,E> rel = new MapSetRelation<E,E>();
	//	for (Node<E> node : elmap.values()) {
	//	    Node<E> root = _find_root(node);
	//	    if (node != root) rel.add(root.element, node.element);
	//	}
	//	return rel.toString();
	// }

	// node representation.
	private static class Node<E> {
		Node<E> parent;
		final E element;
		// the depth of the subtree rooted in this node
		int rank = 0;
		// if this node is the root of a tree, this is number of
		// elements from the tree.  It's meaningless otherwise.
		int nbElems = 1;
		Node(E element) {
			this.parent = this;
			this.element = element;
		}
	}

}
