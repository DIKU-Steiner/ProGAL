package ProGAL.steiner.bnb.nodePrioritizers;

import java.util.Comparator;

import ProGAL.steiner.bnb.Node;

public class BreadthFirst implements Comparator<Node>{

	@Override
	public int compare(Node n1, Node n2) {
		return n1.siteInserted-n2.siteInserted;
	}

}
