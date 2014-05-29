package ProGAL.steiner.bnb.nodePrioritizers;

import java.util.Comparator;

import ProGAL.steiner.bnb.Node;

public class DepthFirst implements Comparator<Node>{

	@Override
	public int compare(Node n1, Node n2) {
		return n2.siteInserted-n1.siteInserted;
	}

}
