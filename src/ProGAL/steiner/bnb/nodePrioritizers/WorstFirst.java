package ProGAL.steiner.bnb.nodePrioritizers;

import java.util.Comparator;

import ProGAL.steiner.bnb.Node;

public class WorstFirst implements Comparator<Node>{

	@Override
	public int compare(Node n1, Node n2) {
		return -Double.compare(n1.lowerBound, n2.lowerBound);
	}

}
