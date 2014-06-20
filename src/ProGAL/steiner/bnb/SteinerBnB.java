package ProGAL.steiner.bnb;

import java.util.Comparator;
import java.util.PriorityQueue;

import ProGAL.math.Constants;
import ProGAL.steiner.bnb.branchers.FollowSiteOrder;
import ProGAL.steiner.bnb.lowerBounds.PartialMST;
import ProGAL.steiner.bnb.nodePrioritizers.DepthFirst;

public class SteinerBnB {
	private final LowerBound lowerBound;
	private final Brancher brancher;
	private final PriorityQueue<Node> nodePool;
	
	public int nodesVisited;
	
	public SteinerBnB(LowerBound lowerBound, Brancher brancher, Comparator<Node> nodePriority){
		this.lowerBound = lowerBound;
		this.brancher = brancher;
		this.nodePool = new PriorityQueue<Node>(10000, nodePriority);
	}
	
	public SteinerBnB(ProGAL.geomNd.Point[] sites){
		this(new PartialMST(sites, 0.0001), new FollowSiteOrder(sites.length), new DepthFirst());
	}
	
	public Topology solve(int N, double upperBound){
		upperBound+=0.00001;
		nodesVisited = 0;
		
		Node bestNode = null;
		
		Node root = new Node();
		root.lowerBound = lowerBound.lowerBound(root);
		nodePool.add(root);
		
		while(!nodePool.isEmpty()){
			Node n = nodePool.poll();
			nodesVisited++;
			
			for(Node child: brancher.branch(n)){
				child.lowerBound = lowerBound.lowerBound(child);
				if(child.lowerBound<upperBound){
					if(child.depth+3==N){
						bestNode = child;
						upperBound = child.lowerBound;
//						System.out.println("Improved upper bound: "+upperBound);
					}else{
						nodePool.add(child);
					}
				}
			}
		}
		
		Topology top = new Topology(N);
		top.updateFromNode(bestNode);
		return top;
	}
}
