package ProGAL.steiner.bnb.branchers;

import java.util.Collection;
import java.util.LinkedList;

import ProGAL.steiner.bnb.Brancher;
import ProGAL.steiner.bnb.Node;

public class FollowSiteOrder implements Brancher{
	private final int N;
	
	public FollowSiteOrder(int N){
		this.N = N;
	}
	
	@Override
	public Collection<Node> branch(Node n) {
		Collection<Node> children = new LinkedList<Node>();
//		System.out.println("FollowSiteOrder.branch("+n+") .. "+n.siteInserted+" .. "+(N-1));
		if(n.siteInserted==N-1) return children;
		
		for(int e=0; e<2*n.siteInserted-1; e++){
			children.add(new Node(n,e,n.siteInserted+1));
		}
		
		return children;
	}

}
