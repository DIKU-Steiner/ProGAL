package ProGAL.steiner.bnb;

import java.util.Arrays;

public class Topology {
	public final int N;
	public int sites;
	private final int[] kVector, sVector; 
	public int[][] edges;
	public int[][] steinerAdjacencies;
	
	private Node currentNode = null;
	
	/**
	 * Construct a topology allowing for up to maxN sites (and maxN-2 Steiner points).  
	 * @param maxN 
	 */
	public Topology(int maxN){
		assert maxN>0: "maxN can not be negative";
		N = maxN;
		edges = new int[N+N-3][2];
		steinerAdjacencies = new int[N-2][3];
		kVector = new int[N-3];
		sVector = new int[N-3];
		sites = 3;
	}
	
	/**
	 * Update the siteAdjacencies and steinerAdjacencies fields to reflect the specified node. 
	 * @param n a node in a bnb tree
	 */
	public void updateFromNode(Node n){
		assert n!=null;

		//TODO: Following commented code might be correct and would speed up the method in some cases. Test.
		//Node prevNode = currentNode;
		currentNode = n;
		sites = n.depth+3;
		
		while(n.parent!=null){ // && n!=prevNode ){
			kVector[n.depth-1] = n.edgeSplit;
			sVector[n.depth-1] = n.siteInserted;
			n = n.parent;
		}
		
		//Create edges
		int newSteiner = N;
		int[][] edge = edges; //Edges
		edge[0][0] = 0; edge[0][1] = newSteiner;
		edge[1][0] = 1; edge[1][1] = newSteiner;
		edge[2][0] = 2; edge[2][1] = newSteiner++;

		int e = 3;
	    for(int i=0; i<currentNode.depth; i++){ 
	        int s = kVector[i], v2 = edge[s][1]; //Split edge
	        edge[e][0] = sVector[i]; 	edge[e++][1] = newSteiner;
	        edge[e][0] = v2; 			edge[e++][1] = newSteiner;
	        edge[s][1] = newSteiner++;
	    }
		
	    int[] count = new int[N-2]; //Count is used to keep track of how full each separate Steiner-node adjacency entry is

	    int[][] adj = steinerAdjacencies;
		adj[0][0] = 0; adj[0][1] = 1; adj[0][2] = 2;  

		for(int i=0;i<e;i++){
	    	int s1 = edge[i][0]-N; if(s1>=0) adj[s1][count[s1]++] = edge[i][1];
	    	int s2 = edge[i][1]-N; if(s2>=0) adj[s2][count[s2]++] = edge[i][0];
	    }
	    
	}
	
	public Topology clone(){
		Topology ret = new Topology(N);
		ret.sites = sites;
		for(int i=0;i<edges.length;i++){
			ret.edges[i][0] = edges[i][0];
			ret.edges[i][1] = edges[i][1];
		}
		for(int i=0;i<steinerAdjacencies.length;i++){
			ret.steinerAdjacencies[i][0] = steinerAdjacencies[i][0];
			ret.steinerAdjacencies[i][1] = steinerAdjacencies[i][1];
			ret.steinerAdjacencies[i][2] = steinerAdjacencies[i][2];
		}
		Arrays.fill(ret.kVector, -1);
		Arrays.fill(ret.sVector, -1);
		return ret;
	}
	
	public String toString(){
		StringBuilder sb = new StringBuilder();
		sb.append("Topology:\n");
		sb.append("> k-vector: "+Arrays.toString(kVector)+"\n");
		sb.append("> edges: ");
		for(int e=0;e<sites+sites-3;e++){
			sb.append(String.format("%3d",edges[e][0]));
		}
		sb.append("\n>        ");
		for(int e=0;e<sites+sites-3;e++){
			sb.append(String.format("%3d",edges[e][1]));
		}
		sb.append("\n> steinerAdjacencies: ");
		for(int s=0;s<sites-2;s++){
			sb.append(String.format("%3d",steinerAdjacencies[s][0]));
		}
		sb.append("\n>                     ");
		for(int s=0;s<sites-2;s++){
			sb.append(String.format("%3d",steinerAdjacencies[s][1]));
		}
		sb.append("\n>                     ");
		for(int s=0;s<sites-2;s++){
			sb.append(String.format("%3d",steinerAdjacencies[s][2]));
		}
		return sb.toString();
	}
	
	public int hashCode(){
		int[] feats = new int[edges.length];
		for(int e=0;e<edges.length;e++){
			feats[e] = edges[e][0]+edges[e][1];
		}
		Arrays.sort(feats);
		return Arrays.hashCode(feats);
	}
	
	public boolean equals(Topology top){
		if(this==top) return true;
		if(this.sites!=top.sites) return false;
		int[][] adj = steinerAdjacencies, tAdj = top.steinerAdjacencies;
		for(int a=0;a<sites-2;a++){
			if( adj[a][0]!=tAdj[a][0] && adj[a][0]!=tAdj[a][1] && adj[a][0]!=tAdj[a][2] ) return false;
			if( adj[a][1]!=tAdj[a][0] && adj[a][1]!=tAdj[a][1] && adj[a][1]!=tAdj[a][2] ) return false;
			if( adj[a][2]!=tAdj[a][0] && adj[a][2]!=tAdj[a][1] && adj[a][2]!=tAdj[a][2] ) return false;
		}
		return true;
		
	}
}
