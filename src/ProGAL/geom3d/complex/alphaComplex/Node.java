package ProGAL.geom3d.complex.alphaComplex;

import java.awt.Color; 
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import java.util.Stack;

import javax.swing.JFrame;

//import ProGAL.dataStructures.viewer.InteractiveBinaryTree;
import ProGAL.dataStructures.viewer.InteractiveBinaryTree;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;

public class Node implements InteractiveBinaryTree {
	public Node left = null;
	public Node right = null;
	private double alpha;
	private double death;
	private LinkedList<CTetrahedron> tetra;
	
	public Node(double a, LinkedList<CTetrahedron> t){
		alpha = a;
		death = -1;
		tetra = t;
	}
	
	public double getAlpha(){
		return alpha;
	}
	
	public void changeAlpha(double newA){
		alpha = newA;
	}
	
	public double getDeath(){
		return death;
	}
	
	public void changeDeath(double newD){
		death = newD;
	}
	
	public LinkedList<CTetrahedron> getTetra(){
		return tetra;
	}
	
	public void setChild(Node n, int index){
		if (index == 0){
			left = n;
		} else {
			right = n;
		}
	}

	private LinkedList<Node> getLeaves(Node n){
		LinkedList<Node> nodes = new LinkedList<Node>();
		if (n == null){
			return null;
		}
		if (n.left == null){
			nodes.add(n);
			return nodes;
		}
		else { 
			nodes.addAll(getLeaves(n.right));
			nodes.addAll(getLeaves(n.left));
			return nodes;
		}
	}
	
	// Calculates persistence of a void (node) from birth time of that void and the latest death time of the leaves below the void.
	//Does not work on "rest" nodes!
	private double maxPersist(){
		LinkedList<Node> leaves = getLeaves(this);
		double b = this.alpha;
		double d = 0;
		for (int i = 0; i<leaves.size();i++) {
			Node leaf = leaves.get(i);
			if (d<leaf.death){
				d=leaf.death;
			}
		}
		return d-b;
	}
	
	
	public Node left() {return left;}
	public Node right() {return right;}
	public Color leftLegColor() {return Color.BLACK;}
	public Color rightLegColor() {return Color.BLACK;}
	public Color nodeColor() {
		if(tetra.size()==0){
			return Color.RED;
		}
		return Color.BLACK;}
	public String label() { return String.format("");
	}
	public void click() {
		System.out.println(String.format("a: %.20f, tets: %d, death: %.20f, persist: %.20f", alpha, tetra.size(), death, maxPersist()));
		
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		scene.frame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		
		//Collect all vertices using iterative DFS
		Set<Point> vertices = new HashSet<Point>();//A set ensures no duplicates
		Set<CTetrahedron> closed = new HashSet<CTetrahedron>();
		Stack<CTetrahedron> fringe = new Stack<CTetrahedron>();
		fringe.add(tetra.get(0));
		while(!fringe.isEmpty()){
			CTetrahedron t = fringe.pop();
			vertices.add(t.getCorner(0));
			vertices.add(t.getCorner(1));
			vertices.add(t.getCorner(2));
			vertices.add(t.getCorner(3));
			closed.add(t);
			CTetrahedron n0 = t.getNeighbour(0);
			CTetrahedron n1 = t.getNeighbour(1);
			CTetrahedron n2 = t.getNeighbour(2);
			CTetrahedron n3 = t.getNeighbour(3);
			if(!n0.containsBigPoint() && !closed.contains(n0)) fringe.add(n0);
			if(!n1.containsBigPoint() && !closed.contains(n1)) fringe.add(n1);
			if(!n2.containsBigPoint() && !closed.contains(n2)) fringe.add(n2);
			if(!n3.containsBigPoint() && !closed.contains(n3)) fringe.add(n3);
		}
		
		//Now paint vertices and void-tetrahedra
		for(Point p: vertices)		scene.addShape(new Sphere(p,0.3), Color.BLUE, 8);
		for(CTetrahedron t: tetra)	scene.addShape(t, Color.GREEN);
	}
}
