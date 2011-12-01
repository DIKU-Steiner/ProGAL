package ProGAL.geom3d.complex.alphaComplex;

import java.util.LinkedList;

import ProGAL.geom3d.complex.CTetrahedron;

public class Node {
	public Node left = null;
	public Node right = null;
	private double alpha;
	private LinkedList<CTetrahedron> tetra;
	
	public Node(double a, LinkedList<CTetrahedron> t){
		alpha = a;
		tetra = t;
	}
	
	public double getAlpha(){
		return alpha;
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
}
