package ProGAL.geom3d.complex.alphaComplex;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Simplex;
import ProGAL.math.Randomization;

public class VoidTree {
	private AlphaFiltration alphaFil;
	public Node root = null;
//	J3DScene scene = J3DScene.createJ3DSceneInFrame();
	
	public VoidTree(List<Point> points, double interval) {
		this.alphaFil = new AlphaFiltration(points);
		createTree(interval);
	}
	
	private ArrayList<LinkedList<CTetrahedron>> getVoid(CTriangle m, LinkedList<CTriangle> tris){
		ArrayList<LinkedList<CTetrahedron>> v = new ArrayList<LinkedList<CTetrahedron>>();
		LinkedList<CTetrahedron> vOne = new LinkedList<CTetrahedron>();
		LinkedList<CTetrahedron> vTwo = new LinkedList<CTetrahedron>();
		v.add(vOne);
		v.add(vTwo);
		//List<CTriangle> tris = alphaFil.getTriangles(alpha);
		/*for (CTriangle tri: tris){
			scene.addShape(tri, java.awt.Color.YELLOW);
		}*/
		//List<CTetrahedron> tetra = alphaFil.getTetrahedra(); Might not be needed after all
		beginning:
			for (int i=0; i<2; i++){
				LinkedList<CTetrahedron> tempVoid = v.get(i);

				LinkedList<CTriangle> tempTris = new LinkedList<CTriangle>();
				CTetrahedron inTet = m.getAdjacentTetrahedron(i);
				if (inTet.containsBigPoint()){
					continue beginning;
				}
				tempVoid.addFirst(inTet);
				for (int j=0; j<4; j++){
					tempTris.addFirst(inTet.getTriangle(j));
				}
				tempTris.remove(m);
				while (tempTris.size()>0){
					CTriangle tri = tempTris.removeFirst();
					if (tris.contains(tri)){
						continue;
					}
					for (int k=0; k<2; k++){
						CTetrahedron adTetra = tri.getAdjacentTetrahedron(k);
						//scene.addShape(adTetra, new java.awt.Color(200,0,0,100));
						if (tempVoid.contains(adTetra)){
							continue;
						}
						if (adTetra.containsBigPoint()){
							tempVoid.clear();
							continue beginning;
						}
						tempVoid.addFirst(adTetra);
						for (int l=0; l<4; l++){
							tempTris.addFirst(adTetra.getTriangle(l));
						}
						tempTris.remove(tri);
					}
				}
			}
		return v;
	}
	
	public Node setRoot(){
		if (root == null){
			LinkedList<CTetrahedron> orgTetra = new LinkedList<CTetrahedron>(alphaFil.getTetrahedra());
		//LinkedList<CTetrahedron> emptyList = new LinkedList<CTetrahedron>();
			root = new Node(0, orgTetra);
		}
		return root;
	}
	
	public LinkedList<Node> getLeaves(Node n){
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
	
	public Node getRest(){
		Node node = root;
		while (node.left != null){
			node = node.left;
		}
		return node;
	}
	
	public Node find(LinkedList<CTetrahedron> list){
		
		LinkedList<Node> nodes = getLeaves(root);
		while (!nodes.isEmpty()){
			Node node = nodes.removeFirst();
			if (node.getTetra().contains(list.getFirst())){
				return node;
			}
		}
		return null;
	}
	
	public void createTree(double interval){
		root = setRoot();
		LinkedList<CTriangle> tris = new LinkedList<CTriangle>();
		int[][] table = alphaFil.getBettiNumbers();
		List<Simplex> simplices = alphaFil.getSimplices();
		for (int i=0; i<table[0].length;i++){
			if(table[5][i]==2){
//				scene.addShape(simplices.get(i), new java.awt.Color(100,200,100,255));
				tris.add((CTriangle)simplices.get(i));
			}
			/*if (table[5][i]==3){
				LinkedList<CTetrahedron> tetraL = new LinkedList<CTetrahedron>();
				double alphaTetra = alphaFil.getInAlpha(simplices.get(i));
				CTetrahedron tetra = (CTetrahedron)simplices.get(i);
				tetraL.add(tetra);
				Node node = find(tetraL);
				if (!(alphaTetra-node.getAlpha()>=interval)){
					
				}
			}*/
			if(table[5][i]==0) 
//				scene.addShape(
//						new Sphere( (Point)simplices.get(i),0.1 ), 
//						java.awt.Color.BLACK );
			
			if (table[5][i]==3){
				LinkedList<CTetrahedron> tets = new LinkedList<CTetrahedron>();
				CTetrahedron tetra = (CTetrahedron) simplices.get(i);
				double d = alphaFil.getInAlpha(simplices.get(i));
				tets.add(tetra);
				Node m = find(tets);
				m.changeDeath(d);
			}
			
			if (table[4][i]==1 && table[5][i]==2){
//				scene.addShape(simplices.get(i), new java.awt.Color(0,200,0,255));
				double alpha = alphaFil.getInAlpha(simplices.get(i));
				CTriangle marked = (CTriangle) simplices.get(i);
				ArrayList<LinkedList<CTetrahedron>> newVoids = getVoid(marked, tris);
				if (newVoids.get(0).isEmpty() && newVoids.get(1).isEmpty()){
					return; //Error!
				}
				if (newVoids.get(0).isEmpty()){
					Node rest = getRest();
					Node newNode = new Node(alpha, newVoids.get(1));
					rest.getTetra().removeAll(newNode.getTetra());
					Node newRest = new Node(0, rest.getTetra());
					rest.setChild(newRest, 0);
					rest.setChild(newNode, 1);
				} else {
					if (newVoids.get(1).isEmpty()){
						Node rest = getRest();
						Node newNode = new Node(alpha, newVoids.get(0));
						rest.getTetra().removeAll(newNode.getTetra());
						Node newRest = new Node(0, rest.getTetra());
						rest.setChild(newRest, 0);
						rest.setChild(newNode, 1);
					} else {					
						LinkedList<CTetrahedron> v1 = newVoids.get(0);
						LinkedList<CTetrahedron> v2 = newVoids.get(1);
						Node match = find(v1);
						Node newNode1 = new Node(alpha, v1);
						Node newNode2 = new Node(alpha, v2);
						match.setChild(newNode1, 0);
						match.setChild(newNode2, 1);
					}
				}
			}
		}
	}
	
	public static void main(String[] args){
		Randomization.seed(0);
//		ArrayList<Point> points = new ArrayList<Point>();
//		//Kube: 
//		points.add(new Point(0,1,0));
//		points.add(new Point(1,1,0));
//		points.add(new Point(1,2,0));
//		points.add(new Point(0,2,0));
//		points.add(new Point(0,1,1));
//		points.add(new Point(1,1,1));
//		points.add(new Point(1,2,1));
//		points.add(new Point(0,2,1));
		
//		points.add(new Point());
		/*tetrahedra:
		points.add(new Point(1,3,0));
		points.add(new Point(1,6,1));
		points.add(new Point(1,0,1));
		points.add(new Point(2,3,2));
		points.add(new Point(0,3,2));*/
		List<Point> points = ProGAL.geom3d.PointList.generatePointsInCube(300);
		List<Point> points2 = new ArrayList<Point>();
		for(Point p: points){
			if(p.x()>0.9 || p.x()<-0.9) points2.add(p);
			else if(p.y()>0.9 || p.y()<-0.9) points2.add(p);
			else if(p.z()>0.9 || p.z()<-0.9) points2.add(p);
		}
//		List<Point> points = ProGAL.geom3d.PointList.generatePointsOnSphere(30);
		
		System.out.println(points2.size());
		VoidTree vt = new VoidTree(points2, 0);
		new ProGAL.datastructures.viewer.BinaryTreePainter(vt.root);
	}
}
