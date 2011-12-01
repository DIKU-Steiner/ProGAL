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
	private Node root = null;
	J3DScene scene = J3DScene.createJ3DSceneInFrame();
	
	public VoidTree(List<Point> points) {
		this.alphaFil = new AlphaFiltration(points);
		createTree();
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
				System.out.println("v: "+v.get(i));
				//System.out.print("inTetra");
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
				//while (v.size() > 0){
				while (tempTris.size()>0){
					System.out.println(tempTris.size());
					CTriangle tri = tempTris.removeFirst();
					if (tris.contains(tri)){
						System.out.println("Tri. is in complex");
						continue;
					}
					for (int k=0; k<2; k++){
						//System.out.println("before fail");
						CTetrahedron adTetra = tri.getAdjacentTetrahedron(k);
						//System.out.print("after fail");
						//scene.addShape(adTetra, new java.awt.Color(200,0,0,100));
						if (tempVoid.contains(adTetra)){
							continue;
						}
						if (adTetra.containsBigPoint()){
							System.out.print("Indeholder BigPoints");
							tempVoid.clear();
							System.out.println("Go-go Beginning!");
							continue beginning;
						}
						tempVoid.addFirst(adTetra);
						for (int l=0; l<4; l++){
							System.out.println("Adding triangle");
							tempTris.addFirst(adTetra.getTriangle(l));
						}
						/*Collection<CTriangle> rem = new ArrayList<CTriangle>();
						rem.add(tri);
						tempTris.removeAll(rem);*/
						tempTris.remove(tri);
						//System.out.println(tempTris.size());
					}
				}
				//}
			System.out.println("v1: "+v.get(0)+"\nv2: "+v.get(1));
			}
		return v;
	}
	
	public Node setRoot(){
		//Empty root (rest)
		LinkedList<CTetrahedron> orgTetra = new LinkedList<CTetrahedron>(alphaFil.getTetrahedra());
		//LinkedList<CTetrahedron> emptyList = new LinkedList<CTetrahedron>();
		root = new Node(0, orgTetra);
		return root;
	}
	
	public LinkedList<Node> getLeafs(Node n){
		LinkedList<Node> nodes = new LinkedList<Node>();
		/*LinkedList<Node> newNodes1 = new LinkedList<Node>();
		LinkedList<Node> newNodes2 = new LinkedList<Node>();*/
		if (n == null){
			return null;
		}
		if (n.left == null){
			nodes.add(n);
			return nodes;
		}
		else { 
			nodes.addAll(getLeafs(n.right));
			nodes.addAll(getLeafs(n.left));
			
			/*newNodes1 = getLeafs(n.left);
			newNodes2 = getLeafs(n.right);
			while (!newNodes1.isEmpty()){
				nodes.add(newNodes1.removeFirst());
			}
			while (!newNodes2.isEmpty()){
				nodes.add(newNodes2.removeFirst());
			}*/
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
		
		LinkedList<Node> nodes = getLeafs(root);
		while (!nodes.isEmpty()){
			Node node = nodes.removeFirst();
			if (node.getTetra().contains(list.getFirst())){
				return node;
			}
		}
		return null; //TODO: return "rest"-node
	}
	
	public void createTree(){
		root = setRoot();
		LinkedList<CTriangle> tris = new LinkedList<CTriangle>();
		int[][] table = alphaFil.getBettiNumbers();
		List<Simplex> simplices = alphaFil.getSimplices();
		//int first = 0;
		for (int i=0; i<table[0].length;i++){
			/*if (first==1){
				break;
			}*/
			//System.out.println("dimension: "+table[5][i]);
			if(table[5][i]==2){
				scene.addShape(simplices.get(i), new java.awt.Color(100,200,100,255));
				tris.add((CTriangle)simplices.get(i));
			}
			if(table[5][i]==0) 
				scene.addShape(
						new Sphere( (Point)simplices.get(i),0.1 ), 
						java.awt.Color.BLACK );
			
			if (table[4][i]==1 && table[5][i]==2){
				//first=1;
				//System.out.println("index: "+i);
				scene.addShape(simplices.get(i), new java.awt.Color(0,200,0,255));
				double alpha = alphaFil.getInAlpha(simplices.get(i));
				CTriangle marked = (CTriangle) simplices.get(i);
				//System.out.print("Kald til voids.");
				ArrayList<LinkedList<CTetrahedron>> newVoids = getVoid(marked, tris);
				System.out.print(newVoids.get(0).size()+"-");
				System.out.print(newVoids.get(1).size());
				if (newVoids.get(0).isEmpty() && newVoids.get(1).isEmpty()){
					System.out.print("Begge voids er tomme.");
					return; //Error!
				}
				if (newVoids.get(0).isEmpty()){
					System.out.println("0 Void er tomt.");
					Node rest = getRest();
					Node newNode = new Node(alpha, newVoids.get(1));
					rest.getTetra().removeAll(newNode.getTetra());
					Node newRest = new Node(0, rest.getTetra());
					rest.setChild(newRest, 0);
					rest.setChild(newNode, 1);
				} else {
					if (newVoids.get(1).isEmpty()){
						System.out.println("1 Void er tomt.");
						Node rest = getRest();
						System.out.print(rest==root);
						Node newNode = new Node(alpha, newVoids.get(0));
						rest.getTetra().removeAll(newNode.getTetra());
						Node newRest = new Node(0, rest.getTetra());
						rest.setChild(newRest, 0);
						rest.setChild(newNode, 1);
					} else {
						System.out.print("Begge voids findes.");					
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
		ArrayList<Point> points = new ArrayList<Point>();
		//Kube: 
		points.add(new Point(0,1,0));
		points.add(new Point(1,1,0));
		points.add(new Point(1,2,0));
		points.add(new Point(0,2,0));
		points.add(new Point(0,1,1));
		points.add(new Point(1,1,1));
		points.add(new Point(1,2,1));
		points.add(new Point(0,2,1));
		/*tetrahedra:
		points.add(new Point(1,3,0));
		points.add(new Point(1,6,1));
		points.add(new Point(1,0,1));
		points.add(new Point(2,3,2));
		points.add(new Point(0,3,2));*/
		//List<Point> points = ProGAL.geom3d.PointList.generatePointsOnSphere(4);
		VoidTree vt = new VoidTree(points);
		new ProGAL.datastructures.viewer.BinaryTreePainter(vt.root);
	}
}
