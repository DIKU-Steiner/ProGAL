package ProGAL.dataStructures.viewer;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.swing.*;

public class BinaryTreePainter extends JFrame{
	private static final long serialVersionUID = 1L;
	public PNode root;
	private PaintPane pp;
	
	public BinaryTreePainter(InteractiveBinaryTree treeRoot){
		super("Binary tree painter");
		super.setSize(600,500);
		
		setup();
		setRoot(treeRoot);
		super.setVisible(true);
	}
	
	void setup(){
		pp = new PaintPane();
		super.getContentPane().add(pp);
	}
	
	public void setRoot(InteractiveBinaryTree treeRoot){
		this.root = new PNode(treeRoot,0);
		pp.updateComponents();
	}
	
	
	private class PaintPane extends JLayeredPane{
		private static final long serialVersionUID = 1L;

		Map<PNode,JComponent> componentMap = new HashMap<PNode,JComponent>();
		Map<JComponent,PNode> nodeMap = new HashMap<JComponent,PNode>();
		PaintPane(){	
			super();
		}
		double[] treeBorder = new double[4];
		int gap = 20;
		
		private void updateComponents(){
			List<PNode> nodes = new LinkedList<PNode>();
			root.collect(nodes);
			for(PNode n: nodes) {
				JComponent comp = new JComponent(){
					private static final long serialVersionUID = 1L;
				};
				super.add(comp,2000);
				Point p = toPoint(n.x,n.y);
				comp.setBounds(p.x-5,p.y-5,10,10);
				componentMap.put(n, comp);
				nodeMap.put(comp,n);
				comp.addMouseListener(new MouseAdapter(){
					public void mousePressed(MouseEvent evt){
						InteractiveBinaryTree bt = nodeMap.get(evt.getSource()).n;
						if(bt instanceof InteractiveBinaryTree){
							((InteractiveBinaryTree)bt).click();
							repaint();
						}
					}
				});
			}
		}

		public void paint(Graphics g){
			Graphics2D g2d = (Graphics2D)g;
			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
			g2d.setColor(Color.white);
			g2d.fillRect(0, 0, getWidth(), getHeight());
			g2d.setColor(Color.BLACK);
			super.paint(g);
			setTreeBorder();
			paintNode(root,g2d);
		}
		
		private Point paintNode(PNode n, Graphics2D g){
			Point p = toPoint(n.x,n.y);
			InteractiveBinaryTree bt = n.n;
			g.setColor(Color.black);
			if(n.left!=null) {
				Point leftP = paintNode(n.left,g);
				g.setColor( bt.leftLegColor() );
				g.drawLine(p.x, p.y, leftP.x, leftP.y);
			}
			if(n.right!=null) {
				Point rightP = paintNode(n.right,g);
				g.setColor( bt.rightLegColor() );
				g.drawLine(p.x, p.y, rightP.x, rightP.y);
			}
			
//			g.setColor( bt.nodeColor() );
//			g.fillOval(p.x-3, p.y-3, 6,6);
			
			if(bt.label()!=null && bt.label().length()>0){
				g.setColor(new Color(255,255,255,200));
				g.fillRect(p.x, p.y+3, 60, 15);
				g.setColor(Color.BLACK);
				g.drawString( bt.label(), p.x, p.y+15);
			}
			g.setColor(Color.black);
			
			componentMap.get(n).setBounds(p.x-5, p.y-5, 10,10);
			
			return p;
		}
		
		void setTreeBorder(){
			List<PNode> nodes = new LinkedList<PNode>();
			root.collect(nodes);
			for(PNode n: nodes) {
				if(n.x<treeBorder[0]) treeBorder[0]=n.x;
				if(n.x>treeBorder[1]) treeBorder[1]=n.x;
				if(n.y<treeBorder[2]) treeBorder[2]=n.y;
				if(n.y>treeBorder[3]) treeBorder[3]=n.y;
			}
		}
		
		Point toPoint(double x, double y){
			double canonicalX = (x-treeBorder[0])/(treeBorder[1]-treeBorder[0]);
			double canonicalY = (y-treeBorder[2])/(treeBorder[3]-treeBorder[2]);
			return new Point((int)(canonicalX*(getWidth()-gap-gap))+gap,(int)(canonicalY*(getHeight()-gap-gap))+gap);
		}
	}

	public static class PNode{
		InteractiveBinaryTree n;
		PNode left, right;
		double x,y;
		int height, depth;
		public PNode(InteractiveBinaryTree n, int depth){
			this.n = n;
			this.depth = depth;
			if(n.left()!=null) left = new PNode(n.left(), depth+1);
			if(n.right()!=null) right = new PNode(n.right(), depth+1);
			
			if(n.left()==null){
				if(n.right()==null) height=0;
				else height = right.height+1;
			}else {
				if(n.right()==null) height = left.height+1;
				else height = Math.max(left.height, right.height)+1; 
			}
			y = depth;
			x = 0;
			
			if(left!=null && right!=null){
				right.translate(1-horizontalDistance(left, right));
				double centerDX = (left.x+right.x)/2;
				left.translate(-centerDX);
				right.translate(-centerDX);
			}else if(left!=null){
				left.translate(-left.x-1);
			}else if(right!=null){
				right.translate(-right.x+1);
			}
		}

		private void translate(double dX){
			x+=dX;
			if(left!=null) left.translate(dX);
			if(right!=null) right.translate(dX);
		}
		
		public static double horizontalDistance(PNode left, PNode right){
			List<PNode> leftNodes = new LinkedList<PNode>();
			List<PNode> rightNodes = new LinkedList<PNode>();
			left.collect(leftNodes);
			right.collect(rightNodes);
			double minDist = Double.POSITIVE_INFINITY;
			for(int d=Math.min(left.depth,right.depth);;d++){
				double leftExtremeRight = Double.NEGATIVE_INFINITY;
				double rightExtremeLeft = Double.POSITIVE_INFINITY;
				boolean foundLeft= false, foundRight = false;
				for(PNode n: leftNodes){
					if(n.depth==d){
						foundLeft=true;
						if(n.x>leftExtremeRight) leftExtremeRight = n.x;
					}
				}
				for(PNode n: rightNodes){
					if(n.depth==d){
						foundRight = true;
						if(n.x<rightExtremeLeft) rightExtremeLeft = n.x;
					}
				}
				if(!foundLeft||!foundRight) break;
				minDist = Math.min(minDist, rightExtremeLeft-leftExtremeRight);
			}
			
			return minDist;
//			if(left.depth!=right.depth) throw new Error("Ahem!");
//			double dist = right.x-left.x;
//			PNode leftChild = left.right;
//			if(leftChild==null) leftChild = left.left;
//			PNode rightChild = right.left;
//			if(rightChild==null) rightChild = right.right;
//			if(rightChild==null||leftChild==null) return dist;
//			
//			return Math.min(dist, horizontalDistance(leftChild,rightChild));
		}
		private void collect(List<PNode> collection){
			collection.add(this);
			if(left!=null) left.collect(collection);
			if(right!=null) right.collect(collection);
		}
		public String toString(){
			return "PNode["+x+","+y+"]";
		}
	}
	
	public static void main(String[] args){
		class MyNode implements InteractiveBinaryTree{
			MyNode right,left;
			Object o;
			MyNode(MyNode left, MyNode right, Object o){
				this.left = left;
				this.right = right;
				this.o = o;
			}
			MyNode(Object o){this(null,null,o); }
			public InteractiveBinaryTree left() {	return left;		}
			public InteractiveBinaryTree right() {	return right;		}
			public Color leftLegColor() {return Color.BLACK;}
			public Color rightLegColor() {return Color.BLACK;}
			public Color nodeColor() {return Color.BLACK;}
			public String label() { return o.toString();	}
			public void click() {}
		}
		
		MyNode n0 = new MyNode(0);
		MyNode n1 = new MyNode(1);
		MyNode n2 = new MyNode(2);
		MyNode n3 = new MyNode(3);
		MyNode n4 = new MyNode(4);
		MyNode n5 = new MyNode(5);
		
		MyNode n01 = new MyNode(n0,n1,"01");
		MyNode n23 = new MyNode(n2,n3,"23");
		MyNode n234 = new MyNode(n23,n4,"23,4");
		MyNode n01234 = new MyNode(n01,n234,"01,234");
		MyNode root = new MyNode(n01234,n5,"01234,5");
		
		new BinaryTreePainter(root);
	}
}
