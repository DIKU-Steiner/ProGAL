package ProGAL.proteins.beltaStructure.bnb;

import java.awt.BorderLayout;
import java.awt.Color;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeSelectionModel;

import ProGAL.dataStructures.viewer.InteractiveBinaryTree;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.PrimaryStructure;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SheetTopology;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.Atom;
import ProGAL.proteins.structure.generators.CABAminoAcidGenerator;

public class BnBTreeExplorer {
	private final PrimaryStructure primaryStructure;
	private final SecondaryStructure secondaryStructure;
	private final BetaTopology betaTopology;
	private final AminoAcidChain chain;
	private final Branchable[] parts;
	private final JFrame frame;
	private final J3DScene scene;
	
	public BnBTreeExplorer(BetaTopology bt){
		this.betaTopology = bt;
		this.secondaryStructure = bt.secondaryStructure;
		this.primaryStructure = secondaryStructure.primaryStructure;
		this.chain = new AminoAcidChain(primaryStructure.sequence, new CABAminoAcidGenerator());
	
		int loopSegments = 0;
		for(SSSegment seg: secondaryStructure.segments){
			if( seg.type!=SSType.STRAND ) loopSegments++;
		}
		List<SheetTopology> sheets = bt.getSheets();
		parts = new Branchable[sheets.size()+loopSegments];
		int i;
		for(i=0;i<sheets.size();i++)
			parts[i] = new SheetStructure(4, bt.getSheets().get(i), chain);
		for(SSSegment seg: secondaryStructure.segments){
			if( seg.type!=SSType.STRAND ) parts[i++] = new SegmentStructure(12, 4, secondaryStructure, seg, chain);
		}
		
		this.scene = J3DScene.createJ3DSceneInFrame();
		this.frame = scene.frame;
		this.scene.setAntialiasing(true);
		
		createTree();
		frame.validate();
	}
	
	
	private void createTree(){
		DefaultMutableTreeNode root = new DefaultMutableTreeNode("BnB Tree");
		
		
		DefaultMutableTreeNode n = root;
		MyNode bnRoot = new MyNode(null,-1);
		
		MyNode bn = new MyNode(bnRoot, 0);
		root.add(n = new DefaultMutableTreeNode(bn));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 0)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 1)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 2)));
		n.add(n = new DefaultMutableTreeNode(bn = new MyNode(bn, 3)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 0)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 1)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 2)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 3)));
		
		bn = new MyNode(bnRoot,1);
		root.add(n = new DefaultMutableTreeNode(bn));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 0)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 1)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 2)));
		n.add(n = new DefaultMutableTreeNode(bn = new MyNode(bn, 3)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 0)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 1)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 2)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 3)));

		bn = new MyNode(bnRoot,3);
		root.add(n = new DefaultMutableTreeNode(bn));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 0)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 1)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 2)));
		n.add(n = new DefaultMutableTreeNode(bn = new MyNode(bn, 3)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 0)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 1)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 2)));
		n.add(new DefaultMutableTreeNode(new MyNode(bn, 3)));
		

		JTree tree = new JTree(root);
		tree.getSelectionModel().setSelectionMode
        (TreeSelectionModel.SINGLE_TREE_SELECTION);
		tree.addTreeSelectionListener(new TreeSelectionListener() {
			public void valueChanged(TreeSelectionEvent e) {
				JTree t = (JTree)e.getSource();
				
				MyNode n = (MyNode)((DefaultMutableTreeNode)t.getSelectionPath().getLastPathComponent()).getUserObject();
				displayNode(n);
			}
		});
		
		JScrollPane treeView = new JScrollPane(tree);
		frame.getContentPane().add(treeView, BorderLayout.WEST);
	}
	
	private List<List<Shape>> partShapes = new LinkedList<List<Shape>>();
	private int currentPart = -1;
	void displayNode(MyNode n){
		System.out.printf("displayNode(%s) .. current = %d\n",n,currentPart);
		while(currentPart>n.part){
			for(Shape s: partShapes.get(currentPart)) scene.removeShape(s);
			partShapes.remove(currentPart);
			currentPart--;
		}
		while(currentPart<n.part){
			currentPart++;
			List<Shape> shapes = generateShapes(parts[currentPart]); 
			partShapes.add(shapes);
			for(Shape s: shapes) scene.addShape(s, Color.GRAY, 5);
		}
		while(!(n.parent==null)){
			parts[n.part].setStructure(n.structure);
			n=(MyNode)n.parent;
		}
//		for(Atom a: chain.atoms()){
//			System.out.println(a);
//		}
		scene.repaint();
	}
	
	List<Shape> generateShapes(Branchable b){
		List<Shape> shapes = new LinkedList<Shape>();
		if(b instanceof SegmentStructure){
			SegmentStructure seg = (SegmentStructure)b;
			for(int r=seg.seg.start+1;r<=seg.seg.end;r++){
				shapes.add(new Sphere(chain.atom(r, 0), 0.5));
				shapes.add(new Sphere(chain.atom(r, 1), 0.5));
				shapes.add(new LSS(chain.atom(r,0),chain.atom(r,1),0.1));
				if(r>seg.seg.start+1){
					shapes.add(new LSS(chain.atom(r-1,0),chain.atom(r,0),0.2));
				}
			}
		}else if(b instanceof SheetStructure){
			SheetStructure sstruc = (SheetStructure)b;
			for(int st: sstruc.sheetTop.strands){
				SSSegment strand =sstruc.sheetTop.secondaryStructure.getStrands()[st];
				for(int r=strand.start;r<=strand.end;r++){
					shapes.add(new Sphere(chain.atom(r, 0), 0.5));
					shapes.add(new Sphere(chain.atom(r, 1), 0.5));
					shapes.add(new LSS(chain.atom(r,0),chain.atom(r,1),0.1));
					if(r>strand.start)
						shapes.add(new LSS(chain.atom(r-1,0),chain.atom(r,0),0.2));
				}	
			}
		}
		return shapes;
	}
	
	class MyNode extends BnBNode{

		protected MyNode(BnBNode parent, int structure) {
			super(parent, structure);
		}
		
		public String toString(){
			return String.format("Node[part=%d, struc=%d]",super.part, super.structure);
		}
	}
	
	public static void main(String[] args){
		PrimaryStructure ps = new PrimaryStructure(        "AAEEKTEFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK");//1CTF
		SecondaryStructure ss = new SecondaryStructure(ps, "       EEEEEEE    HHHHHHHHHHHH   HHHHHHHHHT SEEEEEEE HHHHHHHHHHHHHHT EEEE ");
		BetaTopology bt = new BetaTopology(ss,new boolean[][]{
				new boolean[]{false, true,false},
				new boolean[]{false,false,false},
				new boolean[]{false, true,false}
				});
		BnBTreeExplorer bte = new BnBTreeExplorer(bt);
	}
}
