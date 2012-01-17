package ProGAL.proteins.beltaStructure.sheetLoop;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.surface.ParametricParaboloid;
import ProGAL.geom3d.viewer.ClickListener;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.PrimaryStructure;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.beltaStructure.loop.LoopStructure;
import ProGAL.proteins.beltaStructure.sheet.SurfaceSheetStructure;
import ProGAL.proteins.structure.Atom;

/** 
 * A viewer and editor of SheetLoopStructures.
 *  
 * @author R.Fonseca
 */
public class SheetLoopGUI implements ClickListener{

	private final SheetLoopStructure struc;
	private final JFrame frame;
	private final J3DScene scene;
	private final InfoPanel infoPanel;
	private final ElementListPanel listPanel; 
	private final ElementPropertiesPanel propertiesPanel;
	
	public SheetLoopGUI(SheetLoopStructure struc){
		this.struc = struc;
		this.scene = J3DScene.createJ3DSceneInFrame();
		this.frame = scene.frame;
		
		scene.addClickListener(this);
		setupScene();

		frame.getContentPane().add(infoPanel = new InfoPanel(), BorderLayout.NORTH);
		JPanel p = new JPanel(new BorderLayout());
		p.add(listPanel = new ElementListPanel(), BorderLayout.WEST);
		p.add(propertiesPanel = new ElementPropertiesPanel(), BorderLayout.CENTER);
		frame.getContentPane().add(p, BorderLayout.SOUTH);
		frame.validate();
	}
	
	private void setupScene(){
		scene.removeAllShapes();
		struc.updateAtoms(struc);
		
		Point o = new Point(0,0,0);
		for(Atom a: struc.atoms()){
			if(!a.equals(o)) {
				scene.addShape(new Sphere(a, a.radius()/2), color(a), 9);
			}
		}
	}
	
	private final static Color color(Atom a){
		switch(a.element()){
		case 'C': return Color.GRAY;
		case 'O': return Color.RED;
		case 'N': return Color.BLUE;
		case 'S': return Color.YELLOW;
		case 'H': return Color.WHITE;
		}
		
		return Color.GREEN.darker();
	}
	

	@Override
	public void shapeClicked(Shape shape, MouseEvent event) {
		if(!(shape instanceof Sphere)) return;
		Sphere clickedSphere = (Sphere)shape;
		
		//Locate the corresponding atom
		for(Atom a: struc.atoms()){
			if(a==clickedSphere.getCenter()){//When found, activate the interface
				infoPanel.atomClicked(a);
				PartialStructure ps = getPartialStructure(a);
				listPanel.elementClicked(ps);
				break;
			}
		}
	}
	private PartialStructure getPartialStructure(Atom a){
		SSSegment seg = struc.secondaryStructure.getSegmentContainingResidue(a.aminoAcid().index());
		if(seg.type==SSType.STRAND) return struc.sheetStructures.get(0);//TODO: Support more than one sheet
		for(LoopStructure ls: struc.loopStructures){
			if(seg.segmentIndex>=ls.segment1.segmentIndex && seg.segmentIndex<=ls.segment2.segmentIndex)
				return ls;
		}
		return null;
	}

	
	private class InfoPanel extends JPanel{
		private static final long serialVersionUID = 1L;
		
		JTextField textField = new JTextField();
		InfoPanel(){
			super(new BorderLayout());
			add(textField);
		}
		
		void atomClicked(Atom a){
			textField.setText(a.toString());
		}
	}
	private class ElementListPanel extends JPanel implements ListSelectionListener{
		private static final long serialVersionUID = 1L;
		private JList list;
		ElementListPanel(){
			super(new BorderLayout());
			super.setBorder(BorderFactory.createTitledBorder("Structure elements"));
			Object[] data = new Object[struc.loopStructures.size()+1];
			int c = 0;
			data[c++] = struc.sheetStructures.get(0);
			for(LoopStructure ls: struc.loopStructures){
				data[c++] = ls;
			}
			list = new JList(data);
			add(new JScrollPane(list));
			list.addListSelectionListener(this);
		}

		@Override
		public void valueChanged(ListSelectionEvent arg0) {
			elementClicked((PartialStructure)list.getSelectedValue());
		}
		
		void elementClicked(PartialStructure ps){
			list.setSelectedValue(ps, true);
			propertiesPanel.elementSelected(ps);
		}

	}
	private class ElementPropertiesPanel extends JPanel{
		private static final long serialVersionUID = 1L;
		JPanel currentPanel;
		ElementPropertiesPanel(){
			super(new BorderLayout());
			super.setBorder(BorderFactory.createTitledBorder("Structure element properties"));
			currentPanel = new JPanel(new GridLayout(3,4));
			add(currentPanel);
		}
		
		void elementSelected(PartialStructure ps){
			remove(currentPanel);
			JPanel newPanel = new JPanel(new GridLayout(3,4));
			if(ps instanceof SurfaceSheetStructure) newPanel = new SheetPanel((SurfaceSheetStructure)ps);
			if(ps instanceof LoopStructure) newPanel = new LoopPanel((LoopStructure)ps);
			
			currentPanel = newPanel;
			add(newPanel);
			frame.validate();
		}
		class LoopPanel extends JPanel{
			private static final long serialVersionUID = 1L;
			LoopStructure ls;
			LoopPanel(LoopStructure loopStruc){
				super(new GridLayout(3,4));
				this.ls = loopStruc;
				JButton b;
				b = new JButton("Enforce closure (CCD)");
				b.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent arg0) {
						ls.enforceClosureCCD();
						ls.updateAtoms(struc);
						scene.repaint();
					}});
				add(b);
				b = new JButton("Rebuild (CCD)");
				b.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent arg0) {
						ls.rebuildCCD();
						ls.updateAtoms(struc);
						scene.repaint();
					}});
				add(b);
			}
		}
		class SheetPanel extends JPanel{
			private static final long serialVersionUID = 1L;
			SurfaceSheetStructure sss;
			SheetPanel(SurfaceSheetStructure sss){
				super(new GridLayout(3,4));
				this.sss = sss;
				
				//slider for c
				ParametricParaboloid pp = (ParametricParaboloid)SheetPanel.this.sss.getSurface();
				int c = (int)(pp.getC()*1000);
				JSlider slider = new JSlider(JSlider.HORIZONTAL, -100,100, c);
				slider.addChangeListener(new ChangeListener(){
					public void stateChanged(ChangeEvent arg0) {
						double newVal = ((JSlider)arg0.getSource()).getValue()/1000.0;
						ParametricParaboloid pp = (ParametricParaboloid)SheetPanel.this.sss.getSurface();
						pp.setC(newVal);
						struc.updateAtoms(struc);
//						for(LoopStructure ls: struc.loopStructures)
//							ls.enforceClosureCCD();
//						struc.updateAtoms(struc);
						scene.repaint();
					}});
				add(slider);
			}
			
		}
	}
	
	public static void main(String[] args) {
		PrimaryStructure ps = new PrimaryStructure("AAEEKTEFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK");//1ctf
		SecondaryStructure ss = new SecondaryStructure(ps, "       EEEEEEE GGGHHHHHHHHHHHH   HHHHHHHHHT SEEEEEEE HHHHHHHHHHHHHHT EEEE ");
		BetaTopology bt = new BetaTopology(ss, new boolean[][]{{false,false,false},{true,false,false},{true,false,false}});
		SheetLoopStructure struc = new SheetLoopStructure(bt);
		new SheetLoopGUI(struc);
	}

}
