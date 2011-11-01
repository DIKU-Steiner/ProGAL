package ProGAL.proteins;

import javax.swing.JPanel;

import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.proteins.PDBFile.PDBModel;

public class PDBFileViewer {
	private final J3DScene scene;
	
	private PDBFileViewer(PDBFile f){
		scene = J3DScene.createJ3DSceneInFrame();

		for(PDBModel m: f.getModels()){
			
		}
			
	}
	
	
	private class ViewPanel extends JPanel{
		
	}
	
	
	
	public static J3DScene showPDBFile(PDBFile f){
		PDBFileViewer fv = new PDBFileViewer(f);
		return fv.scene;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		showPDBFile(PDBFile.downloadPDBFile("2CRO"));
	}

}
