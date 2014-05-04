package ProGAL.proteins.viewer;

import java.awt.Color;
import java.util.List;

import javax.swing.JFrame;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.proteins.PDBWebReader;
import ProGAL.proteins.PDBFile.PDBChain;
import ProGAL.proteins.PDBFile.PDBModel;
import ProGAL.proteins.belta.PDBFile;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;

public class PDBFileViewer{
	public static final int STICKS = 0x1;
	public static final int CARTOON = 0x2;
	public static final int SPHERES = 0x4;

	public final PDBFile pdbFile;
	public final J3DScene scene;
	private final JFrame frame;
	private int displayMode = CARTOON;


	public PDBFileViewer(PDBFile f){
		this.pdbFile = f;

		scene = J3DScene.createJ3DSceneInFrame();
		this.frame = scene.frame;
		frame.setTitle("ProGAL viewer: "+f.name);

		displayStructure(f);
	}

	public void displayStructure(PDBFile f){
		scene.removeAllShapes();

		if( (displayMode&CARTOON)!=0 ){
//			for(int m = 0;m<f.getModels().size();m++){
//				for(int c=0;c<f.getModels().get(m).getChains().size();c++){
			int m = 0, c = 0;
			{ 
				{
					List<Point> cAlphas = f.getCACoords(m, c);
					int end = cAlphas.size();
//					if(end<3) continue;

					SecondaryStructure ss = f.getSecondaryStructure();
					for(SSSegment seg: ss.segments){
						switch(seg.type){
						case HELIX: 
							HelixSurface hSurf = new HelixSurface(cAlphas.subList(seg.start, Math.min(seg.end+1,end)));
							scene.addSurface(hSurf, new Color(200,10,10), 0, seg.length-0.001, seg.length*10, 0,1,7);
							break;
						case STRAND:
							StrandSurface sSurf = new StrandSurface(cAlphas.subList(seg.start, Math.min(seg.end+1,end)));
							scene.addSurface(sSurf, new Color(200,200,10), 0,seg.length-1, seg.length*10, 0,1,7);
							break;
						case COIL:
							CoilSurface cSurf = new CoilSurface(cAlphas.subList(seg.start, Math.min(seg.end+1,end)));
							scene.addSurface(cSurf, new Color(200,200,200), 0,seg.length, seg.length*10, 0,1,4);
							break;
						}
					}
				}
			}
		}
		
		scene.centerCamera();
		scene.autoZoom();
	}
	public static void main(String[] args){
		PDBFile f = new PDBFile(PDBWebReader.downloadPDBFile("2CRO"));
		new PDBFileViewer(f);
	}
}
