package ProGAL.geom3d.viewer;

import javax.media.j3d.*;

import ProGAL.geom3d.*;
import ProGAL.geom3d.surface.ParametricSurface;
import ProGAL.math.Constants;

/** 
 * @author R.Fonseca
 */
class Surface3D extends Shape3D {
	private final ParametricSurface surface;
	private final double uMin, uMax, vMin, vMax;
	private final int uDiv, vDiv;

	public Surface3D(ParametricSurface surf, double uMin, double uMax, int uDiv, double vMin, double vMax, int vDiv, Appearance app) {
		super();

		this.surface = surf;
		this.uMin = uMin;
		this.uMax = uMax;
		this.vMin = vMin;
		this.vMax = vMax;
		this.uDiv = uDiv;
		this.vDiv = vDiv;
		
		setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
		update();

		if(app==null)
			setAppearance(new Appearance());
		else
			setAppearance(app);
	}

	void update(){
		Point[][] points = new Point[uDiv+1][vDiv+1];
		Vector[][] norms = new Vector[uDiv+1][vDiv+1];

		double du = (uMax-uMin)/uDiv-Constants.EPSILON;
		double dv = (vMax-vMin)/vDiv-Constants.EPSILON;

		double u, v = vMin;
		int uInt = 0, vInt = 0;
		for(u = uMin; u<=uMax-du;u+=du){
			vInt=0;
			for(v = vMin; v<=vMax-dv;v+=dv){
				points[uInt][vInt] = surface.getPoint(u,v);
				norms[uInt][vInt] = surface.getNormal(u,v);

				vInt++;
			}
			v = vMax;
			points[uInt][vInt] = surface.getPoint(u,v);
			norms[uInt][vInt] = surface.getNormal(u,v);

			uInt++;
		}
		u = uMax;
		vInt = 0;
		for(v = vMin; v<=vMax-dv;v+=dv){
			points[uInt][vInt] = surface.getPoint(u,v);
			norms[uInt][vInt] = surface.getNormal(u,v);

			vInt++;
		}
		v = vMax;
		points[uInt][vInt] = surface.getPoint(u,v);
		norms[uInt][vInt] = surface.getNormal(u,v);

		float[] vertArr = new float[(uDiv)*(vDiv)*4*3];
		float[] normArr = new float[(uDiv)*(vDiv)*4*3];
		int c = 0;
		for(int i=0;i<uDiv;i++){
			for(int j=0;j<vDiv;j++){
				vertArr[c++] = (float)points[i][j].x();
				vertArr[c++] = (float)points[i][j].y();
				vertArr[c++] = (float)points[i][j].z();
				vertArr[c++] = (float)points[i][j+1].x();
				vertArr[c++] = (float)points[i][j+1].y();
				vertArr[c++] = (float)points[i][j+1].z();
				vertArr[c++] = (float)points[i+1][j+1].x();
				vertArr[c++] = (float)points[i+1][j+1].y();
				vertArr[c++] = (float)points[i+1][j+1].z();
				vertArr[c++] = (float)points[i+1][j].x();
				vertArr[c++] = (float)points[i+1][j].y();
				vertArr[c++] = (float)points[i+1][j].z();

			}
		}
		c=0;
		for(int i=0;i<uDiv;i++){
			for(int j=0;j<vDiv;j++){
				normArr[c++] = (float)norms[i][j].x();
				normArr[c++] = (float)norms[i][j].y();
				normArr[c++] = (float)norms[i][j].z();
				normArr[c++] = (float)norms[i][j+1].x();
				normArr[c++] = (float)norms[i][j+1].y();
				normArr[c++] = (float)norms[i][j+1].z();
				normArr[c++] = (float)norms[i+1][j+1].x();
				normArr[c++] = (float)norms[i+1][j+1].y();
				normArr[c++] = (float)norms[i+1][j+1].z();
				normArr[c++] = (float)norms[i+1][j].x();
				normArr[c++] = (float)norms[i+1][j].y();
				normArr[c++] = (float)norms[i+1][j].z();

			}
		}

		QuadArray quads = new QuadArray(vertArr.length/3, QuadArray.COORDINATES | QuadArray.NORMALS);
		quads.setCoordinates(0,vertArr);
		quads.setNormals(0,normArr);
		quads.setCapability(Geometry.ALLOW_INTERSECT);
		setGeometry(quads);

	}

}