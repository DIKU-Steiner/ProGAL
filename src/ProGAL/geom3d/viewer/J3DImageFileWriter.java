package ProGAL.geom3d.viewer;

import java.awt.image.BufferedImage;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.media.j3d.Canvas3D;
import javax.media.j3d.GraphicsContext3D;
import javax.media.j3d.ImageComponent;
import javax.media.j3d.ImageComponent2D;
import javax.media.j3d.Raster;
import javax.vecmath.Point3f;

import org.sourceforge.jlibeps.epsgraphics.EpsGraphics2D;

import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGEncodeParam;
import com.sun.image.codec.jpeg.JPEGImageEncoder;

/** 
 * A wrapper class for exporting image files from a J3DScene. A scene is written by 
 * passing the Canvas3D object to one of the methods in this class. The Canvas3D object 
 * is accessed using the getCanvas() method in J3DScene. 
 * @author R. Fonseca
 */
public class J3DImageFileWriter {


	/** Writes the current view in a <code>Canvas3D</code> object to an JPG file */
	public static void writeJPEGFile(String fName, Canvas3D canvas){
		GraphicsContext3D  ctx = canvas.getGraphicsContext3D();
		// The raster components need all be set!
		Raster ras = new Raster(
				new Point3f(-1.0f,-1.0f,-1.0f),
				Raster.RASTER_COLOR,
				0,0,
				canvas.getWidth(),canvas.getHeight(),
				new ImageComponent2D( ImageComponent.FORMAT_RGB, new BufferedImage(canvas.getWidth(), canvas.getHeight(), BufferedImage.TYPE_INT_RGB)),
				null);

		ctx.readRaster(ras);

		// Now strip out the image info
		BufferedImage img = ras.getImage().getImage();

		// write that to disk....
		try {
			FileOutputStream out = new FileOutputStream(fName);
			JPEGImageEncoder encoder = JPEGCodec.createJPEGEncoder(out);
			JPEGEncodeParam param = encoder.getDefaultJPEGEncodeParam(img);
			param.setQuality(0.95f,false); // 75% quality for the JPEG
			encoder.setJPEGEncodeParam(param);
			encoder.encode(img);
			out.close();
		} catch ( IOException e ) {
			e.printStackTrace();
		}
	}
	
	/** Writes the current view in a <code>Canvas3D</code> object to an EPS file */
	public static void writeEPSFile(String fName, Canvas3D canvas){
		GraphicsContext3D  ctx = canvas.getGraphicsContext3D();
		// The raster components need all be set!
		Raster ras = new Raster(
				new Point3f(-1.0f,-1.0f,-1.0f),
				Raster.RASTER_COLOR,
				0,0,
				canvas.getWidth(),canvas.getHeight(),
				new ImageComponent2D( ImageComponent.FORMAT_RGB, new BufferedImage(canvas.getWidth(), canvas.getHeight(), BufferedImage.TYPE_INT_RGB)),
				null);

		ctx.readRaster(ras);

		// Now strip out the image info
		BufferedImage img = ras.getImage().getImage();

		// write that to disk....
		try {
			FileOutputStream finalImage = new FileOutputStream(fName);
			EpsGraphics2D g = new EpsGraphics2D("Title", finalImage, 0, 0, canvas.getWidth(), canvas.getHeight());
			g.drawImage(img, 0, 0, null);
			g.flush();
			g.close();
			finalImage.close();
		} catch ( IOException e ) {
			e.printStackTrace();
		}
	}
}
