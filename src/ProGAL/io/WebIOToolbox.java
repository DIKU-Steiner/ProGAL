package ProGAL.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;

import javax.swing.JFrame;
import javax.swing.JLabel;

import ProGAL.io.IOToolbox;

/** 
 * A wrapper class for static web input/output methods.
 * @author R. Fonseca
 */
public class WebIOToolbox {
	/** 
	 * Downloads an online document to a temporary file and returns its path. While the file 
	 * is being downloaded a dialog box is shown indicating the URL. When the thread exits 
	 * the temporary file is removed.
	 * 
	 * If anything goes wrong a stack-trace is printed and a RuntimeException is thrown.
	 */
	public static String downloadFile(String urlString){
		JFrame dialog = new JFrame("Downloading .. ");
		dialog.setSize(500,100);
		dialog.setLocation(200, 200);
		dialog.getContentPane().add(new JLabel("    Downloading "+urlString));
		dialog.setVisible(true);
		dialog.setResizable(false);
		try{

			URL url = new URL(urlString);
			BufferedReader in = new BufferedReader(	new InputStreamReader(url.openStream()));

			String inputLine;
			StringBuilder sb = new StringBuilder();

			while ((inputLine = in.readLine()) != null){
				sb.append(inputLine);
				sb.append("\n");
			}

			in.close();
			File f = File.createTempFile("ProGAL_download_"+System.currentTimeMillis(), ".tmp");
			f.deleteOnExit();
			IOToolbox.writeToFile(sb.toString(), f.getAbsolutePath(), false);
			dialog.setVisible(false);
			dialog.dispose();
			return f.getAbsolutePath();
		}catch(MalformedURLException e){} catch (IOException e) {
			e.printStackTrace();
			dialog.setVisible(false);
			dialog.dispose();
			throw new RuntimeException("Error retrieving "+urlString);
		}
		dialog.setVisible(false);
		dialog.dispose();
		return null;
	}
}
