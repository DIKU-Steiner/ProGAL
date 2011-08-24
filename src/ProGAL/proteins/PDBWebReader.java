package ProGAL.proteins;

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
public class PDBWebReader {


	public static String readFasta(String pdbId){
		try{
			URL yahoo = new URL("http://www.pdb.org/pdb/explore/sequenceText.do?structureId="+pdbId.toUpperCase()+"&chainId=A");
			BufferedReader in = new BufferedReader(	new InputStreamReader(yahoo.openStream()));

			String inputLine;
			StringBuilder sb = new StringBuilder();

			while ((inputLine = in.readLine()) != null){
				sb.append(inputLine);
				sb.append("\n");
			}

			in.close();

			String ret1 = "";
			String ret2 = "";
			String[] lines = sb.toString().split("\n");
			for(int l=0;l<lines.length;l++){
				if(lines[l].contains("width=\"85%\"")) {
					if(lines[l].contains("entity_poly"))
						ret1+=lines[l];
					else
						ret2+=lines[l];
					//ret+="\n";
				}
			}
			ret1 = ret1.replaceAll("\\<.*?>","");
			ret2 = ret2.replaceAll("\\<.*?>","");
			ret1 = ret1.replaceAll("\\&nbsp;", " ");
			ret2 = ret2.replaceAll("\\&nbsp;", " ");
			StringBuilder sb1 = new StringBuilder(ret1);
			StringBuilder sb2 = new StringBuilder(ret2);
			int i=10;
			while(i<sb1.length()) { 
				sb1.replace(i, i+1, ""); 
				sb2.replace(i, i+1, ""); 
				i+=10; 
			}
			return sb1.toString()+"\n"+sb2.toString();
		}catch(MalformedURLException e){} catch (IOException e) {}

		return null;
	}

	public static String downloadPDBFile(String pdbId){
		return downloadPDBFile(pdbId, true);
	}
	
	/**
	 * Downloads a PDB-file to a temporary file and returns the path
	 */
	public static String downloadPDBFile(String pdbId, boolean displayStatus){
		JFrame dialog = null;
		if(displayStatus){
			dialog = new JFrame("Downloading .. ");
			dialog.setSize(500,100);
			dialog.setLocation(200, 200);
			dialog.getContentPane().add(new JLabel("    Downloading http://www.pdb.org/pdb/files/"+pdbId.toUpperCase()+".pdb"));
			dialog.setVisible(true);
			dialog.setResizable(false);
		}
		try{

			URL url = new URL("http://www.pdb.org/pdb/files/"+pdbId.toUpperCase()+".pdb");
			BufferedReader in = new BufferedReader(	new InputStreamReader(url.openStream()));

			String inputLine;
			StringBuilder sb = new StringBuilder();

			while ((inputLine = in.readLine()) != null){
				sb.append(inputLine);
				sb.append("\n");
			}

			in.close();
			File f = File.createTempFile(pdbId.toUpperCase(), ".pdb");
			f.deleteOnExit();
			IOToolbox.writeToFile(sb.toString(), f.getAbsolutePath(), false);
			if(displayStatus)
				dialog.setVisible(false);
			return f.getAbsolutePath();
		}catch(MalformedURLException e){} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException("Error retrieving "+pdbId);
		}
		dialog.setVisible(false);
		return null;
	}

	public static void main(String[] args){
		String fileName = downloadPDBFile("3DEX");
		System.out.println(fileName);
		System.out.println(new PDBFile(fileName).getSequence());
	}
}
