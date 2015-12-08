package ProGAL.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * A wrapper class for static file input/output methods. All of the methods have a policy of not
 * throwing IOExceptions (and similar) if an error occurs. Instead they throw RuntimeExceptions 
 * such that try-catch-statements can be excluded if desired. 
 * 
 * All methods (TODO: Only implemented for writeToFile and readFromFile) supports reading (writing)
 * from (to) gzipped files. This feature is implicitly activated when the specified file name ends 
 * with '.gz'.
 * 
 * @author R.Fonseca
 */
public abstract class IOToolbox {
	
	public static boolean createDirectory(String path){
		  File theDir = new File(path);
		  if (!theDir.exists()) {
		    return theDir.mkdir();
		  }
		  return true;
	}
	
	/** 
	 * Writes a string to the specified file-name. 
	 * 
	 * If append is true and the file already exists then the contents are appended to 
	 * the end of the file. If anything goes wrong (e.g. insufficient permissions or bad 
	 * file-name) a stack-trace is written and a RuntimeException is thrown (so that no 
	 * try-catch is required).
	 * 
	 * If fileName ends with '.gz' then the file is zipped. 
	 */
	public static void writeToFile(String contents, String fileName, boolean append){
		boolean zip = fileName.endsWith(".gz");
		try{
			BufferedWriter out;
			if(zip){
				String appendString = null;
				if(append && new java.io.File(fileName).exists()) appendString = readFromFile(fileName);

				out = new BufferedWriter(
						new OutputStreamWriter( 
								new GZIPOutputStream(
										new FileOutputStream(fileName, false) 
										) 
								)
						);

				if(appendString!=null) out.write(appendString);
			} else {
				out = new BufferedWriter(new FileWriter(fileName, append));
			}

			out.write(contents);
			out.close();
		}catch (Exception e){//Catch exception if any
			e.printStackTrace();
			throw new RuntimeException("Error writing file "+fileName);
		}
	}

	/** 
	 * Read the contents of the specified file and return it in a string. 
	 * 
	 * If anything goes wrong (e.g. insufficient permissions or bad file-name) 
	 * a stack-trace is written and a RuntimeException is thrown (so that no 
	 * try-catch is required).
	 * 
	 * If fileName ends with '.gz' then the file is assumed to be zipped. 
	 */
	public static String readFromFile(String fName){
		boolean zip = fName.endsWith(".gz");
		try{
			BufferedReader in;
			if(zip){
				in = new BufferedReader(
						new InputStreamReader(	
								new GZIPInputStream(new FileInputStream( fName )	
										) 
								) 
						);
			}else{
				in = new BufferedReader(new FileReader(fName));
			}


			StringBuilder sb = new StringBuilder();
			String line;
			while( (line=in.readLine())!=null ){
				sb.append(line);
				sb.append('\n');
			}

			in.close();

			return sb.toString();
		}catch(Exception exc){
			exc.printStackTrace();
			throw new RuntimeException("Error reading file "+fName);
		}
	}

	/** 
	 * Read the contents of the specified file and return it in a string. This file is 
	 * assumed to be located within the project, e.g. wrapped in a jar-file with the 
	 * main method.
	 * 
	 * If something goes wrong (e.g. insufficient permissions or bad file-name) a 
	 * stack-trace is written and a RuntimeException is thrown (so that no try-catch 
	 * is required).
	 */
	public static String readFromInternalFile(String fName) throws Error{ 

		try{
			StringBuilder sb = new StringBuilder();
			InputStream input = IOToolbox.class.getResourceAsStream("/"+fName);
			if(input==null) throw new Error("No such resource");
			BufferedReader in = new BufferedReader(new InputStreamReader(input));
			String line;
			while( (line=in.readLine())!=null ){
				sb.append(line);
				sb.append('\n');
			}
			in.close();
			return sb.toString();
		}catch(Exception exc){
			exc.printStackTrace();

			throw new RuntimeException("Error reading file "+fName);
		}
	}



	/** 
	 * Writes a serialized object to the specified file-name. If something goes wrong 
	 * (e.g. insufficient permissions or bad file-name) a stack-trace is written and a 
	 * RuntimeException is thrown (so that no try-catch is required).
	 */
	public static void writeSerializedFile(Serializable obj, String fName){
		try {
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(fName));
			oos.writeObject(obj);
			oos.close();
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException("Error writing file "+fName);
		}
	}

	/** 
	 * Reads a serialized object from the specified file-name.  If something goes wrong 
	 * (e.g. insufficient permissions or bad file-name) a stack-trace is written and a 
	 * RuntimeException is thrown (so that no try-catch is required).
	 */
	public static Object readSerializedFile(String fName){
		Object o = null;
		try{
			ObjectInputStream ois = new ObjectInputStream(new FileInputStream(fName));
			o = ois.readObject();
			ois.close();
		}catch(Exception ex){
			ex.printStackTrace();
			throw new RuntimeException("Error reading file "+fName);			
		}
		return o;
	}

}
