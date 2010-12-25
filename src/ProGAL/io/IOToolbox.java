package ProGAL.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

/**
 * A wrapper class for static file input/output methods. All of the methods have a policy of not
 * throwing IOExceptions (and similar) if an error occurs. Instead they throw RuntimeExceptions 
 * such that try-catch-statements can be excluded if desired. 
 * @author R.Fonseca
 */
public abstract class IOToolbox {
	
	/** 
	 * Writes a string to the specified file-name. If append is true and the file already 
	 * exists then the contents are appended to the end of the file. If something goes wrong 
	 * (e.g. insufficient permissions or bad file-name) a stack-trace is written and a 
	 * RuntimeException is thrown (so that no try-catch is required).
	 */
	public static void writeToFile(String contents, String fName, boolean append){
		try{
			BufferedWriter out = new BufferedWriter(new FileWriter(fName, append));
			out.write(contents);
			out.close();
		}catch (Exception e){//Catch exception if any
			e.printStackTrace();
			throw new RuntimeException("Error writing file "+fName);
		}

	}

	/** 
	 * Read the contents of the specified file and return it in a string. If something goes wrong 
	 * (e.g. insufficient permissions or bad file-name) a stack-trace is written and a 
	 * RuntimeException is thrown (so that no try-catch is required).
	 */
	public static String readFromFile(String fName){
		try{
			StringBuilder sb = new StringBuilder();

			BufferedReader in = new BufferedReader(new FileReader(fName));
			String line;
			while( (line=in.readLine())!=null ){
				sb.append(line);
				sb.append('\n');
			}

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
			throw new RuntimeException("Error writing file "+fName);			
		}
		return o;
	}
}
