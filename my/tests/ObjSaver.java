package my.tests;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

public class ObjSaver<T>{
	public static String folderPath;
	
	public static void setFolderPath(String folderPath){
		ObjSaver.folderPath = folderPath;
	}
	
	public static String getClassNameForFile(Class<?> clazz){
		return clazz.getSimpleName().replaceAll("(\\[\\])+", "_arrays");
	}
	
	private static String getFileName(Object obj, String name){
		return getFileName(getClassNameForFile(obj.getClass()), name);
	}
	
	private static String getFileName(String className, String name){
		return className+"_"+name.replaceAll("\\s+", "_").toLowerCase()+".obj";
	}
	
	public static <E> E readObj(E obj, String name){
		try{
			String filePath = folderPath+"/"+getFileName(obj, name);
			ObjectInputStream input = new ObjectInputStream(new FileInputStream(filePath));
			ObjS os = (ObjS)input.readObject();
			input.close();
			return (E)os;
		}catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	public static Object readObjWhenNeed(Class clazz, Object tagObj, String name){
		if(tagObj != null){
			return tagObj;
		}
		
		return readObj(clazz, name);
	}

	public static Object readObj(Class<?> clazz, String name){
		try{
			String filePath = folderPath+"/"+getFileName(getClassNameForFile(clazz), name);
			ObjectInputStream input = new ObjectInputStream(new FileInputStream(filePath));
			Object obj = input.readObject();
			input.close();
			return obj;
		}catch (Exception e) {
			System.out.println("stored " + clazz.getName() + " obj " + name + " not found");
			return null;
		}
	}

	public T readObjT(Class clazz, String name){
		try{
			String filePath = folderPath+"/"+getFileName(getClassNameForFile(clazz), name);
			ObjectInputStream input = new ObjectInputStream(new FileInputStream(filePath));
			ObjS os = (ObjS)input.readObject();
			input.close();
			return (T)os;
		}catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	public static boolean saveObj(Object obj, String name){
		if(obj==null){
			return false;
		}
		
		String filePath = folderPath+"/"+getFileName(obj, name);
		ObjectOutputStream output;
		try {
			output = new ObjectOutputStream(new FileOutputStream(filePath));
			output.writeObject(obj);
			output.close();
			return true;
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException, ClassNotFoundException {
		String folderPath = "/home/chen/searchprojets/resource";
		ObjSaver.setFolderPath(folderPath);
		
//		ObjS os = new ObjS();
//		os.a = 100.123;
//		os.b = 200.456;
//		os.c = new int[]{1,5,9,3,5,7};
//		os.sub = new ObjSub();
//		os.sub.a = "substring....";
//
//		saver.saveObj(os, "os a");
		
		ObjS os = (ObjS)ObjSaver.readObj(ObjS.class, "os a");

		System.out.println(os);
		
	}

	public static class ObjSub implements Serializable{
		String a;
	}
	
	
	public static class ObjS implements Serializable{
		private double a;
		double b;
		int[] c;
		ObjSub sub;
		
		private static final long serialVersionUID = 1L;
		
//		public double getA() {
//			return a;
//		}
//		public void setA(double a) {
//			this.a = a;
//		}
//		public void setB(double b) {
//			this.b = b;
//		}
//		public int[] getC() {
//			return c;
//		}
//		public void setC(int ... c) {
//			this.c = c;
//		}
		
		public String toString(){
			String sc = "null";
			if(c!=null){
				sc="[";
				String sep = "";
				for(int ci : c){
					sc += sep+ci;
					sep = ",";
				}
				sc+="]";
			}
			
			
			return "a="+a+"\nb="+b+"\nc="+sc+"\nsub="+sub.a;
		}
	}
}
