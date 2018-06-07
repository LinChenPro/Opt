package my.tests;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;

public class TestSerialisable {
	
	public static void main(String[] args) throws FileNotFoundException, IOException, ClassNotFoundException {
//		String filePath = "C:/.....";
		
//		ObjS os = new ObjS();
//		os.a = 100.123;
//		os.b = 200.456;
//		os.c = new int[]{1,5,9,3,5,7};
//		os.sub = new ObjSub();
//		os.sub.a = "substring....";
//		
//		ObjectOutputStream output = new ObjectOutputStream(new FileOutputStream(filePath));
//		output.writeObject(os);

//		ObjectInputStream input = new ObjectInputStream(new FileInputStream(filePath));
//		ObjS os = (ObjS)input.readObject();
//		
//		System.out.println(os);
		
		String[] s= new String[2000*2000];
		System.out.println(s.length);
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
