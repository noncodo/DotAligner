import java.io.*; 
import java.lang.Math ; 
public class ps2pp {
    public static void main(String[] Args) throws Exception  {
        File Folder = new File( Args[ 0 ] );
        String InputDIR = Folder.getAbsolutePath() ; 
        System.err.println( "[[ ps2pp ]] reading .ps files in dir: "+ InputDIR ); 
        File[] ListOfFiles = Folder.listFiles();
        BufferedReader ReadFile ; 
        BufferedWriter WriteFile ; 
        double [][] PP ; 

        for (File File : ListOfFiles) {
            String FName = File.getName() ; 
            if ( File.isFile() 
                 && FName.substring( FName.length() -5 ).equals("dp.ps" ) ) {
    		
                ReadFile = new BufferedReader(new FileReader( File ));
                WriteFile = new BufferedWriter( new FileWriter(InputDIR + "/" + FName.substring(0, FName.lastIndexOf(".") ) +".pp")) ;

        		String Line = "" ; 
        		String Sequence = "" ; 
        		PP = new double [1][1]; 
        		read : while ( (Line = ReadFile.readLine()) != null )  {
        			
        			if (Line.length() > 8 && Line.substring(0,9).equals("/sequence")) {
        				Sequence = ReadFile.readLine() ;
        				Sequence = Sequence.substring(0, Sequence.length()-1 ); 
        				while ( (Line = ReadFile.readLine()).charAt(0) != ')' ) {
        					Sequence = Sequence + Line ;
        					Sequence = Sequence.substring(0, Sequence.length()-1 ) ; 
        				}
        				WriteFile.write( Sequence.length() + "\n" ); 
                        WriteFile.write(  Sequence.toUpperCase() + "\n"); 
        				PP = new double[ Sequence.length() ][ Sequence.length() ] ; 
        			}
        			else if ( Line.length() != 0 && Line.matches(".*ubox") && Line.charAt(0) != '%' && Line.charAt(0) != '/'){
        				String [] Fields = Line.split(" "); 
        				//System.out.println( Line + Sequence.length() ) ; 
        				int x = Integer.parseInt( Fields[ 0 ]); 
        				int y = Integer.parseInt( Fields[ 1 ]); 
        				PP [ x-1 ][ y-1 ] = Math.pow( Double.parseDouble( Fields[ 2 ] ), 2) ;
        				PP [ y-1 ][ x-1 ] = Math.pow( Double.parseDouble( Fields[ 2 ] ), 2) ; 
        			}
        		}
        		ReadFile.close(); 
                // write it out
                // N.B. this creates pretty large files. Encode differently? 
        		double [] UP = new double [ PP.length ]; 
        		for ( int x = 0 ; x != PP.length ; x++ ) {
        			for ( int y = 0 ; y != PP.length ; y++ ) {
                        WriteFile.write( PP[ x ][ y ] +" ");
        				UP[x] = UP[x]+PP[ x ][ y ] ; 
        			}
                    WriteFile.write( "\n" );
    			}
                WriteFile.write( "\n" );

        		for ( int x = 0 ; x != UP.length ; x++ ) 
        			//System.out.print( 1-UP[ x ] +" " );
                    WriteFile.write( 1-UP[ x ] +" " ); 
                WriteFile.write( "\n" );
                WriteFile.close(); 
            }
    	}
	}
}

