package holefilling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

class Model
{
    ArrayList<Vertex> vertices;
    ArrayList<int[]> faces;
    private String header;
    private int oldCountV, oldCountF;
    
    Model()
    {
        vertices = new ArrayList<>();
        faces = new ArrayList<>();
        header = "";
        oldCountV = 0; oldCountF = 0;
    }
    
    boolean ReadFromPLY(String FileName) throws FileNotFoundException, IOException
    {
        File file = new File( "C://Users/Fedorov/Desktop/" + FileName + ".ply");
        if(!file.exists())
            return false;
        
        BufferedReader in = new BufferedReader(new FileReader(file.getAbsoluteFile()));
        String s;
        String[] words;
        while((s = in.readLine()) != null)
        {
            header += s + "\n";
            words = s.split(" ");
            if(words.length > 2 && words[1].equals("vertex"))
                oldCountV = Integer.parseInt(words[2]);
            if(words.length > 2 && words[1].equals("face"))
                oldCountF = Integer.parseInt(words[2]);
            if(words[0].equals("end_header"))
            {
                for(int i = 0; i < oldCountV; i++)
                {
                    words = in.readLine().split(" ");                
                    vertices.add(new Vertex(Double.parseDouble(words[0]),
                            Double.parseDouble(words[1]),
                            Double.parseDouble(words[2])));
                }
               
                for(int i = 0; i < oldCountF; i++)
                {
                       words = in.readLine().split(" ");
                       int[] tmp = {Integer.parseInt(words[1]),
                            Integer.parseInt(words[2]),
                            Integer.parseInt(words[3])};
                    faces.add(tmp);
                }
            }
        }
        in.close();
        return true;
    }
    
    void WriteFromPLY(String FileName) throws FileNotFoundException, IOException
    {
        File file = new File( "C://Users/Fedorov/Desktop/" + FileName + ".ply");
        
        if(!file.exists()) 
            file.createNewFile();
        
        PrintWriter out = new PrintWriter(file.getAbsoluteFile());

        String[] str = header.split("\n");
        for(String s : str)
        {
            if(s.split(" ").length > 2 &&s.split(" ")[1].equals("vertex"))
                s = s.replace(String.valueOf(oldCountV), String.valueOf(vertices.size()));
            if(s.split(" ").length > 2 &&s.split(" ")[1].equals("face"))
                s = s.replace(String.valueOf(oldCountF), String.valueOf(faces.size()));
            out.println(s);
        }
        
        for(Vertex v : vertices)
        {
            out.println(v.x + " " + v.y + " " + v.z);
        }
        
        for(int[] f : faces)
        {
            out.println("3 " + f[0] + " " + f[1] + " " + f[2]);
        }
        
        out.close();
    }
}
