package holefilling;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import com.vividsolutions.jts.triangulate.DelaunayTriangulationBuilder;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

class Hole
{
    Vertex[] boundary;
    //ArrayList<ArrayList<Vertex>> close;
    Vertex[] close;
    ArrayList<Vertex> repair;
    int[] num;
}

public class HoleFilling {
    
    public static void main(String[] args) throws IOException
    {
        Model model = new Model();
        model.ReadFromPLY("ModelWithHole");
        //model.WriteFromPLY("myPhone");
        //MethodDetectingHoles(model);
        MethodFillingHoles(model);
    }
    
    static void MethodDetectingHoles(Model model)
    {
        Neighborhood(model);
    }
    
    static void Neighborhood(Model model)
    {
        boolean flag = true;
        ArrayList<Vertex>[] q = new ArrayList[model.vertices.size()];
        for(int i = 0; i < q.length; i++)
        {
            q[i] = new ArrayList();
            for(int[] face : model.faces)
            for(int j = 0; j < face.length; j++)
            {
                if(face[j] == i)
                {
                    for(int jj = 0; jj < face.length; jj++)
                        if(face[jj] != i)
                        {
                            for(int k = 0; k < q[i].size(); k++)
                            {
                                if(CompareP(q[i].get(k), model.vertices.get(face[jj])))
                                {
                                    flag = false;
                                    break;
                                }
                            }
                            if(flag)
                                q[i].add(model.vertices.get(face[jj]));
                            flag = true;
                        }
                }
            }
        }
       
        double[] R = new double[model.vertices.size()];
        double l = 0;
        
       // for(int i = 0; i < model.vertices.size(); i++)
       // {
         //   for(Vertex v : q[i])
           //     l += LengthEdge(model.vertices.get(i), v);
        //    R[i] = l / 
       // }
        
    }
      
    static void MethodFillingHoles(Model model) throws IOException
    {
        //Hole hole = ReadHole("hole");
        
        Model h = new Model();
        h.ReadFromPLY("Hole");
        Model border = new Model();
        border.ReadFromPLY("HoleBorder");
        
        //border.vertices = sortVertex(border, h);
        //border.WriteFromPLY("newBorder");
        
        //Hole hole = ReadHole(model, h, border);
        Hole hole = ReadHole("hole", model);
        //vibor(hole);
        
       // Model hh = new Model();
        
        //for(Vertex v : hole.boundary)
           // hh.vertices.add(v);
        
        //hh.WriteFromPLY("My");
        
        Vertex[] fp = new Vertex[4];
        Vertex tmpV;
        double st;
        int t = 0;
        hole.repair = new ArrayList();
        //for(int i = 0; i < hole.boundary.length/2; i++)
        for(int i = 0; i < hole.boundary.length / 2; i = i + 1)
            //int  i = 0;
        //for(int j = hole.boundary.length/2; j < hole.boundary.length; j++)
        {
           // if(i != j && i + 1 != j && i - 1 != j)
            //{
            //for(int i1 = 0; i1 < hole.close.get(i).size(); i1++)
            //for(int i2 = 0; i2 < hole.close.get
                //(i + hole.boundary.length / 2).size(); i2++)
            //for(int i2 = 0; i2 < hole.close.get(j).size(); i2++)
            //int i1, i2;
            //if(hole.close.get(i).size() > 0 &&
                  //  hole.close.get(i + hole.boundary.length / 2).size() > 0)
            {
                //i1 = 0; i2 = 0;
            
                //fp[0] = hole.close.get(i).get((int)((hole.close.get(i).size() - 1)/1.0));
                fp[0] = hole.close[i];
                fp[1] = hole.boundary[i];
                //fp[2] = hole.boundary[j];
                //fp[3] = hole.close.get(j).get(i2);
                fp[2] = hole.boundary[i + hole.boundary.length / 2];
                fp[3] = hole.close[i + hole.boundary.length / 2];
                //fp[3] = hole.close.get(i + hole.boundary.length / 2).
                        //get((int)((hole.close.get(i + hole.boundary.length / 2).size() - 1)/1.0));
                
                //double ll = LengthEdge(fp[1], fp[2]);
                
               // fp[0].x *= ll; fp[0].y *= ll;
                //fp[0].z *= ll;
               // fp[3].x *= ll; fp[3].y *= ll;
                //fp[3].z *= ll;
                t = CountT(hole.boundary, fp) + 1;
                double step = LengthEdge(fp[1], fp[2])/t;                
                st = 0;
                for(int k = 0; k < t; k++)
                {
                    //st += 1.0/(t + 1);
                    st += step;
                    tmpV = CatmullRom2(st, fp);
                   // tmpV.x /= ll; tmpV.y /= ll;
                    //tmpV.z /= ll;
                    hole.repair.add(tmpV);
                    model.vertices.add(tmpV);
                    //hh.vertices.add(tmpV);
                }
            }
            //}
        }
        
        Vertex[] tp = new Vertex[hole.boundary.length + hole.repair.size()];        
        
        for(int i = 0; i < hole.boundary.length; i++)
            tp[i] = hole.boundary[i];
        for(int i = 0; i < hole.repair.size(); i++)
            tp[i + hole.boundary.length] = hole.repair.get(i);
                
        int[][] outF = BestTriangular(tp);
        int[] tmpFace;
        
        boolean norm = true;
        int tmp;
        
        for(int[] face : outF)
        {
            tmpFace = new int[3];
            for(int i = 0; i < face.length; i++)
            {
                if(face[i] < hole.num.length)
                    tmpFace[i] = hole.num[face[i]]; 
                
                if(face[i] >= hole.num.length)
                    tmpFace[i] = (face[i] - hole.boundary.length)
                            + (model.vertices.size() - hole.repair.size());
            }
            if(norm)
            {
                tmp = tmpFace[0];
                tmpFace[0] = tmpFace[2];
                tmpFace[2] = tmp;   
            }
            model.faces.add(tmpFace);
        }
        
        model.WriteFromPLY("NewModel");
       // hh.WriteFromPLY("My");
    }
    
   /* static void vibor(Hole hole)
    {
        double min;
        Vertex n;
        for(int i = 0; i < hole.boundary.length; i++)
        {
            min = LengthEdge(hole.boundary[i], hole.close.get(i).get(0));
            n = hole.close.get(i).get(0);
            for(int j = 0; j < hole.close.get(i).size(); j++)
            {
                if(hole.close.get(i).size() < 2)
                    break;
                else
                {
                    if(LengthEdge(hole.boundary[i], hole.close.get(i).get(j)) < min)
                    {
                        min = LengthEdge(hole.boundary[i], hole.close.get(i).get(j));
                        n = hole.close.get(i).get(j);
                    }
                }
            }
            hole.close.get(i).clear();
            hole.close.get(i).add(n);
        }
    }*/
    
    static boolean CompareP(Vertex a, Vertex b)
    {
        return (a.x == b.x && a.y == b.y && a.z == b.z);
    }
    
    static int[][] BestTriangular(Vertex[] tp)
    {
        double[][] v = new double[2][tp.length];
        Vertex[] tmp_tp = new Vertex[tp.length];
        
        for(int i = 0; i < tp.length; i++)
            tmp_tp[i] = Projection(0, 0, 1, 0, tp[i]);
        
        for(int i = 0; i < tmp_tp.length; i++)
        {
            v[0][i] = tmp_tp[i].x;
            v[1][i] = tmp_tp[i].y;
        }
        
        int[][] outF1 = Triangular(v);
        
        for(int i = 0; i < tp.length; i++)
            tmp_tp[i] = Projection(0, 1, 0, 0, tp[i]);
        
        for(int i = 0; i < tmp_tp.length; i++)
        {
            v[0][i] = tmp_tp[i].x;
            v[1][i] = tmp_tp[i].z;
        }
        
        int[][] outF2 = Triangular(v);
        
        for(int i = 0; i < tp.length; i++)
            tmp_tp[i] = Projection(1, 0, 0, 0, tp[i]);
        
        for(int i = 0; i < tmp_tp.length; i++)
        {
            v[0][i] = tmp_tp[i].y;
            v[1][i] = tmp_tp[i].z;
        }
        
        int[][] outF3 = Triangular(v);
        
        if(outF1.length > outF2.length)
        {
            if(outF1.length > outF3.length)
                return outF1;
            else return outF3;
        }
        else
        {
            if(outF2.length > outF3.length)
                return outF2;
            else return outF3;
        }
    }
    
    static int[][] Triangular(double[][] v)
    {
        ArrayList<Coordinate> coords = new ArrayList();
        for(int i = 0; i < v[0].length; i++)
            coords.add(new Coordinate(v[0][i], v[1][i]));
       
        DelaunayTriangulationBuilder builder = new DelaunayTriangulationBuilder();
        builder.setSites(coords);
        Geometry geometry = builder.getTriangles(new GeometryFactory());
        
        int[][] faces = new int[geometry.getNumGeometries()][3];
        double[] tmp = new double[2];
        for(int i = 0; i < faces.length; i++)
        {
            for(int j = 0; j < faces[i].length; j++)
            {
                //tmp = new Vertex(geometry.getGeometryN(i).getCoordinates()[j].x,
                   // geometry.getGeometryN(i).getCoordinates()[j].y, 0);
                tmp[0] = geometry.getGeometryN(i).getCoordinates()[j].x;
                tmp[1] = geometry.getGeometryN(i).getCoordinates()[j].y;
                for(int k = 0; k < v[0].length; k++)
                {
                    if(tmp[0] == v[0][k] && tmp[1] == v[1][k])
                        faces[i][j] = k;
                }
            }
        }
        return faces;
    }
    
    static double dot(Vertex a, Vertex b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }
    
    static Vertex Projection(double a, double b, double c, double d, Vertex v)
    {
        Vertex N = new Vertex(a, b, c);
        double k = -(dot(N, v) + d)/dot(N, N);
        return new Vertex(v.x + k*N.x, v.y + k*N.y, v.z + k*N.z);
    }
    
    static int CountT(Vertex[] bp, Vertex[] fp)
    {
        double average = LengthBoundary(bp) / bp.length;
        return (int)Math.round(LengthEdge(fp[1], fp[2]) / average);
    }
    
    static Vertex CatmullRom(double t, Vertex[] fp)
    {
        double t1 = 0;
        
        Vertex v;
        double s = 1 - t;
        double t2 = t*t;
        double t3 = t2*t;
        v = new Vertex(0.5*((-t)*s*s*fp[0].x + (2 - 5*t2 + 3*t3)*fp[1].x + t*(1 + 4*t - 3*t2)*fp[2].x - t2*s*fp[3].x),
                0.5*((-t)*s*s*fp[0].y + (2 - 5*t2 + 3*t3)*fp[1].y + t*(1 + 4*t - 3*t2)*fp[2].y - t2*s*fp[3].y),
                0.5*((-t)*s*s*fp[0].z + (2 - 5*t2 + 3*t3)*fp[1].z + t*(1 + 4*t - 3*t2)*fp[2].z - t2*s*fp[3].z));
        /*v = new Vertex(0.5*((-t)*Math.pow(1-t, 2)*fp[0].x+(2 - 5*t*t + 3*t*t*t)*fp[1].x + 
                t*(1 + 4*t - 3*t*t)*fp[2].x - t*t*(1-t)*fp[3].x),
                0.5*((-t)*Math.pow(1-t, 2)*fp[0].y+(2 - 5*t*t + 3*t*t*t)*fp[1].y + 
               t*(1 + 4*t - 3*t*t)*fp[2].y - t*t*(1-t)*fp[3].y),
                0.5*((-t)*Math.pow(1-t, 2)*fp[0].z+(2 - 5*t*t + 3*t*t*t)*fp[1].z + 
                t*(1 + 4*t - 3*t*t)*fp[2].z - t*t*(1-t)*fp[3].z));*/
        return v;
    }
    
    
    static Vertex CatmullRom2(double t, Vertex[] fp)
    {
        double[] ti = new double[4];
        ti[0] = 0;
        for(int i = 0; i < 3; i++)
            ti[i+1] = ti[i] + LengthEdge(fp[i + 1], fp[i]);
        
        Vertex A1 = new Vertex((ti[1]-t)/(ti[1]-ti[0])*fp[0].x + (t-ti[0])/(ti[1]-ti[0])*fp[1].x,
                               (ti[1]-t)/(ti[1]-ti[0])*fp[0].y + (t-ti[0])/(ti[1]-ti[0])*fp[1].y,
                               (ti[1]-t)/(ti[1]-ti[0])*fp[0].z + (t-ti[0])/(ti[1]-ti[0])*fp[1].z);
        
        Vertex A2 = new Vertex((ti[2]-t)/(ti[2]-ti[1])*fp[1].x + (t-ti[1])/(ti[2]-ti[1])*fp[2].x,
                               (ti[2]-t)/(ti[2]-ti[1])*fp[1].y + (t-ti[1])/(ti[2]-ti[1])*fp[2].y,
                               (ti[2]-t)/(ti[2]-ti[1])*fp[1].z + (t-ti[1])/(ti[2]-ti[1])*fp[2].z);
        
        Vertex A3 = new Vertex((ti[3]-t)/(ti[3]-ti[2])*fp[2].x + (t-ti[2])/(ti[3]-ti[2])*fp[3].x,
                               (ti[3]-t)/(ti[3]-ti[2])*fp[2].y + (t-ti[2])/(ti[3]-ti[2])*fp[3].y,
                               (ti[3]-t)/(ti[3]-ti[2])*fp[2].z + (t-ti[2])/(ti[3]-ti[2])*fp[3].z);

        Vertex B1 = new Vertex((ti[2]-t)/(ti[2]-ti[0])*A1.x + (t-ti[0])/(ti[2]-ti[0])*A2.x,
                               (ti[2]-t)/(ti[2]-ti[0])*A1.y + (t-ti[0])/(ti[2]-ti[0])*A2.y,
                               (ti[2]-t)/(ti[2]-ti[0])*A1.z + (t-ti[0])/(ti[2]-ti[0])*A2.z);
        
        Vertex B2 = new Vertex((ti[3]-t)/(ti[3]-ti[1])*A2.x + (t-ti[1])/(ti[3]-ti[1])*A3.x,
                               (ti[3]-t)/(ti[3]-ti[1])*A2.y + (t-ti[1])/(ti[3]-ti[1])*A3.y,
                               (ti[3]-t)/(ti[3]-ti[1])*A2.z + (t-ti[1])/(ti[3]-ti[1])*A3.z);

        Vertex C  = new Vertex((ti[2]-t)/(ti[2]-ti[1])*B1.x + (t-ti[1])/(ti[2]-ti[1])*B2.x,
                               (ti[2]-t)/(ti[2]-ti[1])*B1.y + (t-ti[1])/(ti[2]-ti[1])*B2.y,
                               (ti[2]-t)/(ti[2]-ti[1])*B1.z + (t-ti[1])/(ti[2]-ti[1])*B2.z);
        
        return C;
    }
    
    static double LengthEdge(Vertex a, Vertex b)
    {
        return Math.sqrt(Math.pow((a.x - b.x), 2) + 
                    Math.pow((a.y - b.y), 2) + 
                    Math.pow((a.z - b.z), 2));
    }
    
    static double LengthBoundary(Vertex[] v)
    {
        double length = 0;
        int j = 0;
        for(int i = 0; i < v.length; i++)
        {
            if(i == v.length - 1)
                j = 0;
            else j = i + 1;
            length += LengthEdge(v[i], v[j]);
        } 
        return length;
    }
    
   /* static Point3D CalculationNormal(Model m)
    {
        Point3D N = new Point3D(0, 0, 0);
        for(Point3D face : m.Face)
        {
            Point3D p = new Point3D(m.Vertex.get((int)(face.getY())).getX() - m.Vertex.get((int)(face.getX())).getX(),
            m.Vertex.get((int)(face.getY())).getY() - m.Vertex.get((int)(face.getX())).getY(),
            m.Vertex.get((int)(face.getY())).getZ() - m.Vertex.get((int)(face.getX())).getZ());
            
            Point3D q = new Point3D(m.Vertex.get((int)(face.getZ())).getX() - m.Vertex.get((int)(face.getX())).getX(),
            m.Vertex.get((int)(face.getZ())).getY() - m.Vertex.get((int)(face.getX())).getY(),
            m.Vertex.get((int)(face.getZ())).getZ() - m.Vertex.get((int)(face.getX())).getZ());
        
            Point3D n = new Point3D(p.getX()*q.getX(), p.getY()*q.getY(), p.getZ()*q.getZ());
            n = n.normalize();
            N = new Point3D(N.getX() + n.getX(), N.getY() + n.getY(), N.getZ() + n.getZ());
        }
        return N.normalize();
    }*/
    
    
    public static Hole ReadHole(String FileName, Model model) throws FileNotFoundException, IOException
    {
        Hole hole = new Hole();
        
        File file = new File("C://Users/Fedorov/Desktop/" + FileName);
        if(!file.exists())
            return null;
        
        BufferedReader in = new BufferedReader(new FileReader(file.getAbsoluteFile()));
        String s;
        String[] words;
        int i = 0;
        int n = Integer.parseInt(in.readLine());
        hole.boundary = new Vertex[n];
        hole.close = new Vertex[n];
        hole.num = new int[n];
        
        while((s = in.readLine()) != null)
        {
            words = s.split(" ");
            if(i < n)
            {
                hole.boundary[i] = new Vertex(Double.parseDouble(words[0]), 
                Double.parseDouble(words[1]),
                Double.parseDouble(words[2]));
                hole.num[i] = FindNum(model.vertices, hole.boundary[i]);//Integer.parseInt(words[3]);
            }
            if(i >= n && i < 2*n)
            {
                hole.close[i - n] = new Vertex(Double.parseDouble(words[0]), 
                Double.parseDouble(words[1]),
                Double.parseDouble(words[2]));
            }
            i++;
        }
        in.close();
        return hole;
    }
    
    public static int FindNum(ArrayList<Vertex> vs, Vertex v)
    {
        for(int i = 0; i < vs.size(); i++)
        {
            if(CompareP(vs.get(i), v))
                return i;
        }
        return -1;
    }
    
    public static boolean isHave(ArrayList<Vertex> vs, Vertex v)
    {
        boolean flag = false;
        for(int i = 0; i < vs.size(); i++)
        {
            if(CompareP(vs.get(i), v))
            {
                flag = true;
                break;
            }
        }
        return flag;
    }
    
    /*public static Hole tmpReadHole(Model model, Model hole, Model border)
    {
        Hole h = new Hole();
        int numInHole;
        
        h.boundary = new Vertex[border.vertices.size()];
        h.close = new ArrayList();
        h.num = new int[border.vertices.size()];
        
        for(int ii = 0; ii < border.vertices.size(); ii++)
        {
            h.close.add(new ArrayList());
            h.boundary[ii] = border.vertices.get(ii);
            numInHole = FindNum(hole.vertices, border.vertices.get(ii));
            for(int[] f : hole.faces)
            {
                for(int i = 0; i < 3; i++)
                {
                    if(f[i] == numInHole)
                    {
                        for(int j = 0; j < 3; j++)
                            if(j != i && !isHave(h.close.get(ii), hole.vertices.get(f[j]))
                                    && !isHave(border.vertices, hole.vertices.get(f[j])))
                                h.close.get(ii).add(hole.vertices.get(f[j]));
                    }
                }
            }
            h.num[ii] = FindNum(model.vertices, border.vertices.get(ii));
        }
        return h;
    }*/
    
    static ArrayList<Vertex> sortVertex(Model border, Model hole)
    {
        ArrayList<Vertex> newV = new ArrayList<>();
        int numInHole;
        Vertex curV = border.vertices.get(0);
        newV.add(curV);
        boolean flag = false;
        numInHole = FindNum(hole.vertices, curV);
        do{
            for(int[] f : hole.faces)
            {
                for(int i = 0; i < 3; i++)
                {
                    if(f[i] == numInHole)
                    {
                        for(int ii = 0; ii < 3; ii++)
                    if(isHave(border.vertices, hole.vertices.get(f[ii])) &&
                            !CompareP(curV, hole.vertices.get(f[ii]))
                            && !isHave(newV, hole.vertices.get(f[ii])))
                    {
                        //System.out.println(numInHole);
                        curV = hole.vertices.get(f[ii]);
                        newV.add(curV);
                        numInHole = FindNum(hole.vertices, curV);
                        flag = true;
                    }
                    }
                }
                if(flag)
                {
                    flag = false;
                    break;
                }
            }
        }while(newV.size() != border.vertices.size());
        
        return newV;
    }
    
}
