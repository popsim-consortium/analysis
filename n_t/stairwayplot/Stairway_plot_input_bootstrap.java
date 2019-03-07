/**
 *
 * Copyright (c) 2014 @author Xiaoming Liu, Ph.D.
 * Assistant Professor,
 * Human Genetics Center,
 * School of Public Health,
 * The University of Texas Health Science Center at Houston
 * 
 * This source code is distributed under the RECEX SHARED SOURCE LICENSE
 * 
 * You are free to download, copy, compile, study, and refer to the source code for any personal use of yours.
 * You are free to make any modifications to the source covered by this license.
 * You may NOT under any circumstance copy, redistribute and/or republish the source or a work based on it (which
 * includes binary or object code compiled from it) in part or whole.
 * If you intend to incorporate the source code, in part or whole, into any free or proprietary program, you need to explicitly
 * write to the original author(s) to ask for permission.
 * The source code licensed under this license is shared "as is".
 * 
 * You shall already get a copy of the license, if not, you can obtain a copy at 
 * https://raw.github.com/Recex/Licenses/master/SharedSourceLicense/LICENSE.txt
 */
import java.util.*;
import java.io.*;
import simpleGUI.*;
import cern.jet.random.Binomial;
public class Stairway_plot_input_bootstrap {
    
    public static void main(String[] args)throws Exception{
        FileChooser fc = new  FileChooser(".");
        String infile =fc.chooseFile("input file:");
        int nbootstrap=Dialog.getInt("number of bootstrap samples", 199);
        BufferedReader in=new BufferedReader(new FileReader(infile));
        while(in.ready()){
            String paraline=in.readLine();
            if(paraline.isEmpty())continue;
            StringTokenizer t=new StringTokenizer(paraline);
            String popid=t.nextToken();
            int nseq=Integer.parseInt(t.nextToken());
            int L=Integer.parseInt(t.nextToken());
            String temp=in.readLine();
            
            int totale=0;
            int n=nseq;
            
            int[] exi=new int[n];
            t=new StringTokenizer(temp);
            for(int i=1;i<=n-1;i++){
                
                exi[i]=Integer.parseInt(t.nextToken());
                totale+=exi[i];
            }
            if(L-totale>0){
                exi[0]=(int)(L-totale);
            }
            else {
                System.out.println("1-sum(theta)>1");
                exi[0]=0;
            }
            //the first file is the original file
            PrintWriter out=new PrintWriter(new FileWriter(infile+".bootstrap"+1),true);
            out.println(paraline);
            for(int i=1;i<nseq;i++) out.print("\t"+exi[i]);
            out.println();
            out.close();
            //create bootstrap files
            for(int ii=2;ii<=nbootstrap+1;ii++){
                out=new PrintWriter(new FileWriter(infile+".bootstrap"+ii),true);
                int[] tmp=Multinomial_rn(exi,L);
                out.println(paraline);
                for(int j=1;j<nseq;j++) out.print("\t"+tmp[j]);
                out.println();
                out.close();
            }
        }
        in.close();
        
    }
    private static int[] Multinomial_rn(int[] Count, int N){
        int n=Count.length;
        double[] p=new double[n];
        for(int i=0;i<n;i++) p[i]=Count[i];
        int[] dist=new int[n];
        Arrays.fill(dist,0);
        double[] cp=new double[n];
        cp[0]=p[0];
        for(int i=1;i<n;i++)cp[i]=cp[i-1]+p[i];
        int remain=N;
        int category=n-1;
        while(remain>0){
            int deviate = Binomial.staticNextInt(remain,p[category]/cp[category]);
            dist[category]=deviate;
            remain=remain-deviate;
            if(remain<0)System.out.println("remain="+remain);
            category--;
            if(category==0){
                dist[0]=remain;
                remain=0;
            }
        }
        return dist;
    }

}
