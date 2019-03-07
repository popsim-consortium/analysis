/**
 *
 * Copyright (c) @author Xiaoming Liu, Ph.D.
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
public class ms_output_to_stairway_plot_input {
    public static void main(String[] args)throws Exception{
        FileChooser fc = new  FileChooser(".");
        String infile =fc.chooseFile("ms output file:");
        String popid=Dialog.getString("population id:");
        String L=Dialog.getString("seq length simulated:", "10000000");
        int nrep=Dialog.getInt("number of replications performed", 200);
        String outfile =fc.chooseFile("stairway plot input file stem");
        double corr=0.99;//Dialog.getDouble("autocorr:", 0.99);
        
        BufferedReader in=new BufferedReader(new FileReader(infile));
        boolean ready=true;
        boolean begin=false;
        String line=in.readLine();
        StringTokenizer t=new StringTokenizer(line);
        t.nextToken();
        int nseq=Integer.parseInt(t.nextToken());
        
        int[] xi=new int[nseq];
        Arrays.fill(xi,0);
        while(ready){
            line="";
            try{
                    line = in.readLine();
            }catch(java.io.IOException e){
                    //System.out.println(e);
            }
            if(line==null){

                ready=false;
                line="";
                continue;
            }
            if(line.equalsIgnoreCase("//")) {
                PrintWriter out=new PrintWriter(new FileWriter(outfile+nrep),true);
                nrep--;
                begin=true;
                int nsite=0;
                String time="time:";
                String segsites="segsites:";
                while(begin){
                    line=in.readLine();
                    if(line.isEmpty())continue;
                    t=new StringTokenizer(line);
                    String tmp=t.nextToken();
                    if(tmp.equalsIgnoreCase("time:")) time=line;
                    else if(tmp.equalsIgnoreCase("segsites:")) {
                        segsites=line;
                        nsite=Integer.parseInt(t.nextToken());
                        
                    }
                    else if(tmp.equalsIgnoreCase("positions:")){
                        out.println(popid+"\t"+nseq+"\t"+L+"\t"+1+"\t"+(nseq-1));
                        Arrays.fill(xi,0);
                        String[] seq=new String[nseq];
                        for(int i=0;i<nseq;i++) seq[i]=in.readLine();
                        int count=0;
                        for(int i=0;i<nsite;i++){
                            count=0;
                            for(int j=0;j<nseq;j++) if(seq[j].charAt(i)=='1')count++;
                            xi[count]++;
                        }
                        count=0;
                        for(int i=1;i<=nseq-1;i++) {
                            count+=xi[i];
                            if(i<nseq-1)out.print(xi[i]+"\t");
                            else out.print(xi[i]);
                        }
                        out.println();
                        if(xi[0]>0)System.out.println("xi[0]="+xi[0]);
                        if(count!=nsite) System.out.println("count="+count+" nsite="+nsite);
                        
                        in.readLine();
                        begin=false;
                    }
                }
                out.close();
            }

        }
        in.close();
        
    }
}
