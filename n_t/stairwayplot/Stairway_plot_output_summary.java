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
public class Stairway_plot_output_summary {
    public static void main(String[] args)throws Exception{
        FilesChooser mfc = new  FilesChooser(".");
        String[] infiles =mfc.chooseFiles("Stairway plot output (.addTheta) files:");
        Arrays.sort(infiles);
        FileChooser fc = new  FileChooser(mfc.getDir());
        double mu=Dialog.getDouble("assumed mutation rate per site per generation:", 1.2e-8);
        double year_per_generation=Dialog.getDouble("assumed generation time (in years)", 24);
        String outfile =fc.chooseFile("output summary file:");
        PrintWriter out=new PrintWriter(new FileWriter(outfile),true);
        BufferedReader in=new BufferedReader(new FileReader(infiles[0]));
        String line=in.readLine();
        out.println(line);

        StringTokenizer t=new StringTokenizer(line,"\t");
        int nt=t.countTokens();

        String[] mspara=new String[nt];
        for(int i=0;i<nt;i++) mspara[i]=t.nextToken();
        int nseq=Integer.parseInt(mspara[1]);
        int nrep=infiles.length;
        long L=Long.parseLong(mspara[2]);
        
        System.out.println("nrep="+nrep);
        out.println("nrep="+nrep);
        double[][] allest=new double[nseq-1][nrep];
        out.println("final fitness:");   
        for(int ii=0;ii<nrep-1;ii++)out.print(infiles[ii]+"\t");
        out.println(infiles[nrep-1]);
        for(int ii=0;ii<nrep;ii++){
            in=new BufferedReader(new FileReader(infiles[ii]));
            in.readLine();
            in.readLine();
            System.out.println((ii+1)+"\t"+infiles[ii]);
            line=in.readLine();
            double fitness=0;
            while(line.indexOf("final model:")==-1)line=in.readLine();
            t=new StringTokenizer(line,"\t");
            t.nextToken();
            fitness=Double.parseDouble(t.nextToken());
            if(ii==nrep-1) out.print(fitness);
            else out.print(fitness+"\t");
            line=in.readLine();
            t=new StringTokenizer(line);
            int dim=t.countTokens();
            int[][] group=new int[dim][];
            for(int i=0;i<dim;i++){
                String tmp=t.nextToken();
                StringTokenizer t2=new StringTokenizer(tmp,",");
                int nxi=t2.countTokens();
                group[i]=new int[nxi];
                for(int j=0;j<nxi;j++)group[i][j]=Integer.parseInt(t2.nextToken());
                
            }
            
            line=in.readLine();
            t=new StringTokenizer(line);
            for(int i=0;i<dim;i++){
                double theta=Double.parseDouble(t.nextToken());
                for(int j=0;j<group[i].length;j++) {
                    allest[group[i][j]-2][ii]=theta;
                    //allest0[group[i][j]-2][ii]=theta;
                }
            }
            
            in.close();
        }
        out.println();
        int dim=nseq-1;
        double[] time=new double[dim];
        for(int i=0;i<dim;i++) {
            Arrays.sort(allest[i]);
            time[i]=allest[i][nrep/2]/L/(i+2)/(i+1);//median time
        }
        for(int i=dim-2;i>=0;i--) time[i]=time[i]+time[i+1];//cummulative time
        
        out.println("mutation_per_site\ttheta\ttheta_per_site_median\ttheta_per_site_2.5%\ttheta_per_site_97.5%\tyear\tNe_median\tNe_2.5%\tNe_97.5%");
        for(int i=dim-1;i>=0;i--){
            if(i==dim-1){
                out.print(0+"\t"+(i+2)+"\t"+allest[i][nrep/2]/L+"\t"+allest[i][(int)(nrep*0.025)]/L+"\t"+allest[i][(int)(nrep*0.975)]/L+"\t"+0/mu*year_per_generation+"\t"+allest[i][nrep/2]/L/mu/4+"\t"+allest[i][(int)(nrep*0.025)]/L/mu/4+"\t"+allest[i][(int)(nrep*0.975)]/L/mu/4);
                out.println();
                out.print(time[i]+"\t"+(i+2)+"\t"+allest[i][nrep/2]/L+"\t"+allest[i][(int)(nrep*0.025)]/L+"\t"+allest[i][(int)(nrep*0.975)]/L+"\t"+time[i]/mu*year_per_generation+"\t"+allest[i][nrep/2]/L/mu/4+"\t"+allest[i][(int)(nrep*0.025)]/L/mu/4+"\t"+allest[i][(int)(nrep*0.975)]/L/mu/4);
                out.println();
            }
            else{
                out.print(time[i+1]+"\t"+(i+2)+"\t"+allest[i][nrep/2]/L+"\t"+allest[i][(int)(nrep*0.025)]/L+"\t"+allest[i][(int)(nrep*0.975)]/L+"\t"+time[i+1]/mu*year_per_generation+"\t"+allest[i][nrep/2]/L/mu/4+"\t"+allest[i][(int)(nrep*0.025)]/L/mu/4+"\t"+allest[i][(int)(nrep*0.975)]/L/mu/4);
                out.println();
                out.print(time[i]+"\t"+(i+2)+"\t"+allest[i][nrep/2]/L+"\t"+allest[i][(int)(nrep*0.025)]/L+"\t"+allest[i][(int)(nrep*0.975)]/L+"\t"+time[i]/mu*year_per_generation+"\t"+allest[i][nrep/2]/L/mu/4+"\t"+allest[i][(int)(nrep*0.025)]/L/mu/4+"\t"+allest[i][(int)(nrep*0.975)]/L/mu/4);
                out.println();
            }
        }
        
        out.close();
    }
}
