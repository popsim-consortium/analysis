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
import swarmops.*;
import swarmops.optimizers.*;

public class Stairway_plot_theta_estimation02 {
     public static void main(String[] args)throws Exception{
        boolean verbose=false;
        boolean fitnesstrail=true;
        int numRuns = 1;
        int dimFactor = 5000;//dimfactor;
        int totalRuns=0;
        if(args.length!=3){
            System.out.println("Usage: java Stairway_plot_theta_estimation2 <input file name> <numRuns> <dimFactor>");
            System.exit(0);
        }
        
        String infile =args[0];
        numRuns=Integer.parseInt(args[1]);
        dimFactor=Integer.parseInt(args[2]);
        PrintWriter out=new PrintWriter(new FileWriter(infile+".addTheta"),true);
        BufferedReader in=new BufferedReader(new FileReader(infile));
        int nline=0;
        while(in.ready()){
            in.readLine();
            nline++;
        }
        in.close();
        in=new BufferedReader(new FileReader(infile));
        int nsfs=nline/2;
        for(int ii=0;ii<nsfs;ii++){
            String line=in.readLine();
            out.println(line);
            StringTokenizer t=new StringTokenizer(line,"\t");
            String popid=t.nextToken();
            int nseq=Integer.parseInt(t.nextToken());
            double pautocorr=0.99;//Double.parseDouble(t.nextToken());
            long L=Long.parseLong(t.nextToken());
            int xibegin=Integer.parseInt(t.nextToken());//begin of observed Xi, for example 2 for xi(2)
            int xiend=Integer.parseInt(t.nextToken());//end of observed Xi, for example n-2 for xi(n-2)
            
            line=in.readLine();
            out.println(line);
            t=new StringTokenizer(line,"\t");
            long[] c=new long[nseq];//count xi(i)
            long total=0;
            for(int i=1;i<=nseq-1;i++) {
                c[i]=Integer.parseInt(t.nextToken());
                total+=c[i];
            }
            c[0]=L-total;
            int nx=nseq;
            SFS_lnL02 problem = new SFS_lnL02(nx);
            boolean[] obs=new boolean[nx];
            Arrays.fill(obs,false);
            for(int i=xibegin;i<=xiend;i++) obs[i]=true;
            obs[0]=true;//not necessory
            problem.setObsXi(obs);
            Optimizer optimizer = new DE(problem);
            double[] parameters = optimizer.getDefaultParameters();
            boolean keepgoing=true;
            boolean[] splitbefore=new boolean[nx+1];//split before theta i
            Arrays.fill(splitbefore,false);
            splitbefore[2]=true;// dim=1, all thetas are equal, no split
            int dim=1;
            int[][] group=new int[dim][];
            int current=0;
            int count=1;
            for(int i=3;i<=nx;i++){
                if(!splitbefore[i])count++;
                else {
                    group[current]=new int[count];
                    count=1;
                    current++;
                }
            }
            group[current]=new int[count];
            current=0;
            group[0][0]=2;
            count=1;
            for(int i=3;i<=nx;i++){
                if(!splitbefore[i]){
                    group[current][count]=i;
                    count++;
                }
                else {
                    current++;
                    group[current][0]=i;
                    count=1;
                }
            }
            
            int numIterations = dimFactor * dim;
            problem.setThetaGroup(group);
            problem.setData(c);
            if(pautocorr<1)problem.setAutoCorr(pautocorr);
            Globals.random = new swarmops.random.MersenneTwister();
            problem.maxIterations = numIterations;
            double bestfit=Double.MAX_VALUE;
            double[] bestest=new double[dim];
            int nfeasible=0;
            for(int i=0;i<numRuns;i++){
                Result result=optimizer.optimize(parameters);
                totalRuns++;
                if(result.feasible){
                    if(result.fitness<bestfit){
                        bestfit=result.fitness;
                        bestest=new double[dim];
                        System.arraycopy(result.parameters, 0, bestest, 0, dim);
                        if(verbose){
                            System.out.println(i+" "+result.iterations+" "+result.fitness+" "+result.feasible+" "+result.parameters.length);
                            for(int j=0;j<dim;j++)System.out.print(bestest[j]+" ");
                            System.out.println();
                        }
                    }
                    nfeasible++;
                }
                else i--;
            }
            boolean[] currentsplit=(boolean[])splitbefore.clone();
            boolean[] bestsplit=(boolean[])splitbefore.clone();
            if(fitnesstrail){
                out.println("fitness trail:");
                out.print("dim:\t"+dim+"\tbestfit:\t"+bestfit);
                for(int i=2;i<=nx;i++){
                    if(bestsplit[i])out.print("\t"+i);
                    else out.print(","+i);
                }
                for(int i=0;i<dim;i++) out.print("\t"+bestest[i]);
                out.println();
            }
            int theta_begin=3;
            while(keepgoing){
                double previous_bestfit=bestfit;
                dim++;
                numIterations = dimFactor * dim;
                problem.maxIterations = numIterations;
                keepgoing=false;
                currentsplit=(boolean[])bestsplit.clone();
                for(int j=theta_begin;j<=nx;j++){
                    splitbefore=(boolean[])currentsplit.clone();
                    if(!splitbefore[j])splitbefore[j]=true;
                    else continue;
                    group=new int[dim][];
                    current=0;
                    count=1;
                    for(int i=3;i<=nx;i++){
                        if(!splitbefore[i])count++;
                        else {
                            group[current]=new int[count];
                            count=1;
                            current++;
                        }
                    }
                    group[current]=new int[count];
                    current=0;
                    group[0][0]=2;
                    count=1;
                    for(int i=3;i<=nx;i++){
                        if(!splitbefore[i]){
                            group[current][count]=i;
                            count++;
                        }
                        else {
                            current++;
                            group[current][0]=i;
                            count=1;
                        }
                    }
                    problem.setThetaGroup(group);
                    nfeasible=0;
                    for(int i=0;i<numRuns;i++){
                        Result result=optimizer.optimize(parameters);
                        totalRuns++;
                        if(result.feasible){
                            if(verbose)System.out.println(dim+" "+j+" "+result.fitness);
                            if(result.fitness<bestfit){
                                boolean theta0=false;
                                for(int k=0;k<dim;k++) if(result.parameters[k]==0) theta0=true;
                                if(!theta0){//avoid any theta=0
                                    bestfit=result.fitness;
                                    bestest=new double[dim];
                                    System.arraycopy(result.parameters, 0, bestest, 0, dim);
                                    bestsplit=(boolean[])splitbefore.clone();
                                    if(verbose){
                                        System.out.println(i+" "+result.iterations+" "+result.fitness+" "+result.feasible+" "+result.parameters.length);
                                        for(int jj=0;jj<dim;jj++)System.out.print(bestest[jj]+" ");
                                        System.out.println();
                                    }
                                }
                            }
                            nfeasible++;
                        }
                        else i--;
                    }
                }
                if(bestfit+3.841/2<previous_bestfit)keepgoing=true;//Chi square test with 1 df, alpha=0.05
                if(fitnesstrail&&keepgoing){
                    out.print("dim:\t"+dim+"\tbestfit:\t"+bestfit);
                    for(int i=2;i<=nx;i++){
                        if(bestsplit[i])out.print("\t"+i);
                        else out.print(","+i);
                    }
                    for(int i=0;i<dim;i++) out.print("\t"+bestest[i]);
                    out.println();
                }
            }
            dim=bestest.length;
            double[] solution=new double[dim];
            out.println("final model:\t"+bestfit);
            for(int i=0;i<dim;i++)solution[i]=bestest[i]*L;
            for(int i=2;i<=nx;i++){
                if(bestsplit[i]){
                    if(i==2)out.print(i);
                    else out.print("\t"+i);
                }
                else out.print(","+i);
            }
            out.println();
            for(int i=0;i<dim-1;i++) out.print(solution[i]+"\t");
            out.println(solution[dim-1]);

        }
        
        in.close();
        out.close();
    }

}
