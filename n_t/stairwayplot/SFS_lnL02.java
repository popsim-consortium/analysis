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
import swarmops.Problem;
import swarmops.Tools;

public class SFS_lnL02 extends Problem{
    int n;
    long L;
    long K=0;
    double[] logfac;
    double[] logexi;
    long[] xi;
    boolean initialized=false;
    double maxfit;
    double[] lowerBound;
    double[] upperBound;
    String[] parameterName;
    int[][] group;
    boolean[] obs;
    int nobs;
    boolean setdimpenalty=false;
    double dimpenalty=0;
    boolean setautocorr=false;
    double pautocorr=1;
    public SFS_lnL02(int N){
        super();
        n=N;
        logfac= new double[n+1];
        logfac[0]=0;
	for(int i=1;i<=n;i++){
            logfac[i]=logfac[i-1]+Math.log(1.0*i);
        }
        logexi=new double[n];
        xi=new long[n];
        group=new int[n-1][1];
        for(int i=0;i<n-1;i++)group[i][0]=i+2;//default group: each theta(i) as a group
        obs=new boolean[n];
        for(int i=0;i<n;i++) obs[i]=true;//default observed xi: all xi
        nobs=n;
    }
    public double[] getLogExi(double[] Theta){
        if(!initialized) {
            getInitialize(Theta);
        }    
        return logexi;
    }
    
    public boolean getInitialize(double[] Theta){
        //define group
        double currenttheta=Theta[2];
        ArrayList list=new ArrayList();
        list.add(new Integer(2));
        for(int i=3;i<=n;i++){
            if(Math.abs(Theta[i]-currenttheta)>1e-15){
                currenttheta=Theta[i];
                list.add(new Integer(i-1));
                list.add(new Integer(i));
            }
        }
        list.add(new Integer(n-1));
        Integer[] group0=(Integer[])list.toArray(new Integer[list.size()]);
        //calculate
        double totale=0;
        for(int i=1;i<=n-1;i++){
            list.clear();
            for(int ii=0;ii<group0.length-1;ii+=2){
                int begin=group0[ii].intValue();
                int end=group0[ii+1].intValue();
                if(begin>n-i+1)begin=n-i+1;
                if(end>n-i+1)end=n-i+1;
                list.add(new Integer(begin));
                list.add(new Integer(end));
                if(end==n-i+1)break;
            }
            Integer[] group=(Integer[])list.toArray(new Integer[list.size()]);
            double exii=0;
            if(group.length==2){//single theta value
                exii=1.0/i*Theta[2];
            }
            else if(group.length==4){//two theta values
                int j=group[1].intValue();
                exii=(1-Math.exp(logfac[n-j]+logfac[n-i-1]-logfac[n-i-j]-logfac[n-1]))/i*Theta[j];
                j=group[2].intValue();
                exii+=Math.exp(logfac[n-j+1]+logfac[n-i-1]-logfac[n-i-j+1]-logfac[n-1])/i*Theta[j];
            }
            else{
                int j=group[1].intValue();
                exii=(1-Math.exp(logfac[n-j]+logfac[n-i-1]-logfac[n-i-j]-logfac[n-1]))/i*Theta[j];
                for(int ii=2;ii<=group.length-4;ii+=2){
                    j=group[ii].intValue();
                    int l=group[ii+1].intValue();
                    exii+=(Math.exp(logfac[n-i-1]+logfac[n-j]-logfac[n-i-j+1]-logfac[n-1])*(n-j+1)-Math.exp(logfac[n-i-1]+logfac[n-l-1]-logfac[n-i-l]-logfac[n-1])*(n-l))/i*Theta[j];
                }
                j=group[group.length-2].intValue();
                exii+=Math.exp(logfac[n-j+1]+logfac[n-i-1]-logfac[n-i-j+1]-logfac[n-1])/i*Theta[j];
            }
            totale+=exii;
            if(exii>0) logexi[i]=Math.log(exii);
            else logexi[i]=Math.log(Double.MIN_VALUE);
        }
        if(1-totale>0){
            logexi[0]=Math.log(1-totale);
            initialized=true;
            return true;
        }
        else {
            logexi[0]=Math.log(Double.MIN_VALUE);
            return false;
        }
    }
    public void setData(long[] Xi){
        K=0;
        for(int i=0;i<=n-1;i++){
            xi[i]=Xi[i];
            K+=xi[i];
        }
        double totall=0;
        L=0;
        for(int i=0;i<=n-1;i++){
            L+=xi[i];
        }
        double eother=xi[0]*1.0/L;
        double oother=xi[0];
        for(int i=1;i<n;i++){
            if(xi[i]>0&&obs[i])totall+=Math.log(xi[i]*1.0/L)*xi[i]-gammln(xi[i]+1.0);
            else{
                eother+=xi[i]*1.0/L;
                oother+=xi[i];
            }
        }
        totall+=Math.log(eother)*oother-gammln(oother+1.0);
        totall+=gammln(L+1.0);//add constant
        maxfit=totall;//no panelty
        
    }
    
    public double getLogLikelihood(double[] Theta){
        getInitialize(Theta);
        double totall=0;
        double eother=Math.exp(logexi[0]);
        double oother=xi[0];
        for(int i=1;i<n;i++){
            if(obs[i]) totall+=logexi[i]*xi[i]-gammln(xi[i]+1.0);
            else{
                eother+=Math.exp(logexi[i]);
                oother+=xi[i];
            }
        }
        totall+=Math.log(eother)*oother-gammln(oother+1.0);
        totall+=gammln(L+1.0);//add constant
        double adjustlogL=totall;//no panelty
        if(setdimpenalty)adjustlogL=totall-dimpenalty;
        return adjustlogL;
    }
    public void setDimPenalty(double p){
        dimpenalty=p;
        setdimpenalty=true;
    }
    public void setAutoCorr(double p){
        setautocorr=true;
        pautocorr=p;
    }
    public void setThetaGroup(int[][] Group){
        group=(int[][])Group.clone();
        for(int i=0;i<Group.length;i++)group[i]=(int[])Group[i].clone();
        
    }
    /**
     * set whether xi(i) is used for analysis
     * @param Obs 
     */
    public void setObsXi(boolean[] Obs){
        obs=(boolean[])Obs.clone();
        nobs=1;
        for(int i=1;i<n;i++) if(obs[i])nobs++;
    }
    public String getName() {
            return "SFS_lnL2";
    }
    public int getDimensionality() {
            return group.length;
    }

    public double[] getLowerBound() {
        lowerBound = new double[getDimensionality()];
        Arrays.fill(lowerBound,0);
        return lowerBound;
    }

    public double[] getUpperBound() {
        upperBound = new double[getDimensionality()];
        Arrays.fill(upperBound,0.2);
        return upperBound;
    }
    @Override
    public double[] getLowerInit() {
            return getLowerBound();
    }

    @Override
    public double[] getUpperInit() {
            return getUpperBound();
    }
    public double getMinFitness() {
        return -maxfit;
    }
    @Override
    public double getAcceptableFitness() {
        return -maxfit;
    }
    @Override
    public String[] getParameterName() {
        parameterName =new String[getDimensionality()];
        for(int i=0;i<n-1;i++) parameterName[i]="theta"+(i+2);
        return parameterName;
    }
    @Override
    public double fitness(double[] x) {
        assert x != null && x.length == getDimensionality();
        double[] Theta=new double[n+1];
        for(int i=0;i<group.length;i++){
            for(int j=0;j<group[i].length;j++) Theta[group[i][j]]=x[i];
        }
        return -getLogLikelihood(Theta);
    }
    @Override
    public boolean enforceConstraints(double[] x) {
        // Enforce boundaries.
        Tools.bound(x, getLowerBound(), getUpperBound());
        // Return feasibility.
        return isFeasible(x);
    }
    @Override
    public boolean isFeasible(double[] x) {//isFeasible come first before geting fitness
        assert x != null && x.length == getDimensionality();
        boolean feasible=true;
        double[] Theta=new double[n+1];
        for(int i=0;i<group.length;i++){
            for(int j=0;j<group[i].length;j++) Theta[group[i][j]]=x[i];
        }
        for(int i=2;i<=n;i++) if(Theta[i]==0) return false;
        if(setautocorr){
            double p2=(1-pautocorr)/2;
            for(int i=3;i<=n;i++){//assume theta[i] draw from an exponential distribution with mean Theta[i-1], from Bayesian skyline plot
                double lambda=1/Theta[i-1];
                double cdf=1-Math.exp(-lambda*Theta[i]);
                if(cdf<p2||cdf>1-p2) return false;
            }
            for(int i=2;i<n;i++){//assume theta[i] draw from an exponential distribution with mean Theta[i+1], from Bayesian skyline plot
                double lambda=1/Theta[i+1];
                double cdf=1-Math.exp(-lambda*Theta[i]);
                if(cdf<p2||cdf>1-p2) return false;
            }
        }
        double k=0;
        for(int i=2;i<=n;i++)k+=Theta[i]/(i-1);
        if(k>1)feasible=false;
        return feasible;
    }
    /**
   * Compute LnGamma(x)
   */
    private static double gammln(double xx)
    {
        int j;
        double temp;
        double cof[] = new double[7];
        double stp, half, one, fpf, x, tmp, ser;
        cof[1] = 76.18009172947146;
        cof[2] = -86.50532032941677;
        cof[3] = 24.01409824083091;
        cof[4] = -1.231739572450155;
        cof[5] = 0.001208650973866179;
        cof[6] = -0.000005395239384953;
        stp = 2.5066282746310005;
        half = 0.5;
        one = 1.0;
        fpf = 5.5;
        x = xx ;
        tmp = x + fpf;
        tmp = (x + half) * Math.log(tmp) - tmp;
        ser = 1.000000000190015;
        for (j = 1; j <= 6; j++)
        {
            x = x + one;
            ser = ser + cof[j] / x;
        }
        temp = tmp + Math.log(stp * ser/xx);
        return temp;
    }

}
