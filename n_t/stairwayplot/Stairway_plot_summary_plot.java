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
import simpleGUI.*;
import javax.swing.JFrame;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.plots.*;
import de.erichseifert.gral.plots.axes.*;
import de.erichseifert.gral.plots.lines.*;
import de.erichseifert.gral.ui.InteractivePanel;
import de.erichseifert.gral.util.Insets2D;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.util.*;
import java.io.*;

public class Stairway_plot_summary_plot {
    public static void main(String[] args)throws Exception{
        int beforerow=6;//row before data
        int beforecol=5;//col before the 4 columns for plotting
        FileChooser fc=new FileChooser(".");
        String infile=fc.chooseFile("choose summary file");
        String title=Dialog.getString("Title:");
        String xrange=Dialog.getString("Time (1k year) range (blank for default,format: xmin,xmax");
        String yrange=Dialog.getString("Ne (1k individual) range (blank for default,format: ymin,ymax");
        double xspacing=Dialog.getDouble("X axis spacing",1);
        double yspacing=Dialog.getDouble("Y axis spacing",2);
        int fontsize=Dialog.getInt("Font size", 18);
        double xmin=-1;
        double xmax=-1;
        double ymin=-1;
        double ymax=-1;
        if(!xrange.isEmpty()){
            StringTokenizer t=new StringTokenizer(xrange,", ");
            xmin=Double.parseDouble(t.nextToken());
            xmax=Double.parseDouble(t.nextToken());
        }
        if(!yrange.isEmpty()){
            StringTokenizer t=new StringTokenizer(yrange,", ");
            ymin=Double.parseDouble(t.nextToken());
            ymax=Double.parseDouble(t.nextToken());
        }
        
        //read in data
        DataTable Median=new DataTable(Double.class,Double.class);
        DataTable P25=new DataTable(Double.class,Double.class);
        DataTable P975=new DataTable(Double.class,Double.class);
        BufferedReader in=new BufferedReader(new FileReader(infile));
        for(int i=0;i<beforerow;i++)in.readLine();
        while(in.ready()){
            String line=in.readLine();
            StringTokenizer t=new StringTokenizer(line,"\t");
            for(int i=0;i<beforecol;i++)t.nextToken();
            double year=Double.parseDouble(t.nextToken())/1000;
            double median=Double.parseDouble(t.nextToken())/1000;
            double p25=Double.parseDouble(t.nextToken())/1000;
            double p975=Double.parseDouble(t.nextToken())/1000;
            Median.add(year,median);
            P25.add(year,p25);
            P975.add(year,p975);
        }
        in.close();
        //begin plot
        XYPlot plot = new XYPlot();
        plot.add(0,Median,true);
        plot.add(1,P25,true);
        plot.add(2,P975,true);
        double insetsTop = 20.0,
               insetsLeft = 80.0,
               insetsBottom = 80.0,
               insetsRight = 40.0;
        plot.setInsets(new Insets2D.Double(
                insetsTop, insetsLeft, insetsBottom, insetsRight));
        plot.getTitle().setText(title);
        plot.getTitle().setFont(new Font("SansSerif", Font.BOLD, fontsize));
        plot.getPlotArea().setBorderStroke(new BasicStroke(2f));
        LineRenderer lines1 = new DefaultLineRenderer2D();
        lines1.setColor(new Color(Color.BLACK.getRed(),Color.BLACK.getGreen(),Color.BLACK.getBlue(),255));
        LineRenderer lines2 = new DefaultLineRenderer2D();
        lines2.setColor(new Color(Color.BLUE.getRed(),Color.BLUE.getGreen(),Color.BLUE.getBlue(),255));
        LineRenderer lines5 = new DefaultLineRenderer2D();
        lines5.setColor(new Color(Color.GRAY.getRed(),Color.GRAY.getGreen(),Color.GRAY.getBlue(),255));
        plot.setLineRenderer(Median, lines1);
        plot.setPointRenderer(Median, null);
        plot.setLineRenderer(P25, lines5);
        plot.setPointRenderer(P25, null);
        plot.setLineRenderer(P975, lines5);
        plot.setPointRenderer(P975, null);
        if(!xrange.isEmpty())plot.getAxis(XYPlot.AXIS_X).setRange(xmin, xmax);
        if(!yrange.isEmpty())plot.getAxis(XYPlot.AXIS_Y).setRange(ymin, ymax);
        plot.setAxisRenderer(XYPlot.AXIS_X, new LogarithmicRenderer2D());
        plot.setAxisRenderer(XYPlot.AXIS_Y, new LogarithmicRenderer2D());
        plot.getAxisRenderer(XYPlot.AXIS_X).setLabel("Time (1k year)");
        plot.getAxisRenderer(XYPlot.AXIS_X).setLabelFont(new Font("SansSerif", Font.BOLD, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_Y).setLabel( "Ne (1k individual)");
        plot.getAxisRenderer(XYPlot.AXIS_Y).setLabelFont(new Font("SansSerif", Font.BOLD, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_X).setIntersection(-Double.MAX_VALUE);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setIntersection( -Double.MAX_VALUE);
        plot.getAxisRenderer(XYPlot.AXIS_X).setTickFont(new Font("SansSerif", Font.PLAIN, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTickFont(new Font("SansSerif", Font.PLAIN, fontsize));
        plot.getAxisRenderer(XYPlot.AXIS_X).setTicksAutoSpaced(false);
        plot.getAxisRenderer(XYPlot.AXIS_X).setTickSpacing(xspacing);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTicksAutoSpaced(false);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setTickSpacing(yspacing);
        plot.getAxisRenderer(XYPlot.AXIS_X).setMinorTicksVisible(false);
        plot.getAxisRenderer(XYPlot.AXIS_Y).setMinorTicksVisible(false);
        // Display on screen
        JFrame frame=new JFrame();
        frame.getContentPane().add(new InteractivePanel(plot), BorderLayout.CENTER);
        frame.setDefaultCloseOperation(frame.EXIT_ON_CLOSE);
        frame.setMinimumSize(frame.getContentPane().getMinimumSize());
        frame.setSize(1500,800);
        frame.setVisible(true);
        
    }
}
