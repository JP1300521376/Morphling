/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package options;

/**
 *
 * @author jiadonglin
 */
public class terminalParamLoader {
    
    public int readLen;
    public int fragMean;
    public int fragStd;
    public int cutStd = 3;
    public int clusteringDist;
    public int minMapQ = 10;
    public String chr = null;
    public int givenRegionS;
    public int givenRegionE;
    public String bamFile;
    public String fastaFile;
    public String fastaIndexFile;
    public String superitemOut;
    public String svOut;
    public String abnormalSignalOut = null;
    public String mergedPatternOut = null;
    public String frequentPatternOut = null;
    public boolean hasParamInput = true;
    public terminalParamLoader(){        
        
    }
    
    public void loadParam(String[] args){
        if (args[0].equals("help")){
            printOptionsHelp();
            hasParamInput = false;
        }else{
            loadTerminalInputParams(args);
        }
    }
    private void printOptionsHelp(){
        System.out.println("============ Morphling help info ===============");
        StringBuilder sb = new StringBuilder();
        sb.append("bamFile=   given the path of your BAM file\n");
        sb.append("faFile=   given the path of your reference file\n");     
        sb.append("itemOut=   given the output path of super-item file\n"); 
        sb.append("svOut=   given the output path of discovered SVs\n");
        sb.append("freOut=  given the output path of frequent patterns (optional)\n");
        sb.append("sigOut=  given the output path of abnormal alignments (optional)\n");
        sb.append("patternOut=  given the output path of frequent patterns (optional)\n");
        sb.append("readLen=    given the read length of your BAM file\n");
        sb.append("fragMean=    given the estimated insert size average of your BAM file\n");
        sb.append("fragStd=    given the estimated insert size standard deviation of your BAM file\n");
        sb.append("cutStd=    given the cutoff in unit of the standard deviation (default=3)\n");
        sb.append("maxD=    given the maximum distance to cluster abnormal read pairs (default=readLen)\n");
        sb.append("minQ=    given the minimum mapping quality of read (default=10)\n");
        sb.append("chrom=   given a specific region, process whole chromosome if coordinates are not given. e.g chr1:1000-2000\n");        
        System.out.println(sb.toString());
    }
    
    private void loadTerminalInputParams(String[] args){
        for (int i = 0; i < args.length; i++){                     
            String[] argTokens = args[i].split("=");
            if (argTokens[0].equals("readLen")){
                readLen = Integer.parseInt(argTokens[1]);
                
            }
            if (argTokens[0].equals("fragMean")){
                fragMean = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("fragStd")){
                fragStd = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("cutStd")){
                cutStd = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("maxD")){
                clusteringDist = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("minQ")){
                minMapQ = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("chrom")){
                String givenRegion = argTokens[1];
                if (givenRegion.length() == 1){
                    chr = givenRegion;
                }else{
                    chr = givenRegion.split(":")[0];
                    givenRegionS = Integer.parseInt(givenRegion.split(":")[1].split("-")[0]);
                    givenRegionE = Integer.parseInt(givenRegion.split(":")[1].split("-")[1]);
                }
            }
            if (argTokens[0].equals("bamFile")){
                bamFile = argTokens[1];
            }
            if (argTokens[0].equals("faFile")){
                fastaFile = argTokens[1];
                fastaIndexFile = fastaFile + ".fai";
            }
            if (argTokens[0].equals("itemOut")){
                superitemOut = argTokens[1];
            }
            if (argTokens[0].equals("sigOut")){
                abnormalSignalOut = argTokens[1];
            }
            if (argTokens[0].equals("patternOut")){
                mergedPatternOut = argTokens[1];
            }
            if (argTokens[0].equals("svOut")){
                svOut = argTokens[1];
            }
            if (argTokens[0].equals("freOut")){
                frequentPatternOut = argTokens[1];
            }
        }
    }
}
