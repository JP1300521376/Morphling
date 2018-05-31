/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mophlingv2;

import contiguousfspm.ContiguousFSPM;
import dataStructures.SequenceDatabase;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import superitemGenerator.SignalReader;

/**
 *
 * @author jiadonglin
 */
public class MophlingV2 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException{
        // TODO code application logic here
        String workingDir = "/Users/jiadonglin/SV_data/NA19238/";
//        String workingDir = "/Users/jiadonglin/SV_data/tumors/";
//        String workingDir = "/Users/jiadonglin/SV_data/NA12878/";

        String chr = null;      
        int regionS = 0;
        int regionE = 0;
        
        String bamFile = workingDir + "NA19238.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.bam";
//        String bamFile = workingDir + "CPD.bam";
        String fastaIndexFile = "/Users/jiadonglin/SV_data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai";
        String fastaFile = "/Users/jiadonglin/SV_data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa";
        
//        String bamFile = workingDir + "NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam";
//        String fastaIndexFile = "/Users/jiadonglin/SV_data/1K_Project_Ref/hs37d5.fa.fai";
        String superitemOut = workingDir + "morphlingv3/wgs.all.superitems.v3.txt";
        String mergedPatternOut = workingDir + "morphlingv3/wgs.merged.patterns.debug.txt";
//        String frequentPatternOut = workingDir + "morphling/chr1.frequent.patterns.txt";       
        String svRegionOut = workingDir + "morphlingv3/wgs.svRegion.sup30.mapq20.debug.out";
        String centromeres = "/Users/jiadonglin/SV_data/ref_genome/centromeres.txt";
//        int fragMean = 425;
//        int readLen = 250;
//        int fragStd = 147; 

        int fragMean = 564;
        int readLen = 126;
        int fragStd = 160; 
        
        String signalOut = null;        
        
        int cutStd = 3;
        int clusteringDist = readLen;

        int minMapQ = 10;                
        
        System.out.println("fragMean: " + fragMean + " fragStd: " + fragStd + " readLen: " + readLen + " clustD: " + clusteringDist);

//        SignalReader myReader = new SignalReader(fragMean, fragStd, cutStd, readLen, clusteringDist, minMapQ);
//        myReader.doWork(bamFile, fastaIndexFile, chr, regionS, regionE, superitemOut, signalOut); 
       
                
        SequenceDatabase sequenceDatabase = new SequenceDatabase(); 
        int minSup = 30;
        int patternMaxRegionSpan = fragMean;
        
        BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(svRegionOut));
        
        sequenceDatabase.loadSequencesFromFile(superitemOut, svRegionWriter);
        
        ContiguousFSPM algoContiguousFSPM = new ContiguousFSPM(minSup, patternMaxRegionSpan);
        algoContiguousFSPM.setCenteromeres(centromeres);
        algoContiguousFSPM.runAlgorithm(sequenceDatabase, null, mergedPatternOut, svRegionWriter, fastaFile);
        algoContiguousFSPM.printAlgoStatistics();
    }
    
}
