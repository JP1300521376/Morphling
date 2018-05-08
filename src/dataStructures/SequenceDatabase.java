/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dataStructures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileWriter;
import java.util.*;
import utils.MemoryLogger;
/**
 *
 * @author jiadonglin
 */
public class SequenceDatabase {
    
    
    private long startTime;
    private long endTime;
    
    /** saves all sequences */
    private final List<Sequence> sequences = new ArrayList<Sequence>();
    List<List<SuperItem>> superitemInput = new ArrayList<>();
    Map<Integer, String> seqChromMap = new HashMap<>();
    
//    private double weightRatioThresh = 0.2;                 
    
    /**
     * Load superitems from file and call small indels generated from BWA (pileup)
     * @param superitemFilePath
     * @param svRegionWriter
     * @throws IOException 
     */
    public void loadSequencesFromFile(String superitemFilePath, BufferedWriter svRegionWriter) throws IOException{
//        System.out.println("Start loading mutational database.....");
        System.out.println("\nLoading SuperItem sequence database and formatting...");
        String thisLine;
        startTime = System.currentTimeMillis();
        MemoryLogger.getInstance().reset();                        
        BufferedReader myInput = null;
        
            
        FileInputStream fin = new FileInputStream(new File(superitemFilePath));
        myInput = new BufferedReader(new InputStreamReader(fin));
                    
        // skip header
        myInput.readLine();
        int curSeq = -1;
        while((thisLine = myInput.readLine()) != null){
            String[] tokens = thisLine.split("\t");                  
            String SuperItemType = tokens[0];
            String chrom = tokens[1];
            int ChromIdx = Integer.parseInt(chrom);
            if (!seqChromMap.containsKey(ChromIdx)){
                curSeq = superitemInput.size();
                superitemInput.add(new ArrayList<>());
                seqChromMap.put(superitemInput.size() - 1, chrom);                   
            }else{
                curSeq = ChromIdx;
            }
            double ratio = Double.parseDouble(tokens[7]);
     
            if (SuperItemType.equals("MDM") || SuperItemType.equals("MIM")){
//                outputIndels(svRegionWriter, tokens, ratio);
                continue;
            }
            if (!SuperItemType.contains("ARP")&& (ratio > 0.1)){
                SuperItem superitem = new SuperItem(tokens);
                superitem.setChromName(seqChromMap.get(curSeq));
                superitemInput.get(curSeq).add(superitem);

            }else if (SuperItemType.contains("ARP") && (ratio > 0.1)){   
                SuperItem superitem = new SuperItem(tokens);
                superitem.setChromName(seqChromMap.get(curSeq));
                superitemInput.get(curSeq).add(superitem);
            }               
        }
//        indelWriter.close();
       
        databaseFormatter(superitemInput);
        endTime = System.currentTimeMillis();
        
        printDatabaseStats();
       
    }
    /**
     * This is used to re-format the database, in case of <(AB)(C)> happens.
     * @param superitemSeqs 
     */
    private void databaseFormatter(List<List<SuperItem>> superitemSeqs) {       
        
//        superitemTypeCount(SIsequences);
        for (int i = 0; i < superitemSeqs.size(); i ++){
            Sequence sequence = new Sequence(i);
            List<SuperItem> smallIndels = new ArrayList<>();
            List<SuperItem> SIlist = superitemSeqs.get(i);
            Collections.sort(SIlist);
//            int sequenceSize = 0;
            List<SuperItem> itemset = new ArrayList<>();
            for (SuperItem si : SIlist){
                
                if (si.isSmallIndel()){
                    smallIndels.add(si);
                    continue;
                } 
                if (!smallIndels.isEmpty()){
                    SuperItem indelSuperItem = smallIndels.get(smallIndels.size() - 1);
                    int indelSuperItemStartPos = indelSuperItem.getPos();
                    int indelSuperItemEndPos = indelSuperItem.getSplitAlignPos();
                    
                    if (itemset.isEmpty()){                    
                        if (si.getPos() < indelSuperItemEndPos && si.getPos() > indelSuperItemStartPos){
                            continue;
                        }
                        itemset.add(si);
                    }else{
                        if (si.getPos() < indelSuperItemEndPos && si.getPos() > indelSuperItemStartPos){
                            continue;
                        }
    //                    sequenceSize += 1;
                        sequence.addItemset(itemset);
                        itemset = new ArrayList<>();  
                        itemset.add(si);
                    }
                }else{
                    if (itemset.isEmpty()){                                        
                        itemset.add(si);
                    }else{

    //                    sequenceSize += 1;
                        sequence.addItemset(itemset);
                        itemset = new ArrayList<>();  
                        itemset.add(si);
                    }
                }
                                                                              
            }
            sequence.addItemset(itemset);
            sequence.setSmallIndelSuperitem(smallIndels);
            sequences.add(sequence);
//            System.out.println("Sequence " + sequence.getId() + " processed, " + " size: " + sequenceSize);           
        }
    }
    
    public int size(){
        return sequences.size();
    }
    
    public List<Sequence> getSequences(){
        return sequences;
    }
    
    public Sequence getSequenceByID(int seqID){
        return sequences.get(seqID);
    }
    
    /**
     * Call Indels from SuperItem
     * @param indelOutPath
     * @throws IOException 
     */
    private void outputIndels(BufferedWriter indelWriter, String[] tokens, double ratio) throws IOException{
        StringBuilder sb = new StringBuilder();
        String SuperItemType = tokens[0];
        String chrIdx = tokens[1];
        if (SuperItemType.equals("MDM") && ratio > 0.2){           
            sb.append(chrIdx);
            sb.append("\t");
            sb.append(tokens[3]);
            sb.append("\t");
            sb.append(tokens[4]);
            sb.append("\t");
            sb.append("svType=DEL");
            indelWriter.write(sb.toString());
            indelWriter.newLine();
        }
//        else if (SuperItemType.equals("MIM") && ratio > 0.2){
//            sb.append(chrIdx);
//            sb.append("\t");
//            sb.append(tokens[3]);
//            sb.append("\t");
//            sb.append("-");
//            sb.append("\t");
//            sb.append("svType=INS");
//            indelWriter.write(sb.toString());
//            indelWriter.newLine();
//        }        
    }
    
    private void printDatabaseStats(){
        System.out.println("\n============  Sequence Database STATS ==========");
	System.out.println("Number of sequences : " + sequences.size());
        long size = 0;
        for(Sequence sequence : sequences){
            size += sequence.size();
            
        }
        double meansize = ((float)size) / ((float)sequences.size());
        System.out.println("Average sequence size : " + meansize);
        System.out.println("Time: " + (endTime - startTime) + "ms");
        System.out.println("Memory usage: " + MemoryLogger.getInstance().getMaxMemory());
        System.out.println("================================================\n");
        
    }
}
