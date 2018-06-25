/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package superitemGenerator;
import htsjdk.samtools.*;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;

import java.io.BufferedWriter;
import java.io.FileWriter;

import Channels.*;
import dataStructures.SuperItem;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
/**
 *
 * @author jiadonglin
 */
public class SignalReader {
    
    private long startTime;
    private long endTime;
    
    private BufferedWriter writer = null;
    

    private String[] chromInfoMap;
    private final Map<Integer, String> idxToChromMap = new HashMap<>();
    
    private int[] chromLength;
    
        
    private final int isizeUpper;
    private final int isizeLower;
    private final int readLen;
    private final int fragMean;
    
    private final int minMapQ;
    private int clusteringDist; 
    
//    private int superitemMaxSpanRange;
    
    
    private Map<String, List<SAMRecord>> rpTracker;
    private ChannelParser breakChannel;
    private ChannelParser isizeLargeChannel;
    private ChannelParser isizeSmallChannel;
    private ChannelParser oemChannel;
    private ChannelParser oriChannel;
    
    private int superitemCount;
    // keep normal read per base every 1Mb length.
    private int readDepthContainerBuffer = 1000000;
    private int[] readDepthContainer = new int[readDepthContainerBuffer];
    
    private int numOfARPs;
    private int numOfRPs;
    
    
    
    public SignalReader(int fragMean, int fragStd, int cutStd, int readLen, int maxDist, int minMapQ){
        isizeUpper = fragMean + cutStd * fragStd;
        isizeLower = fragMean - cutStd * fragStd;
        
        this.readLen = readLen;
        this.fragMean = fragMean;
        this.minMapQ = minMapQ;
        this.clusteringDist = maxDist;
        if (maxDist == -1){
            this.clusteringDist = fragMean - 2 * readLen;
        }
//        superitemMaxSpanRange = 2 * fragMean;
                
    }
    /**
     * 
     * @param bamFile the alignment file
     * @param fastaIndexFile the reference index file
     * @param chrom a user specified chromosome (optional)
     * @param chromStart start of a genome region (optional)
     * @param chromEnd end of a genome region (optional)
     * @param superitemOutPath output path of the created superitems
     * @param abnormalSigOut output path of the abnormal signals (optional)
     * @throws IOException 
     */
    public void doWork(String bamFile, String fastaIndexFile, String chrom, int chromStart, int chromEnd, String superitemOutPath, String abnormalSigOut) throws IOException{
        startTime = System.currentTimeMillis();
        
        // the channel constructor need a name, it can be whatever you like, just used for naming some output file.
        breakChannel = new ChannelParser(clusteringDist, "break", abnormalSigOut);
        isizeLargeChannel = new ChannelParser(clusteringDist, "isize_large", abnormalSigOut);
        isizeSmallChannel = new ChannelParser(clusteringDist, "isize_small", abnormalSigOut);
        oemChannel = new ChannelParser(clusteringDist, "oem", abnormalSigOut);
        oriChannel = new ChannelParser(clusteringDist, "ori", abnormalSigOut);        
        
        if (!superitemOutPath.isEmpty()){
            createSuperItemWriter(superitemOutPath);            
        }
        extractSignalsFromBAM(bamFile, fastaIndexFile, chrom, chromStart, chromEnd);
                
        endTime = System.currentTimeMillis();
        writer.close();
        printSuperItemGeneratorStats();
            
    }
    
    private void extractSignalsFromBAM(String bamFile, String fastaIndexFile, String chrom, int regionStart, int regionEnd) throws IOException{        
               
        rpTracker = new HashMap<>();


        int windowStart = 0;
        int windowEnd = 0;
        final SamReader samReader = openBAMReader(bamFile, ValidationStringency.SILENT, false);
        
        CigarOps cigarOper = new CigarOps();
        readFaIdxFile(fastaIndexFile);
        // access user specified region
        if (chrom != null){  
            
            SAMFileHeader samFileHeader = samReader.getFileHeader();
            
            SAMSequenceDictionary sequenceDictionary = samFileHeader.getSequenceDictionary();
            SAMSequenceRecord refSequenceRecord = sequenceDictionary.getSequence(chrom);
            int refSequenceLength = refSequenceRecord.getSequenceLength();           
            int nWindows = refSequenceLength / readDepthContainerBuffer;
            
                        
            
            if (regionStart != 0 && regionEnd != 0){
                refSequenceLength = regionEnd - regionStart + 1;
                if (refSequenceLength <= readDepthContainerBuffer){
                    nWindows = 1;
                    readDepthContainerBuffer = refSequenceLength;
                }else{
                    nWindows = refSequenceLength/readDepthContainerBuffer;
                }
                windowStart = regionStart;
                
                
            }
            
            int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
            
            for (int i = 0; i < nWindows; i++){  
                windowEnd = windowStart + readDepthContainerBuffer;                 
                SAMRecordIterator iterator = samReader.query(chrom, windowStart, windowEnd, false);               
                
                analysisAlignment(iterator, windowStart, cigarOper);                
                processRemainingSignals();
                int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                System.out.println("processed region: [" + windowStart + ", " + windowEnd + "] " + "#superitems: " + curBinSuperitemCount);
                windowStart = windowEnd;
                
//                readDepthPreStepBuffer = copyFromReadDepthBuffer();
                readDepthPreStepBuffer = readDepthContainer;
                
                readDepthContainer = new int[readDepthContainerBuffer];
               
                writeAllSuperItems();
                    
            }
            // process remaining alignment in BAM
            SAMRecordIterator iterator = samReader.query(chrom, windowStart, refSequenceLength, false);
            analysisAlignment(iterator, windowStart, cigarOper);
            processRemainingSignals();
            int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
            
            System.out.println("processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#superitems: " + curBinSuperitemCount);
            
            
            writeAllSuperItems();
           
           
        } 
        // read whole genome
        else{            
            int length = chromLength.length;
            SAMRecordIterator iterator;
            for (int i = 0;i < length; i ++){
                int refSequenceLength = chromLength[i];
                String curChrom = chromInfoMap[i];
                System.out.println("Start processing chrom: " + curChrom + ", chrom length: " + refSequenceLength);
                
                int nWindows = refSequenceLength / readDepthContainerBuffer;

                int[] readDepthPreStepBuffer = new int[readDepthContainerBuffer];
                for (int k = 0; k < nWindows; k++){
                    windowEnd = windowStart + readDepthContainerBuffer;                 
                    iterator = samReader.query(curChrom, windowStart, windowEnd, false);
                    
                    analysisAlignment(iterator, windowStart, cigarOper);                
                    processRemainingSignals();
                    int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                    System.out.println("processed region: [" + windowStart + ", " + windowEnd + "] " + "#superitems: " + curBinSuperitemCount);
                    windowStart = windowEnd;

                    readDepthPreStepBuffer = copyFromReadDepthBuffer();
                    readDepthContainer = new int[readDepthContainerBuffer];
                    
                    writeAllSuperItems();                                                                
                    
                }
                iterator = samReader.query(curChrom, windowStart, refSequenceLength, false);
                analysisAlignment(iterator, windowStart, cigarOper);
                processRemainingSignals();
                int curBinSuperitemCount = assignReadDepthAndCountSuperItem(windowStart, readDepthContainerBuffer, readDepthPreStepBuffer); 
                System.out.println("processed region: [" + windowStart + ", " + refSequenceLength + "] " + "#superitems: " + curBinSuperitemCount);
                
                writeAllSuperItems();
                windowStart = 0;
                windowEnd = 0;
            
                rpTracker.clear();
            }
        }                                               
    }        
    /**
     * Analysis each BAM record through different channels
     * @param iterator
     * @param windowStart
     * @param cigarOper 
     */
    private void analysisAlignment(SAMRecordIterator iterator, int windowStart, CigarOps cigarOper){
//        CigarOps corasenCigar = new CigarOps();
        while(iterator.hasNext()){
            SAMRecord record = iterator.next();
            int mapq = record.getMappingQuality();            
            int recordChrIdx = record.getReferenceIndex();
            String recordChrName = record.getReferenceName();
            if (!idxToChromMap.containsKey(recordChrIdx)){
                idxToChromMap.put(recordChrIdx, recordChrName);
            }
            // Discard reads of low quality and PCR duplicated reads
            if (mapq <= minMapQ || record.getDuplicateReadFlag()){
                continue;
            }            
            List<CigarElement> cigarElements = record.getCigar().getCigarElements();
            
            if (badReads(cigarElements)){
                continue;
            }
            // count the number of normal read per base
            int isGoodAlign = exactAlignment(cigarElements);
            if (isGoodAlign != -1){
                updateReadDepthArray(record.getAlignmentStart(), isGoodAlign, windowStart);
            }                                      
           
            SEClippedParser(record, cigarElements, cigarOper);   
            RPUnmappedParser(record);
            RPisizeParser(record, cigarElements);
            
            if (!rpTracker.containsKey(record.getReadName())){
                List<SAMRecord> records = new ArrayList<>();
                records.add(record);
                rpTracker.put(record.getReadName(), records);
            }else{
                List<SAMRecord> records = rpTracker.get(record.getReadName());
                records.add(record);
                numOfRPs += 1;
                RPoriParser(records, cigarElements);
                rpTracker.remove(record.getReadName());
            }                
            
            
        } 
        iterator.close();
    }
    private int[] copyFromReadDepthBuffer(){
        int[] newBuffer = new int[readDepthContainerBuffer];
//        int startPosToCopy = readDepthContainerBuffer - readLen;
        for (int i = 0; i < readDepthContainerBuffer; i ++){
            int val = readDepthContainer[i];
//            newBuffer[i - startPosToCopy] = val;
            newBuffer[i] = val;
        }
        return newBuffer;
    }
    /**
     * Discard reads of clipped length longer than 70% of the read length 
     * @param cigarElements
     * @return 
     */
    private boolean badReads(List<CigarElement> cigarElements){
        int clippedLength = 0;
        boolean isBad = false;
        for (CigarElement element : cigarElements){
            String operation = element.getOperator().toString();
            int optLength = element.getLength();
            if (operation.equals("S") || operation.equals("H")){
                clippedLength += optLength;
            }
        }
        if (clippedLength > 0.7 * readLen){
            isBad = true;
        }
        return isBad;
    }
    /**
     * Get exact matched read
     * @param cigarElements
     * @return 
     */
    private int exactAlignment(List<CigarElement> cigarElements){
        if (cigarElements.size() == 1){
            String cigarOperation = cigarElements.get(0).getOperator().toString();
            int opLength = cigarElements.get(0).getLength();
            return cigarOperation.equals("M") ? opLength : -1;
        }
        else return -1;
    }
    
    private void updateReadDepthArray(int pos, int length, int windowStart){
        
        for (int i = 0; i < length ; i ++){
            if ( (pos + i - windowStart) >= readDepthContainerBuffer){
                continue;
            }
            else if (pos + i < windowStart) {
                continue;
            }
            else{
                readDepthContainer[pos + i - windowStart] += 1;            
            }
            
        }
        
    }
//    private int isOverlapRP(List<SAMRecord> records){
//        int overlap = -1;
//        SAMRecord leftMostAlign = records.get(0);
//        int leftMostAlignStart = leftMostAlign.getAlignmentStart();
//        SAMRecord mateAlign = records.get(1);
//        int mateAlignStart = mateAlign.getAlignmentStart();
//        int tmp = mateAlignStart - leftMostAlignStart;
//        if (tmp < readLen){
//            overlap = tmp;
//        }
//        return overlap;
//    }   
    /**
     * Process soft clipped reads
     * @param record
     * @param cigarElements
     * @param cigarOper 
     */
    
    private void SEClippedParser(SAMRecord record, List<CigarElement> cigarElements, CigarOps cigarOper) {
        // For a mapped read and read of relatively high mapQ
        if (!record.getReadUnmappedFlag()){

            String firstOperation = cigarElements.get(0).getOperator().toString();
                        
            
            cigarOper.calQueryPosFromCigar(cigarElements, 1, record.getReadNegativeStrandFlag(),readLen);
            int qsPos = cigarOper.getqsPos();
            
            String cigarStr = cigarOper.getCigarStr();
            int mutCoord = record.getAlignmentStart();
            
            if (!cigarStr.equals("M") && !cigarStr.isEmpty()){
                if (firstOperation.equals("M")){                    
                    mutCoord += qsPos;                
                }
                if (cigarOper.isCoIDread() && firstOperation.equals("S")){
                    mutCoord += qsPos;
                }
                
                String ori = record.getReadNegativeStrandFlag() ? "-" : "+";                               
                MutSignal mutSignal = new MutSignal(record, cigarStr, mutCoord, ori);
                mutSignal.setIsizeNormal(isizeUpper, isizeLower);

                breakChannel.addSignals(mutSignal, fragMean, readLen);
            }
        }
        
    }
    
     /**
     * One end unmapped read
     * @param record
     * @return 
     */
    private void RPUnmappedParser(SAMRecord record){
        
        // read unmapped
        if (record.getReadUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
           
            oemChannel.addSignals(mutSignal, fragMean, readLen);
            numOfARPs += 1;

        }else if (record.getMateUnmappedFlag()){
            int mutCoord = record.getMateAlignmentStart();
            String ori = record.getMateNegativeStrandFlag() ? "-": "+";
            
            MutSignal mutSignal = new MutSignal(record, "ARP_OEM", mutCoord, ori);
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
                        
            oemChannel.addSignals(mutSignal, fragMean, readLen);
            numOfARPs += 1;
        }

    }
    /**
     * Process PE of abnormal insert size
     * @param record
     * @param cigarElements 
     */
    private void RPisizeParser(SAMRecord record, List<CigarElement> cigarElements){
        // only process read-pair mapped on the same chrom.
       
        CigarElement leftMostCigarElement = cigarElements.get(0);
        String leftMostCigarOperator = leftMostCigarElement.getOperator().toString();

        int mutCoord = record.getAlignmentStart();
        if (leftMostCigarOperator.equals("M") && !record.getReadNegativeStrandFlag()){
            mutCoord += leftMostCigarElement.getLength();
        }

        int insertSize = record.getInferredInsertSize();

        String ori = record.getReadNegativeStrandFlag() ? "-" : "+";
        if (Math.abs(insertSize) >= isizeUpper){    
//                System.out.println(record.getReadName() + " isize: " +Math.abs(insertSize));

            MutSignal mutSignal = new MutSignal(record, "ARP_LARGE_INSERT", mutCoord, ori);               
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                        record.getInferredInsertSize(), "ARP_LARGE_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);
            isizeLargeChannel.addSignals(mutSignal, fragMean, readLen);            
            numOfARPs += 1;
        }
        else if (Math.abs(insertSize) <= isizeLower && insertSize != 0){

            MutSignal mutSignal = new MutSignal(record, "ARP_SMALL_INSERT", mutCoord, ori);       
//            MutSignal mutSignal = new MutSignal(record.getReadName(), record.getReferenceIndex(), record.getReferenceName(), 
//                record.getInferredInsertSize(), "ARP_SMALL_INSERT", mutCoord, ori, record.getAlignmentStart(), record.getMateAlignmentStart());
            mutSignal.setIsizeNormal(isizeUpper, isizeLower);

            isizeSmallChannel.addSignals(mutSignal, fragMean, readLen); 
            numOfARPs += 1;
            
        }
                
    }
    /**
     * Process PE of abnormal alignment orientation
     * @param records
     * @param cigarElements 
     */
    private void RPoriParser(List<SAMRecord> records, List<CigarElement> cigarElements){
        
        SAMRecord leftMostRecord = records.get(0);
        SAMRecord rightMostRecord = records.get(records.size() - 1);
        // For read-pair, it should be proper paired. Its read and mate are all mapped.
        if (leftMostRecord.getReadPairedFlag() && !leftMostRecord.getReadUnmappedFlag() && !leftMostRecord.getMateUnmappedFlag()){
            int mutCoord = leftMostRecord.getAlignmentStart();
            if (leftMostRecord.getReadNegativeStrandFlag()== leftMostRecord.getMateNegativeStrandFlag()){
                String mutType;
                String ori = leftMostRecord.getReadNegativeStrandFlag() ? "-":"+";
                if (leftMostRecord.getReadNegativeStrandFlag()){
                    mutType = "ARP_RR";
                }else{
                    CigarElement leftMostCigarElement = cigarElements.get(0);
                    String leftMostCigarOperation = leftMostCigarElement.getOperator().toString();

                    if (leftMostCigarOperation.equals("M")){
                        mutCoord += leftMostCigarElement.getLength();
                    }
                    mutType = "ARP_FF";
                }
                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, ori);
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                MutSignal mateMutSignal = new MutSignal(leftMostRecord, mutType, rightMostRecord.getMateAlignmentStart(), ori);
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);   

                oriChannel.addSignals(readMutSignal, fragMean, readLen);
                oriChannel.addSignals(mateMutSignal, fragMean, readLen);
                numOfARPs += 1;
            }
            else if (leftMostRecord.getReadNegativeStrandFlag() && !leftMostRecord.getMateNegativeStrandFlag()){
                String mutType = "ARP_RF";

                MutSignal readMutSignal = new MutSignal(leftMostRecord, mutType, mutCoord, "-");
                readMutSignal.setIsizeNormal(isizeUpper, isizeLower);
                MutSignal mateMutSignal = new MutSignal(leftMostRecord, mutType, rightMostRecord.getMateAlignmentStart(), "+");
                mateMutSignal.setIsizeNormal(isizeUpper, isizeLower);

                oriChannel.addSignals(readMutSignal, fragMean, readLen);
                oriChannel.addSignals(mateMutSignal, fragMean, readLen);

                numOfARPs += 1;
            }
        }       
    }

    
    public Map<Integer, String> getIdxToChromMap(){
        return idxToChromMap;
    } 
    
    private SamReader openBAMReader(String bamFile, ValidationStringency stringency, boolean includeFileSource) throws IOException{
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(stringency).enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
        if(includeFileSource){
            samReaderFactory.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS);            
        }
//        samReaderFactory.samRecordFactory(DefaultSAMRecordFactory.getInstance());
        final SamReader samReader = samReaderFactory.open(new File(bamFile));
        return samReader;
    }
    
    private void readFaIdxFile(String fastaIndexFile) throws IOException{
        chromInfoMap = new String[24];
        chromLength = new int[24];
        FileInputStream fin = new FileInputStream(new File(fastaIndexFile));
        BufferedReader myInput = new BufferedReader(new InputStreamReader(fin));
        String thisLine;
        int lineIdx = 0;
        while ((thisLine = myInput.readLine()) != null){
            String[] tokens = thisLine.split("\t");
            String refSeq = tokens[0];
            // escape "M", "MT", "chrM"
            if (refSeq.contains("M") || refSeq.contains("_")){
                continue;
            }
                                
            chromLength[lineIdx] = Integer.parseInt(tokens[1]);
            chromInfoMap[lineIdx] = refSeq;
            
            lineIdx += 1;
            // Read until Y chromosome, ignore other contigs.
            if (refSeq.equals("Y") || refSeq.equals("chrY")){
                break;
            }
        }
        
    }
       
    private void processRemainingSignals() {
        breakChannel.processFinalSignals(fragMean, readLen);
        isizeLargeChannel.processFinalSignals(fragMean, readLen);
        isizeSmallChannel.processFinalSignals(fragMean, readLen);
        oemChannel.processFinalSignals(fragMean, readLen);
        oriChannel.processFinalSignals(fragMean, readLen);
//        interChromChannel.processFinalSignals(fragMean, readLen);
                             
    }
    
    private void writeAllSuperItems() throws IOException{
        if (writer != null){
            breakChannel.writeSuperItemsInChannel(writer);
            isizeLargeChannel.writeSuperItemsInChannel(writer);
            isizeSmallChannel.writeSuperItemsInChannel(writer);
            oemChannel.writeSuperItemsInChannel(writer);
            oriChannel.writeSuperItemsInChannel(writer);
//            interChromChannel.writeSuperItemsInChannel(writer);
        }
        
    }
    /**
     * calculate normal read aligned at a specific position and the number of superitems that generated within this window.
     * @param windowStart
     * @param windowSize
     * @param preReadDepthBuffer
     * @return 
     */
    private int assignReadDepthAndCountSuperItem(int windowStart, int windowSize, int[] preReadDepthBuffer){
        breakChannel.setSuperitemWeightRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeLargeChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        isizeSmallChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oriChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        oemChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
//        interChromChannel.setARPSuperItemRatio(readDepthContainer, windowStart, windowSize, preReadDepthBuffer);
        
        int curWindowSuperItem = 0;
        curWindowSuperItem += breakChannel.getSuperitemCount();
        curWindowSuperItem += isizeLargeChannel.getSuperitemCount();
        curWindowSuperItem += isizeSmallChannel.getSuperitemCount();
        curWindowSuperItem += oemChannel.getSuperitemCount();
        curWindowSuperItem += oriChannel.getSuperitemCount();
//        curWindowSuperItem += interChromChannel.getSuperitemCount();
        
        
        return curWindowSuperItem;
    }
    private void createSuperItemWriter(String superitemOutPath) throws IOException{
        writer = new BufferedWriter(new FileWriter(superitemOutPath));
        writer.write("type\tchromIdx\tnread\tpos\tsaPos\tori\tweight\tratio\tsumMapQ\tplusRead\tminusRead\tsplitRead\tsplitMapQ\titxRead\tregion\tmateRegion\tqnames\tmConsensus\tcConsensus\n");

    }
    
    public int getWGARPNum(){
        return numOfARPs;
    }    
    public void printSuperItemGeneratorStats(){
        StringBuilder sb = new StringBuilder();
        sb.append("\n==============  SuperItem Generation =============\n");
        sb.append("Time: " + (endTime - startTime) + "ms");  
        sb.append("\nTotal superitems: " + superitemCount);   
        sb.append("\nTotal number of read-pairs: " + numOfRPs);
        sb.append("\nTotal number of doiscordant read-pairs: " + numOfARPs);
        sb.append("\n======================================================");
        
        System.out.println(sb.toString());
    }
    
}
