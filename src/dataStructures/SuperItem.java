/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package dataStructures;

import java.util.*;
import java.util.Map.Entry;
import Channels.MutSignal;
import htsjdk.samtools.QueryInterval;
import java.text.DecimalFormat;
import htsjdk.samtools.util.*;
/**
 *
 * @author jiadonglin
 */
public class SuperItem implements Comparable<SuperItem>{
    

    protected int superItemPos;
    protected int weight;
    protected String type;
    protected String ori;
    protected double entropy;
    protected QueryInterval superitemInterval;
    protected QueryInterval superitemMateInterval;

    protected int chromIdx;
    protected String chromName;
    protected int nread = 0;
    protected double weightRatio = 0;
    protected String[] qnames;
    protected int splitAlignedPos = -1;
    protected int numPlusRead;
    protected int numMinusRead;
    protected int sumMapQ;
    protected int splitReadMapQ;
    protected int splitAlignedRead;
    protected String matchedConsensus = "";
    protected String clippedConsensus = "";
//    protected String readConsensus; 
    
    
    public SuperItem(){
        
    }
    // construct SuperItem from file records
    public SuperItem(String[] tokens){
        
        type = tokens[0];
        chromIdx = Integer.parseInt(tokens[1]);
        
        nread = Integer.parseInt(tokens[2]);
        superItemPos = Integer.parseInt(tokens[3]);
        splitAlignedPos = Integer.parseInt(tokens[4]);
        ori = tokens[5];
        weight = Integer.parseInt(tokens[6]);           
        weightRatio = Double.parseDouble(tokens[7]);
        
        sumMapQ = Integer.parseInt(tokens[8]);
        numPlusRead = Integer.parseInt(tokens[9]);
        numMinusRead = Integer.parseInt(tokens[10]);
        splitAlignedRead = Integer.parseInt(tokens[11]);
        splitReadMapQ = Integer.parseInt(tokens[12]);
        
        superitemInterval = decodeIntervalFromString(tokens[13]);
                
        superitemMateInterval = decodeIntervalFromString(tokens[14]);
                
        String qNameColumn = tokens[15];
        if (!qNameColumn.equals("*")){
            qnames = qNameColumn.split(",");

        }
        matchedConsensus = tokens[16];
        clippedConsensus = tokens[17];
              
    }
    
    // Create SuperItem in memory
    public SuperItem(List<MutSignal> mutSignals){
        createSuperItem(mutSignals);        
    }
    
    @Override
    public int compareTo(SuperItem otherSuperItem){
        return this.superItemPos - otherSuperItem.getPos();
    }
    /**
     * Consider the equality of superitem. If all attribute considered, too much identical superitems.
     * @param obj
     * @return 
     */
    public boolean equals(Object obj){
        if (obj instanceof SuperItem){
            SuperItem si = (SuperItem) obj;
            return (si.type.equals(this.type));
        }else{
            return false;
        }
    }
    
    
    public boolean isEqual(SuperItem other){
        return ((type.equals(other.type)) && (superItemPos == other.superItemPos));
    }
    public boolean isSmallIndel(){
        if (!type.contains("ARP") && (type.contains("I")||type.contains("D"))){
            return true;
        }
        else return false;
    }
    
    // debug useage
    public String toConciseString(){
        StringBuilder sb = new StringBuilder();
        sb.append("id:" + chromIdx);
        sb.append(" t:" + type);
        sb.append(" p:" + superItemPos);
        sb.append(" w:" + weight);
        sb.append(" r:" + weightRatio);
        return sb.toString();
    }    
    
//    public String toItemType(){
//        
//    }
    
    private String qNamesToString(){
        String str = "*";
        if (qnames.length != 0){
            StringBuilder sb = new StringBuilder();
            for (String qname : qnames){
                sb.append(qname);
                sb.append(",");
            }
            str = sb.substring(0, sb.length() - 1);
        }        
        return str;
    }
    
    @Override
    public String toString(){
        StringBuilder sb = new StringBuilder();
        DecimalFormat twoDForm = new DecimalFormat("#.###");
        double roundedRatio = Double.valueOf(twoDForm.format(weightRatio));
//        double roundedEntropy = Double.valueOf(twoDForm.format(entropy));
        sb.append(type);
        sb.append("\t");
        sb.append(chromIdx);
        sb.append("\t");
        sb.append(nread);
        sb.append("\t");
        sb.append(superItemPos);
        sb.append("\t");        
        sb.append(splitAlignedPos);
        sb.append("\t");              
        sb.append(ori);
        sb.append("\t");
        sb.append(weight);
        sb.append("\t");
        sb.append(roundedRatio);
        sb.append("\t");
        sb.append(sumMapQ);
        sb.append("\t");
        sb.append(numPlusRead);
        sb.append("\t");
        sb.append(numMinusRead);
        sb.append("\t");
        sb.append(splitAlignedRead);
        sb.append("\t");
        sb.append(splitReadMapQ);
        sb.append("\t");
        sb.append(superitemInterval.toString());
        sb.append("\t");
        sb.append(superitemMateInterval.toString());
        sb.append("\t");
        sb.append(qNamesToString());
        sb.append("\t");
        if (matchedConsensus.equals("")){
            matchedConsensus = "-";
        }
        sb.append(matchedConsensus);  
        sb.append("\t");
        if (clippedConsensus.equals("")){
            clippedConsensus = "-";
        }
        sb.append(clippedConsensus);
        return sb.toString();
    }
    @Override
    public int hashCode(){
        return type.hashCode();
    }
    public String[] getQNames(){
        return qnames;
    }
    public boolean isARPsuperitem(){
        return type.contains("ARP");
    }
    public QueryInterval getSuperitemRegion(){
        return superitemInterval;
    }
    public QueryInterval getSuperitemMateRegion(){
        return superitemMateInterval;
    }
    public double getWeightRatio(){
        return this.weightRatio;
    }
    public String getType(){
        return this.type;
    }
    public int getPos(){
        return this.superItemPos;
    }
    public int getWeight(){
        return this.weight;
    }
    public double getRatio(){
        return this.weightRatio;
    }
    public int getChromIdx(){
        return chromIdx;
    }
    public void setChromName(String chrom) {
        chromName = chrom;
    }
    public String getChromName(){
        return chromName;
    }
    public String getOri(){
        return ori;
    }
    public int getSplitAlignPos(){
        return splitAlignedPos;
    }
    public double getEntropy(){
        return entropy;
    }
    public int getSumMapQ(){
        return sumMapQ;
    }
    public int getNumMinusRead(){
        return numMinusRead;
    }
    
    public int getNumPlusRead(){
        return numPlusRead;
    }
    
    public String getReadOri(){
        return numPlusRead + "+" + numMinusRead + "-";
    }
    public int getSplitReadCount(){
        return splitAlignedRead;
    }
    
    public int getSplitReadMapQ(){
        return splitReadMapQ;
    }
    public double getMinusPlusRatio(){
        double ratio = Double.MAX_VALUE;
        if (numPlusRead != 0 && numMinusRead >= numPlusRead){
            ratio = (double) numMinusRead / numPlusRead;
        }                
        if (numMinusRead != 0 && numMinusRead < numPlusRead){
            ratio = (double) numPlusRead / numMinusRead;
        }        
        return ratio;
    }
    
    public String getMachedConsensus(){
        return matchedConsensus;
    }
    public String getClippedConsensus(){
        return clippedConsensus;
    }
    
    /**
     * Use the dominant signals
     * @param mutList sorted mutational signals by mutPos
     */
    private void createSuperItem(List<MutSignal> mutList){
        
        Map<MutSignal, Integer> typeCountMap = new HashMap<>();
        
        List<Integer> superitemRecordPos = new ArrayList<>();
        List<Integer> superitemMateRecordPos = new ArrayList<>();
        
        Set<String> queryName = new HashSet<>();
        qnames = new String[mutList.size()];
        
        int longestD = -1;
        int signalIdx = 0;
        // sort the mutational signal to get the SuperItem genome intervals
        for (MutSignal signal : mutList){   
            sumMapQ += signal.getMapQ();            
            if (signal.getMutSignalOri().equals("-")){
                numMinusRead += 1;
            }else{
                numPlusRead += 1;
            }
            if(signal.isSplitAlign()){
                splitAlignedPos = signal.getSplitAlignPos();
                splitAlignedRead += 1;
                splitReadMapQ += signal.getMapQ();
            }
            else if (signal.getMutSignalType().equals("MDM") && longestD == -1){
                longestD = signal.getLongestD();
            }
            String signalQname = signal.getqName();
            superitemRecordPos.add(signal.getRecordPos());
            superitemMateRecordPos.add(signal.getMateRecordPos());                        
            
            if (!queryName.contains(signalQname)){                
                queryName.add(signalQname);
                qnames[signalIdx] = signalQname;
                signalIdx += 1;
            }
            
            if (!typeCountMap.containsKey(signal)){
                typeCountMap.put(signal, 1);
            }else{
                int count = typeCountMap.get(signal);
                count += 1;
                typeCountMap.put(signal, count);
            }
        }
        
        Collections.sort(superitemRecordPos);
        Collections.sort(superitemMateRecordPos);
        
        List<Integer> superItemOutlierRemoved = removePosOutlier(superitemRecordPos);
        List<Integer> superItemMateOutlierRemoved = removePosOutlier(superitemMateRecordPos);
        
        
        Collections.sort(mutList, new Comparator<MutSignal>(){
            @Override
            public int compare(MutSignal o1, MutSignal o2){
                return o1.getRecordPos() - o2.getRecordPos();
            }
        });
        
        MutSignal majoritySignal = new MutSignal();
        int maxCount = 0;
        List<Integer> countList = new ArrayList<>();
        for (Entry<MutSignal, Integer> entry : typeCountMap.entrySet()){
            countList.add(entry.getValue());
            if (entry.getValue() > maxCount){
                maxCount = entry.getValue();
                majoritySignal = entry.getKey();                
            }
        }
                 
        type = majoritySignal.getMutSignalType();        
        ori = majoritySignal.getMutSignalOri();
        chromIdx = majoritySignal.getSignalChromIdx();
        chromName = majoritySignal.getSignalRef();
        superItemPos = majoritySignal.getMutPos();
        weight = mutList.size();   
        
//        if (superItemPos == 258783){
//            System.out.println("ssss");
//        }
        
        if (longestD != -1){
            splitAlignedPos = superItemPos + longestD;
        }        
                
        setSuperItemEntropy(countList); 
        if (!type.equals("SMS")){
            buildConsensusString(mutList);
            if (("SM").equals(type)){
                StringBuilder sb = new StringBuilder(clippedConsensus);
                clippedConsensus = sb.reverse().toString();
            }
        }        
//        buildReadConsensusString(mutList);
        
        // If this is a discordant read-pair based superitem, adjust the position and assign query names for further process.
        if (type.contains("ARP")){
//            qnames = qNameByteList;
            superitemInterval = new QueryInterval(chromIdx, superItemOutlierRemoved.get(0), superItemOutlierRemoved.get(superItemOutlierRemoved.size() - 1));
            superitemMateInterval = new QueryInterval(chromIdx, superItemMateOutlierRemoved.get(0), superItemMateOutlierRemoved.get(superItemMateOutlierRemoved.size() - 1));
                        
            if (ori.equals("+")){
                superItemPos = mutList.get(mutList.size() - 1).getMutPos() - 1;
            }else{
                superItemPos = mutList.get(0).getMutPos() + 1;
            }

        }else{
//            qnames = qNameByteList;
            // make sure there is a potential second majority superitem
            if (type.equals("SMS") && typeCountMap.size() > 1){
                refineSuperItemType(typeCountMap, mutList);
            }    
            superitemInterval = new QueryInterval(chromIdx, superItemPos, superItemPos);
            superitemMateInterval = new QueryInterval(chromIdx, superItemMateOutlierRemoved.get(0), superItemMateOutlierRemoved.get(superItemMateOutlierRemoved.size() - 1));                        
        }        
    }   
    
    private List<Integer> removePosOutlier(List<Integer> coords){
        List<Integer> dedupedCoords = new ArrayList<>();
        dedupedCoords.add(coords.get(0));
        for (Integer val : coords){
            int lastVal = dedupedCoords.get(dedupedCoords.size() - 1);
            if (lastVal != val){
                dedupedCoords.add(val);
            }
        }
        int len = dedupedCoords.size();
        int mid = len / 2;
        List<Integer> beforeMid = new ArrayList<>();
        List<Integer> afterMid = new ArrayList<>();
        for (int i = 0 ; i < len ; i++){
            if (i <= mid){
                beforeMid.add(dedupedCoords.get(i));
            }else{
                afterMid.add(dedupedCoords.get(i));
            }
        }
        
        int firstQuarter = -1;
        int thirdQuarter = -1;
        if (beforeMid.size() > 1 && afterMid.size() > 1){
            firstQuarter = getMedian(beforeMid);
            thirdQuarter = getMedian(afterMid);
        }
        List<Integer> outlierRemoved = new ArrayList<>();
        
        if (firstQuarter != -1 && thirdQuarter != -1){
            for (Integer val : coords){
                if (val >= firstQuarter && val <= thirdQuarter){
                    outlierRemoved.add(val);
                }
            }
        }else{
            outlierRemoved = coords;
        }
        
        return outlierRemoved;
    }
    
    private int getMedian(List<Integer> numList){
        int len = numList.size();
        int median = -1;
        if (len % 2 == 0){
            median = (numList.get(len >> 1) + numList.get((len >> 1) - 1)) / 2;
        }else{
            median = numList.get(len >> 1);
        }
        return median;
    }
    
    public void setSuperitemReadDepth(int rd){
        nread = rd;
        weightRatio = (double) weight / (weight + nread);
    }
    public void setARPsuperitemRatio(double avgCov){
        nread = (int) avgCov;
        weightRatio = weight / (weight + avgCov);
    }
    /**
     * Refine superitem of type SMS
     * @param signalCountMap 
     */
    private void refineSuperItemType(Map<MutSignal, Integer> signalCountMap, List<MutSignal> mutList){
        MutSignal secondMajoritySignal = new MutSignal();
        int maxCount = 0;
        for(Entry<MutSignal, Integer> entry : signalCountMap.entrySet()){
           String signalType = entry.getKey().getMutSignalType();
           if (!signalType.equals("SMS")){
               if (entry.getValue() > maxCount){
                   maxCount = entry.getValue();
                   secondMajoritySignal = entry.getKey();
               }
           }
        }
        
        if (secondMajoritySignal != null){
           type = secondMajoritySignal.getMutSignalType();
           buildConsensusString(mutList);
        }
       
    }
    private QueryInterval decodeIntervalFromString(String str){
        int colonIndex = str.lastIndexOf(':');
        int start = 0;
        int end = 0;
        if (colonIndex == -1){
            start = 1;
            end = Integer.MAX_VALUE;
        }else{
            int dashIndex = str.indexOf('-', colonIndex);
            if (dashIndex == - 1){
                System.err.println("Incorrect interval format !");
            }else{
                start = Integer.parseInt(str.substring(colonIndex + 1, dashIndex).replaceAll(",", ""));
                end = Integer.parseInt(str.substring(dashIndex + 1).replaceAll(",", ""));               
            }
        } 
        return new QueryInterval(chromIdx, start, end);
    }
    
    public void setSuperItemEntropy(List<Integer> countList){      
        double ent = 0;
        int sum = 0;
        for (Integer ele : countList){
            sum += ele;
        }        
        for(int i = 0; i < countList.size(); i++){
            double prob = (double) countList.get(i) / (double) sum;
            double val = prob * (Math.log(prob)/Math.log(2));
            ent-= val;
        }        
        entropy = ent;
    }
    
    private void buildReadConsensusString(List<MutSignal> mutList){
        int mutListStartPos = mutList.get(0).getRecordPos();        
        Map<Integer, Integer[]> readCount = new HashMap<>();

        
        for (MutSignal signal : mutList){
            if (!signal.getMutSignalType().equals(type)){
                continue;
            }            
            String readSeq = signal.getReadSequence();
            int readStartPos = signal.getRecordPos() - mutListStartPos;
            readBaseCount(readSeq, 0, readSeq.length(), readStartPos, readCount);            
        }
//        readConsensus = getReadConsensus(readCount);
    }
    
    private void buildConsensusString(List<MutSignal> mutList){    
        int mutListStartPos = mutList.get(0).getRecordPos();        
        Map<Integer, Integer[]> matchedNucleotide = new HashMap<>();
        Map<Integer, Integer[]> clippedNucleotide = new HashMap<>();
        
        for (MutSignal signal : mutList){
            if (!signal.getMutSignalType().equals(type)){
                continue;
            }
            
//            System.out.println(signal.getqName());
            
            int[] clippedStatus = signal.getClippedStatus();           
            String readSeq = signal.getReadSequence();
            int readSeqLength = readSeq.length();
            
            int readStartPos = signal.getRecordPos() - mutListStartPos;
            int readStringEndIdx = readSeqLength;
            
            if (clippedStatus[0] == 0){
                readStringEndIdx -= clippedStatus[1];                    
                readBaseCount(readSeq, 0, readStringEndIdx, readStartPos, matchedNucleotide);
                readBaseCount(readSeq, readStringEndIdx, readSeqLength, 0, clippedNucleotide);                
            }else{
                readBaseCount(readSeq, clippedStatus[0] + 1, readStringEndIdx, 0, matchedNucleotide);
                clipFirstBaseCount(readSeq, clippedStatus[0], clippedNucleotide);
            }
        }
        matchedConsensus = getConsensusString(matchedNucleotide, false);        
        clippedConsensus = getConsensusString(clippedNucleotide, true);        
    } 
    
    private void clipFirstBaseCount(String readString, int readStringEndAt, Map<Integer, Integer[]> readPosNuc){        
        for (int i = readStringEndAt; i >=0 ; i--){
            int keyInMap = readStringEndAt - i;
            char nucleotide = readString.charAt(i);
            Integer[] nucleotideCount = readPosNuc.get(keyInMap);
            if (nucleotideCount == null){
                Integer[] counts = new Integer[4]; 
                Arrays.fill(counts, 0);
                if (nucleotide == 'a' || nucleotide == 'A'){
                    counts[0] += 1;
                    readPosNuc.put(keyInMap, counts);
                    continue;
                }
                if (nucleotide == 't' || nucleotide == 'T'){
                    counts[1] += 1;
                    readPosNuc.put(keyInMap, counts);
                    continue;
                }
                if (nucleotide == 'c' || nucleotide == 'C'){
                    counts[2] += 1;                   
                    readPosNuc.put(keyInMap, counts);
                    continue;
                }
                if (nucleotide == 'g' || nucleotide == 'G'){
                    counts[3] += 1;
                    readPosNuc.put(keyInMap, counts);
                    
                }else{
                    readPosNuc.put(keyInMap, counts);                    
                }
            }else{
                if (nucleotide == 'a' || nucleotide == 'A'){
                    nucleotideCount[0] += 1;
                    continue;
                }
                if (nucleotide == 't' || nucleotide == 'T'){
                    nucleotideCount[1] += 1;
                    continue;
                }
                if (nucleotide == 'c' || nucleotide == 'C'){
                    nucleotideCount[2] += 1;  
                    continue;
                }
                if (nucleotide == 'g' || nucleotide == 'G'){
                    nucleotideCount[3] += 1;                    
                }
            }
            
        }
    }
    
    private void readBaseCount(String readString, int readStringStartAt, int readStringEndAt, int startPos, Map<Integer, Integer[]> readPosNuc){
        // sort mutational signals according to read alignment position   
        int length = readStringEndAt - readStringStartAt;
        for (int i = 0; i < length; i++){
            int keyInMap = startPos + i;
            char nucleotide = readString.charAt(readStringStartAt + i);                        
            Integer[] nucleotideCount = readPosNuc.get(keyInMap);
            if (nucleotideCount == null){
                Integer[] counts = new Integer[4]; 
                Arrays.fill(counts, 0);
                if (nucleotide == 'a' || nucleotide == 'A'){
                    counts[0] += 1;
                    readPosNuc.put(keyInMap, counts);
                    continue;
                }
                if (nucleotide == 't' || nucleotide == 'T'){
                    counts[1] += 1;
                    readPosNuc.put(keyInMap, counts);
                    continue;
                }
                if (nucleotide == 'c' || nucleotide == 'C'){
                    counts[2] += 1;                   
                    readPosNuc.put(keyInMap, counts);
                    continue;
                }
                if (nucleotide == 'g' || nucleotide == 'G'){
                    counts[3] += 1;
                    readPosNuc.put(keyInMap, counts);
                    
                }else{
                    readPosNuc.put(keyInMap, counts);
                    
                }
                
            }else{
                if (nucleotide == 'a' || nucleotide == 'A'){
                    nucleotideCount[0] += 1;
                    continue;
                }
                if (nucleotide == 't' || nucleotide == 'T'){
                    nucleotideCount[1] += 1;
                    continue;
                }
                if (nucleotide == 'c' || nucleotide == 'C'){
                    nucleotideCount[2] += 1;  
                    continue;
                }
                if (nucleotide == 'g' || nucleotide == 'G'){
                    nucleotideCount[3] += 1;
                    
                }
            }
        }
    }
    
    private String getReadConsensus(Map<Integer, Integer[]> baseMap){
        StringBuilder sb = new StringBuilder();

        for (Entry<Integer, Integer[]> entry : baseMap.entrySet()){
            Integer[] baseCount = entry.getValue();            
                        
            int maxCount = 0;
            int maxCountBaseIndex = -1;
            for (int i = 0; i < 4; i++){
                int count = baseCount[i];
                if (count > maxCount){
                    maxCountBaseIndex = i;
                    maxCount = count;
                }
            }
            if (maxCountBaseIndex == 0){
                sb.append("A");
                continue;
            }
            if (maxCountBaseIndex == 1){
                sb.append("T");
                continue;
            }
            if (maxCountBaseIndex == 2){
                sb.append("C");                
            }
            else{
                sb.append("G");
            }
        }        
        return sb.toString();
    }
    
    private String getConsensusString(Map<Integer, Integer[]> baseMap, boolean clipped){
        StringBuilder sb = new StringBuilder();        
        for (Entry<Integer, Integer[]> entry : baseMap.entrySet()){
            Integer[] baseCount = entry.getValue();            
            int maxCount = 0;
            int maxCountBaseIndex = -1;
            for (int i = 0; i < 4; i++){
                int count = baseCount[i];
                if (count > maxCount){
                    maxCountBaseIndex = i;
                    maxCount = count;
                }
            }
            if (maxCountBaseIndex == 0){
                sb.append("A");
                continue;
            }
            if (maxCountBaseIndex == 1){
                sb.append("T");
                continue;
            }
            if (maxCountBaseIndex == 2){
                sb.append("C");                
            }
            else{
                sb.append("G");
            }
        }               
        return sb.toString();
    }       
}
