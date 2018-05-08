/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package contiguousfspm;

import dataStructures.SequenceDatabase;
import dataStructures.SuperItem;
import utils.*;
import htsjdk.samtools.QueryInterval;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.HashSet;
import utils.svOutInfo;



/**
 *
 * @author jiadonglin
 */
public class pseudoSequentialPattern implements Comparable<pseudoSequentialPattern>{
    
    List<pseudoSuperItem> superitems = new ArrayList<>();
    List<pseudoSuperItem> ARPSuperItems = new ArrayList<>();
    List<pseudoSuperItem> MSSuperItems = new ArrayList<>();
    List<QueryInterval> superitemInterval;
    List<QueryInterval> superitemMateInterval;
    List<Integer> weights;
    List<Integer> postions;    
    List<Double> ratios;
    List<String> oris;   
    Map<QueryInterval, List<Integer>> indexMap = new HashMap<>();
    
    int ChromId = -1;    
    
    int patternLength;
    int numOfBPs;
    int patternLeftMostPos; // position of the first SuperItem in the pattern
    int patternRightMostPos; // position of the last SuperItem in the pattern
    int patternLeftMostIntervalStart; // read position of the left most interval of the pattern
    int[] linkedPatterns = new int[2]; // for a self linked pattern, save its linked coord for boundary estimation
    int numOfLinkedEvidence = 0;
    int splitReadSup = 0;
    int splitReadMapQ = 0;
    
    int splitAlignLeftMatch = -1;
    int splitAlignRightMatch = -1;
    int splitAlignMate = -1;
    int[] splitAlignCoords;
    int[] crossSupInfo;
    
    // some configuration of self-linked patterns
    int[] selfLinkedSuperItemMapQ;
    int[] selfLinkedSuperItemWeight;
    

    @Override
    public int compareTo(pseudoSequentialPattern other){                                          
        return patternLeftMostIntervalStart - other.patternLeftMostIntervalStart;
    }

    public List<QueryInterval> getSortedIntervals(){
        List<QueryInterval> sortedInterval = new ArrayList<>();
        for (QueryInterval val : superitemInterval){
            sortedInterval.add(val);
        }
        return sortedInterval;
    }

    public pseudoSequentialPattern(List<pseudoSuperItem> itemset, SequenceDatabase database, String clipped){
        superitems = itemset;
        patternLength = itemset.size();
        for (int i = 0; i < patternLength ; i ++){
            SuperItem superitem = superitems.get(i).getSuperItem(database);    
            if (ChromId == -1){
                ChromId = superitem.getChromIdx();
                break;
            }
        }  
    }

    public pseudoSequentialPattern (List<pseudoSuperItem> itemset, SequenceDatabase database){

        superitems = itemset;
        patternLength = itemset.size();
        superitemInterval = new ArrayList<>();
        superitemMateInterval = new ArrayList<>();
        postions = new ArrayList<>();
        weights = new ArrayList<>();
        ratios = new ArrayList<>();
        oris = new ArrayList<>();
        
        patternLeftMostPos = itemset.get(0).getSuperItem(database).getPos();
        patternRightMostPos = itemset.get(itemset.size() - 1).getSuperItem(database).getPos();

        for (int i = 0; i < patternLength ; i ++){
            SuperItem superitem = superitems.get(i).getSuperItem(database);   

//                System.out.println(superitem.toConciseString());

            weights.add(superitem.getWeight());            
            postions.add(superitem.getPos());
            ratios.add(superitem.getWeightRatio());
            oris.add(superitem.getReadOri());
            
            if (ChromId == -1){
                ChromId = superitem.getChromIdx();
            }
            if (superitem.isARPsuperitem() && ! superitem.getType().contains("OEM")){                      
                ARPSuperItems.add(superitems.get(i));                                         
                superitemInterval.add(superitem.getSuperitemRegion());
                superitemMateInterval.add(superitem.getSuperitemMateRegion());
            }else{
                MSSuperItems.add(superitems.get(i));                    
                numOfBPs += 1;
            }
        }   
        Collections.sort(ARPSuperItems);
        Collections.sort(superitemInterval);
        if (!superitemInterval.isEmpty()){
            patternLeftMostIntervalStart = superitemInterval.get(0).start;
        }

    } 
    public boolean hasSplitAlign(SequenceDatabase database){
        boolean hasSplit = false;
        if (splitAlignCoords[0] != 0 && splitAlignCoords[1] != 0){
            hasSplit = true;
        }
        return hasSplit;
    }
    
    public int[] getSplitAlignCoords(){
        return splitAlignCoords;
    }
    
    public int[] getCrossSupInfo(){
        return crossSupInfo;
    }
    public List<Integer> getIndex(QueryInterval aInterval){
        return indexMap.get(aInterval);
    }
    public List<Integer> getWeights(){
        return weights;
    }
    
    public List<Double> getWeightRatio(){
        return ratios;
    }
    public List<Integer> getPos(){
        return postions;
    }
    public List<String> getOris(){
        return oris;
    }
    public int getLinkSupport(){
        return numOfLinkedEvidence;
    }

    public int getNumOfBPs(){
        return numOfBPs;
    }
    public List<pseudoSuperItem> getMSSuperItems(){
        return MSSuperItems;
    }
    public boolean hasArpSuperItem(){
        return !ARPSuperItems.isEmpty();
    }
    public boolean hasClippedSuperItem(){
        return !MSSuperItems.isEmpty();
    }
    public pseudoSuperItem getSuperItemAt(int idx){
        return superitems.get(idx);
    }
    
    public int getSplitSupCount(){
        return splitReadSup;
    }
    
    public String getChrName(int dbIdx, String[] dbIdxToChr){
        return dbIdxToChr[dbIdx];
    }
    
    public int[] getSelfLinkedItemMapQ(){
        return selfLinkedSuperItemMapQ;
    }
    
    public int[] getSelfLinkedItemWeight(){
        return selfLinkedSuperItemWeight;
    }
    public int getSplitReadMapQ(){
        return splitReadMapQ;
    }
    
    public String toString(SequenceDatabase database){
        StringBuilder sb = new StringBuilder();
        for (pseudoSuperItem item : superitems){
            SuperItem superitem = item.getSuperItem(database);
//            sb.append('(');
//            sb.append(superitem.toConciseString());
//            sb.append(')');
            sb.append(superitem.getType());
            sb.append(",");
            
        }
        String str = sb.substring(0, sb.length() - 1);
        return str;
    }
    public List<SuperItem> getSuperItemsOfPattern(SequenceDatabase database){
        List<SuperItem> superItemsList = new ArrayList<>();
        for (pseudoSuperItem item : superitems){
            SuperItem superitem = item.getSuperItem(database);
            superItemsList.add(superitem);
        }
        return superItemsList;
    }
    
    public boolean splitAlignLink(SequenceDatabase database, pseudoSequentialPattern matePattern, int[] expectedCoords){
        boolean isLinked = false;               
        List<pseudoSuperItem> mateMsSuperItems = matePattern.getMSSuperItems();
        for (int i = 0; i < mateMsSuperItems.size(); i++){
            SuperItem superItem = mateMsSuperItems.get(i).getSuperItem(database);
            if (superItem.getPos() == expectedCoords[1] || superItem.getPos() == expectedCoords[0]){
                isLinked = true;
                break;
            }
        }               
        return isLinked;
    }
    
    public List<svOutInfo> doLocalAlign(SequenceDatabase database, String refSeqString, int refRegionLeft){
                
        stringMatcher strMatcher = new stringMatcher();
        List<String> mForStrings = new ArrayList<>(patternLength);
        StringBuilder sb;
        for (int i = 0; i < patternLength; i++){
            
            SuperItem superItem = superitems.get(i).getSuperItem(database);            
            sb = new StringBuilder();            
            if (superItem.getType().equals("MS")){
                sb.append(superItem.getMachedConsensus());
                sb.append(superItem.getClippedConsensus());
            }
            else if (superItem.getType().equals("SM")){
                sb.append(superItem.getClippedConsensus());
                sb.append(superItem.getMachedConsensus());
            }                        
            
            mForStrings.add(sb.toString());            
        }        
//        int[] alignCoords = new int[2];
        List<svOutInfo> svFromLocalAlign = new ArrayList<>();
        strMatcher.doStringAlign(mForStrings, superitems, refRegionLeft, refSeqString, svFromLocalAlign, database);               
        
        return svFromLocalAlign;
    }
    
    
    public int getSplitStatus(SequenceDatabase database){
        int splitStatus = -3;
        
        if (splitAlignLeftMatch >= 0 && splitAlignRightMatch >= 0){
            splitStatus = -2;
        }
        else if(splitAlignCoords[0] >= patternLeftMostPos && splitAlignCoords[1] <= patternRightMostPos){
            splitStatus = -2;
            if (splitAlignLeftMatch >= 0){
                pseudoSuperItem leftPSItem = superitems.get(splitAlignLeftMatch);
                String type = leftPSItem.getSuperItem(database).getType();                
                int closestToSplitAlign = Integer.MAX_VALUE;
                int rightSplitAlign = splitAlignCoords[1];               
                int rightClosestSuperItemPos = -1;
                for (int i = splitAlignLeftMatch; i < patternLength; i++){
                    pseudoSuperItem psItem = superitems.get(i);
                    SuperItem superItem = psItem.getSuperItem(database);
                    if (superItem.getType().equals(type)){
                        continue;
                    }
                    int dist = Math.abs(superItem.getPos() - rightSplitAlign);
                    if (dist < closestToSplitAlign){
                        closestToSplitAlign = dist;                        
                        rightClosestSuperItemPos = superItem.getPos();
                    }
                }
                splitAlignCoords[1] = rightClosestSuperItemPos;
            }else if (splitAlignRightMatch >= 0){
                pseudoSuperItem rightPSItem = superitems.get(splitAlignRightMatch);
                String type = rightPSItem.getSuperItem(database).getType();
                int closestToSplitAlign = Integer.MAX_VALUE;
                int leftSplitAlign = splitAlignCoords[0];               
                int leftClosestSuperItemPos = -1;
                for (int i = splitAlignRightMatch; i > 0; i--){
                    pseudoSuperItem psItem = superitems.get(i);
                    SuperItem superItem = psItem.getSuperItem(database);
                    if (superItem.getType().equals(type)){
                        continue;
                    }
                    int dist = Math.abs(superItem.getPos() - leftSplitAlign);
                    if (dist < closestToSplitAlign){
                        closestToSplitAlign = dist;                        
                        leftClosestSuperItemPos = superItem.getPos();
                    }
                }
                splitAlignCoords[0] = leftClosestSuperItemPos;
            }
        }
        else if(splitAlignMate > 0){
            splitStatus = splitAlignMate;
        }else{
            splitStatus = -1;
            if (splitAlignCoords[0] < patternLeftMostPos || splitAlignCoords[1] > patternRightMostPos){
                splitStatus = -3;
            }
        }
        return splitStatus;
    }
    
    public void splitAlignCheck(SequenceDatabase database, List<pseudoSequentialPattern> patternCandidates, Map<Integer, Integer> patternIdx, int[] splitAlignPos){                 
//        boolean leftBPMatched = false;
//        boolean rightBPMatched = false;
        
        // Link status is set to -1 by default, which indicates none linkage.
//        int linkStatus = -1;        
        
        splitAlignCoords = splitAlignPos;
        for (int i = 0; i < superitems.size(); i++){            
            pseudoSuperItem psSuperItem = superitems.get(i);           
            SuperItem superItem = psSuperItem.getSuperItem(database);
            if (superItem.getPos() == splitAlignCoords[0]){
//                leftBPMatched = true;
//                break;
                splitAlignLeftMatch = i;
            }            
        }
        for (int i = superitems.size() - 1; i >= 0; i--){
            pseudoSuperItem psSuperItem = superitems.get(i);
            SuperItem superItem = psSuperItem.getSuperItem(database);

            if (superItem.getPos() == splitAlignCoords[1]){
//                rightBPMatched = true;
//                break;
                splitAlignRightMatch = i;
            }
        }

        if (splitAlignLeftMatch == -1 || splitAlignRightMatch == -1){
            Set<Integer> unlinkablePatternIdx = patternIdx.keySet();
            for (Integer idx : unlinkablePatternIdx){
                pseudoSequentialPattern matePattern = patternCandidates.get(idx);
                
                boolean islinked = splitAlignLink(database, matePattern, splitAlignCoords);
                if (islinked){
                    splitAlignMate = idx;
                    break;
                }
            }
        }
    }
    /**
     * A pattern is self-linked by discordant read-pairs
     * @param database
     * @return 
     */
    public boolean isSelfLinked(SequenceDatabase database){
        // Pattern linked by split alignment
        
        boolean linked = false;
        boolean hasEnoughARPs = false;
        Linkage linker = new Linkage();

        int Arps = ARPSuperItems.size();

        int curSuperItemIdx = -1;
        List<pseudoSuperItem> searchSpace;
        int machtedSuperItemIdx = -1;

        Set<Integer> matchedItem = new HashSet<>();

        for (int i = 0; i < Arps; i++){
            if (matchedItem.contains(i)){
                continue;
            }
            pseudoSuperItem target = ARPSuperItems.get(i);
            Map<Integer, Integer> idxMap = new HashMap<>();
            searchSpace = minusSelf(ARPSuperItems, i, idxMap);

            int mateIndex = linker.mateSuperItemSearch(searchSpace, target, database);
            if (mateIndex != -1){    

                curSuperItemIdx = i;
                machtedSuperItemIdx = mateIndex;
                int originalIdx = idxMap.get(machtedSuperItemIdx);
                matchedItem.add(originalIdx);
                SuperItem superItemOne = ARPSuperItems.get(curSuperItemIdx).getSuperItem(database);
//                    SuperItem superItemTwo = searchSpace.get(machtedSuperItemIdx).getSuperItem(database);
                SuperItem superItemTwo = ARPSuperItems.get(originalIdx).getSuperItem(database);

                boolean isEnoughARPs = linker.supportARPs(superItemOne, superItemTwo);
                numOfLinkedEvidence = linker.getSupLink() > numOfLinkedEvidence ? linker.getSupLink():numOfLinkedEvidence;
                if (isEnoughARPs){
                    
                    hasEnoughARPs = true;
                    int superItemOnePos = superItemOne.getPos();
                    int superItemTwoPos = superItemTwo.getPos();
                    if (superItemTwoPos > superItemOnePos){
                        linkedPatterns[0] = curSuperItemIdx;
                        linkedPatterns[1] = idxMap.get(machtedSuperItemIdx);
                    }else{
                        linkedPatterns[0] = idxMap.get(machtedSuperItemIdx);
                        linkedPatterns[1] = curSuperItemIdx;
                    }
                }
            }
        }
        if (hasEnoughARPs){ 
            linked = true;              
        }
        return linked;
    }
//        private int depthSearch(List<pseudoSuperItem> psSuperItems, int start){
//            
//        }
    private List<pseudoSuperItem> minusSelf(List<pseudoSuperItem> psSuperItems, int idx, Map<Integer, Integer> idxMap){
        List<pseudoSuperItem> newSuperItems= new ArrayList<>();
        int length = psSuperItems.size();
        for (int i = 0; i < length; i ++){
            if (i != idx){
                newSuperItems.add(psSuperItems.get(i));
                idxMap.put(newSuperItems.size() - 1, i);
            }

        }
        return newSuperItems;
    }

    public SuperItem getSuperItemFromOriginal(SequenceDatabase database, int idx){
        pseudoSuperItem item = superitems.get(idx);
        return item.getSuperItem(database);
    }
    public SuperItem getSuperItemOfPatternAtPos(SequenceDatabase database, int idx){
        pseudoSuperItem item = ARPSuperItems.get(idx);
        return item.getSuperItem(database);
    }
    public List<pseudoSuperItem> mergeTwoPattern(pseudoSequentialPattern aPattern, SequenceDatabase database){          

        List<pseudoSuperItem> mergedSuperitems = new ArrayList<>();
        List<SuperItem> patternOneSuperItems = getSuperItemsOfPattern(database);
        List<SuperItem> patternTwoSuperItems = aPattern.getSuperItemsOfPattern(database);

        int lengthOne = patternOneSuperItems.size();
        int lengthTwo = patternTwoSuperItems.size();

        int matchedIndexAtPatternOne = -1;           
        int lastMatchedIndexAtPatternOne = lengthOne; 

        SuperItem patternTwoStartSuperItem = patternTwoSuperItems.get(0);
        SuperItem patternTwoLastSuperItem = patternTwoSuperItems.get(lengthTwo - 1);

        for (int i = 0; i < lengthOne ;i++){
            SuperItem patternOneSuperItem = patternOneSuperItems.get(i);
            if (patternTwoStartSuperItem.isEqual(patternOneSuperItem)){
                matchedIndexAtPatternOne = i;                    
            }
            if (patternTwoLastSuperItem.isEqual(patternOneSuperItem)){
                lastMatchedIndexAtPatternOne = i;
            }
        }

        if (matchedIndexAtPatternOne != - 1){
            if (lastMatchedIndexAtPatternOne == lengthOne){
                List<pseudoSuperItem> subListOfPatternOne = superitems.subList(0, matchedIndexAtPatternOne);
                mergedSuperitems.addAll(subListOfPatternOne);
                mergedSuperitems.addAll(aPattern.superitems);
            }else{
                mergedSuperitems = superitems;
            }                
        }

        return mergedSuperitems;            
    }
    /**
     * Estimate breakpoints of two linked patterns
     * @param matePattern
     * @param database
     * @return 
     */
    public int[] estimateBreakpointPos(pseudoSequentialPattern matePattern, SequenceDatabase database){
        // First check the split alignment info on both side
        int[] tmpCoords = splitAlignForBP(database);
        int[] coords = new int[2];
        
        // Process the patten ahead.
        int rightMostARPsuperitemIdx = -1;
        String ARPsuperItemOri = "";
        int ARPsuperItemPos = 0;
        for (int k = patternLength - 1; k >=0 ;k--){
            SuperItem superitem = getSuperItemFromOriginal(database, k);
            if (superitem.isARPsuperitem()){
                coords[0] = superitem.getPos();
                rightMostARPsuperitemIdx = k;
                ARPsuperItemOri = superitem.getOri();
                ARPsuperItemPos = superitem.getPos();
                break;
            }
        }        
        if (ARPsuperItemOri.equals("+")){
//            int closestBP = 500;
            int closestBP = 100;
            if (rightMostARPsuperitemIdx < patternLength){
                for (int i = rightMostARPsuperitemIdx; i < patternLength; i ++){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){
//                        coords[0] = superitem.getPos();
                        int dist = Math.abs(superitem.getPos() - ARPsuperItemPos);
                        if (dist < closestBP){
                            closestBP = dist;
                            coords[0] = superitem.getPos();
                        }
                        if (coords[0] == tmpCoords[0]){
                            return tmpCoords;
                        }
                        break;
                    }
                }
            }if(rightMostARPsuperitemIdx > 0){
                for (int i = rightMostARPsuperitemIdx; i>=0 ; i--){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){
                        int dist = Math.abs(superitem.getPos() - ARPsuperItemPos);
                        if (dist < closestBP){                            
                            coords[0] = superitem.getPos();
                            if (coords[0] == tmpCoords[0]){
                                return tmpCoords;
                            }
                            break;
                        }
                    }
                }
            }
            
        }else{
//            int closestBP = 500;
            int closestBP = 100;
            if (rightMostARPsuperitemIdx < patternLength){
                for (int i = rightMostARPsuperitemIdx; i < patternLength; i ++){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){
//                        coords[0] = superitem.getPos();
                        int dist = Math.abs(superitem.getPos() - ARPsuperItemPos);
                        if (dist < closestBP){
                            closestBP = dist;
                            coords[0] = superitem.getPos();
                        }
                        if (coords[0] == tmpCoords[0]){
                            return tmpCoords;
                        }
                        break;
                    }
                }
            }if(rightMostARPsuperitemIdx > 0){
                for (int i = rightMostARPsuperitemIdx; i>=0 ; i--){
                    SuperItem superitem = getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){
                        int dist = Math.abs(superitem.getPos() - ARPsuperItemPos);
                        if (dist < closestBP){                            
                            coords[0] = superitem.getPos();
                            if (coords[0] == tmpCoords[0]){
                                return tmpCoords;
                            }
                            break;
                        }
                    }
                }
            }
        }
        
        // Process the mate pattern
        int leftMostARPsuperitemIdx = -1;
        int matePatternLength = matePattern.patternLength;
        String mateARPsuperItemOri = "";
        int mateARPsuperItemPos = 0;
        for (int k = 0; k < matePatternLength; k ++){
            SuperItem superitem = matePattern.getSuperItemFromOriginal(database, k);
            if (superitem.isARPsuperitem()){
                leftMostARPsuperitemIdx = k;
                mateARPsuperItemOri = superitem.getOri();
                coords[1] = superitem.getPos();
                mateARPsuperItemPos = superitem.getPos();
                break;
            }
        }
        if (mateARPsuperItemOri.equals("-")){
//            int closestBP = 500;
            int closestBP = 100;
            if (leftMostARPsuperitemIdx > 0){
                for (int i = leftMostARPsuperitemIdx; i >= 0; i--){
                    SuperItem superitem = matePattern.getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){
//                        coords[1] = superitem.getPos();
                        int dist = Math.abs(mateARPsuperItemPos - superitem.getPos());
                        if (dist < closestBP){
                            closestBP = dist;
                            coords[1] = superitem.getPos();
                        }
                        if (coords[1] == tmpCoords[1]){
                            return tmpCoords;
                        }
                        break;
                    }
                }
            }
            if (leftMostARPsuperitemIdx < matePatternLength){
                for (int i = leftMostARPsuperitemIdx; i < matePatternLength; i++){
                    SuperItem superitem = matePattern.getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("SM")){
                        int dist = Math.abs(mateARPsuperItemPos - superitem.getPos());
                        if (dist < closestBP){
                            coords[1] = superitem.getPos();
                            if (coords[1] == tmpCoords[1]){
                                return tmpCoords;
                            }
                            break;
                        }
                        
                    }
                }
            }
        }else{
//            int closestBP = 500;
            int closestBP = 100;
            if (leftMostARPsuperitemIdx > 0){
                for (int i = leftMostARPsuperitemIdx; i >= 0; i--){
                    SuperItem superitem = matePattern.getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){
//                        coords[1] = superitem.getPos();
                        int dist = Math.abs(mateARPsuperItemPos - superitem.getPos());
                        if (dist < closestBP){
                            closestBP = dist;
                            coords[1] = superitem.getPos();
                        }
                        if (coords[1] == tmpCoords[1]){
                            return tmpCoords;
                        }
                        break;
                    }
                }
            }
            if (leftMostARPsuperitemIdx < matePatternLength){
                for (int i = leftMostARPsuperitemIdx; i < matePatternLength; i++){
                    SuperItem superitem = matePattern.getSuperItemFromOriginal(database, i);
                    if (!superitem.isARPsuperitem() && superitem.getType().equals("MS")){
                        int dist = Math.abs(mateARPsuperItemPos - superitem.getPos());
                        if (dist < closestBP){
                            coords[1] = superitem.getPos();
                            if (coords[1] == tmpCoords[1]){
                                return tmpCoords;
                            }
                            break;
                        }                        
                    }
                }
            }
        }
        if (coords[1] < coords[0]){
            int tmp = coords[0];
            coords[0] = coords[1];
            coords[1] = tmp;
        }
        if (ARPsuperItemPos == tmpCoords[0] || mateARPsuperItemPos == tmpCoords[1]){
            coords = tmpCoords;
        }
        return coords;
    }
//    private int getARPpatternIdx(SuperItem superitem, SequenceDatabase database){
//        int idx = -1;
//        for (int i = 0; i < patternLength; i++){
//            SuperItem curSuperItem = getSuperItemFromOriginal(database, i);
//            if (curSuperItem.isEqual(superitem)){
//                idx = i;
//                break;
//            }
//        }
//        return idx;
//    }
    /**
     * Patterns do not have link info.
     * @param database
     * @return 
     */
    public int[] arpLinkInfer(SequenceDatabase database, int minQ){
        // Use split-alignment to get the BP pos if there exist split align of a SuperItem.         
        int[] pos = new int[]{-1,-1, -1};
        int arpsNum = ARPSuperItems.size();
        
        if (arpsNum > 0){
            pseudoSuperItem psItem = ARPSuperItems.get(arpsNum - 1);
            SuperItem superItem = psItem.getSuperItem(database);
            
            if (superItem.getSumMapQ() < superItem.getWeight() * minQ){
                return pos;
            }
            
            pos[0] = superItem.getPos();
            pos[2] = superItem.getWeight();
            QueryInterval interval = superItem.getSuperitemMateRegion();
            if (superItem.getOri().equals("+")){
                pos[1] = interval.start;
            }else{
                pos[1] = interval.end;                
            }
        }
        if (pos[1] < pos[0]){
            int tmp = pos[0];
            pos[0] = pos[1];
            pos[1] = tmp;
        }
        
        return pos;
    }   
    
    public int[] crossLinkBpEstimate(SequenceDatabase database){    
        // info[0], info[1] start and end pos. 
        // info[2] shared string length
        // info[3] number of reads share string
        int[] info = new int[4];
        
        StringBuilder sb;
        stringMatcher strMatcher = new stringMatcher();
        List<String> mStrings = new ArrayList<>();
        List<String> sForwardStrings = new ArrayList<>();
        List<String> sReverseStrings = new ArrayList<>();

        int allStrLength = 0;
        for (pseudoSuperItem psItem : superitems){
            SuperItem superItem = psItem.getSuperItem(database);
            String matchedStr = superItem.getMachedConsensus();
            String clippedStr = superItem.getClippedConsensus();
            sb = new StringBuilder(clippedStr);
            String clippedRevereStr = sb.reverse().toString();

            allStrLength += matchedStr.length();
            allStrLength += clippedStr.length();
            allStrLength += clippedRevereStr.length();

            mStrings.add(matchedStr);
            sForwardStrings.add(clippedStr);
            sReverseStrings.add(clippedRevereStr);                
        }        
        strMatcher.patternGrowth(mStrings, sForwardStrings, sReverseStrings);
                
        int numOfStrs = superitems.size();
        double avgLen = (double)allStrLength/superitems.size();
        double val = Math.log(numOfStrs * avgLen) / Math.log(4);
        int expect = (int) Math.ceil(val);
        
        int[] linkInfo = strMatcher.isCrossLinked();
        int infoSum = linkInfo[0] + linkInfo[1] + linkInfo[2];
        
        if (infoSum > 0){     
            int observedLen = strMatcher.isForwardExist() ? strMatcher.getForwardSharedStrLength():strMatcher.getReverseSharedStrLength();
            
            if (observedLen > expect){
//                System.out.println("Expect: " + expect + "\tObserve: " + observedLen); 
//                strMatcher.printSharedString(superitems, database);
//                System.out.println("\n");
                List<stringIdentity> identitys = strMatcher.isForwardExist()?strMatcher.getForwardSharedStrIdentitys():strMatcher.getReverseSharedStrIdentitys();
                if (validClipLink(identitys, database)){
                    int[] coords = strMatcher.getEstimateBp(superitems, database);
                    if (coords[0] != coords[1]){
                        info[0] = coords[0];
                        info[1] = coords[1];
                        info[2] = observedLen;
                        info[3] = infoSum;
                    }
                    
                }

            }            
        }
        crossSupInfo = info;
        return info;
    }
    
    private boolean validClipLink(List<stringIdentity> strIdentitys, SequenceDatabase database){
        List<String> types = new ArrayList<>();
        boolean isValid = false;
        Collections.sort(strIdentitys, new Comparator<stringIdentity>(){
        @Override
            public int compare(stringIdentity o1, stringIdentity o2){
                return o1.getSeqId() - o2.getSeqId();
            }        
        });
        
        for (stringIdentity id : strIdentitys){
            pseudoSuperItem psItem = superitems.get(id.getSeqId());
            SuperItem superItem = psItem.getSuperItem(database);
            types.add(superItem.getType());
        }
        if (types.size() > 1){            
            String firstType = types.get(0);
            String lastType = types.get(types.size() - 1);
            if ((firstType.equals("MS") && lastType.equals("SM")) || (firstType.equals("SM") && lastType.equals("MS"))){
                isValid = true;
            }
            isValid = true;
        }        
        return isValid;
    }
    
    public int[] splitAlignForBP(SequenceDatabase database){
        int[] pos = new int[2];
                        
        for (int i = 0;i < superitems.size();i++){
            SuperItem superitem = superitems.get(i).getSuperItem(database);
            int superItemPos = superitem.getPos();
            int splitAlignPos = superitem.getSplitAlignPos();
            if (splitAlignPos != -1 && superItemPos != splitAlignPos){
                
                if (superItemPos < splitAlignPos){
                    pos[0] = superItemPos;
                    pos[1] = splitAlignPos;
                }else{
                    pos[0] = splitAlignPos;
                    pos[1] = superItemPos;
                }
                // If we find the split align within the pattern defined range, return.
                if (pos[0] >= patternLeftMostPos && pos[1] <= patternRightMostPos){
                    splitReadSup = superitem.getSplitReadCount();
                    splitReadMapQ += superitem.getSplitReadMapQ();
                    break;
                }                
            }
            
        }
        return pos;
    }
    
    public boolean selfLinkedSuperItemMapQCheck(int minMapQ){
        boolean isHighQ = true;
        
        int leftExpectedMin = selfLinkedSuperItemWeight[0] * minMapQ;
        int rightExpectedMin = selfLinkedSuperItemWeight[1] * minMapQ;
        
        if (selfLinkedSuperItemMapQ[0] < leftExpectedMin || selfLinkedSuperItemMapQ[1] < rightExpectedMin){
            isHighQ = false;
        }
        
        return isHighQ;
    }
    
    /**
     * Estimate breakpoints of a self-linked pattern
     * @param database
     * @return 
     */
    public int[] selfLinkedPatternBP(SequenceDatabase database){
        int[] pos = new int[]{-1,-1};
//        int[] splitAlignPos = splitAlignForBP(database);
        
        SuperItem leftARPItem = ARPSuperItems.get(linkedPatterns[0]).getSuperItem(database);
        int leftARPpatternSuperItemPos = leftARPItem.getPos();
        SuperItem rigthARPItem = ARPSuperItems.get(linkedPatterns[1]).getSuperItem(database);
        int rightARPpatternSuperItemPos = rigthARPItem.getPos();
        
        selfLinkedSuperItemMapQ = new int[2];
        selfLinkedSuperItemMapQ[0] = leftARPItem.getSumMapQ();
        selfLinkedSuperItemMapQ[1] = rigthARPItem.getSumMapQ();
        
        selfLinkedSuperItemWeight = new int[2];
        selfLinkedSuperItemWeight[0] = leftARPItem.getWeight();
        selfLinkedSuperItemWeight[1] = rigthARPItem.getWeight();
        
        int leftClosestToARP = -1;
        int rightClosestToARP = -1;
        int maxDistToLeftARP = 500;
        int maxDistToRightARP = 500;
        
        // Get the closest clipped SuperItem if it exist.
        for (int i = 0; i < patternLength; i++){
            SuperItem superitem = getSuperItemFromOriginal(database, i);
            
            boolean isARPSuperItem = superitem.isARPsuperitem();
            String SuperItemType = superitem.getType();
            int SuperItemPos = superitem.getPos();
            
            if (!isARPSuperItem && SuperItemType.equals("MS")){
                int dist = Math.abs(SuperItemPos - leftARPpatternSuperItemPos);
                if (dist < maxDistToLeftARP){
                    maxDistToLeftARP = dist;
                    leftClosestToARP = i;
                }
            }
            if (!isARPSuperItem && SuperItemType.equals("SM")){
                int dist = Math.abs(SuperItemPos - rightARPpatternSuperItemPos);
                if (dist < maxDistToRightARP){
                    maxDistToRightARP = dist;
                    rightClosestToARP = i;
                }
            }
        }
        
        if (leftClosestToARP != -1 && rightClosestToARP != -1){
            pos[0] = superitems.get(leftClosestToARP).getSuperItem(database).getPos();
            pos[1] = superitems.get(rightClosestToARP).getSuperItem(database).getPos();
        }        
        else if (leftClosestToARP != -1){            
            SuperItem leftSuperItem = getSuperItemFromOriginal(database, leftClosestToARP);
            pos[0] = leftSuperItem.getPos();
            pos[1] = leftSuperItem.getSplitAlignPos();
            // Left clipped SuperItem dose not have supplementary alignments.
            if (pos[1] == -1 && rightClosestToARP > leftClosestToARP){
                SuperItem rightSuperItem = getSuperItemFromOriginal(database, rightClosestToARP);                
                pos[1] = rightSuperItem.getPos();
                int rightSuperItemSplitPos = rightSuperItem.getSplitAlignPos();
                
                // Refine the position if right clipped SuperItem has supplementary alignments and left
                // position is not decided yet.
                if (rightSuperItemSplitPos != -1 && pos[0] == 0){
                    pos[0] = rightSuperItemSplitPos;
                }
            }else if (rigthARPItem.getWeight() > 5 && pos[1] == -1){
                pos[1] = rightARPpatternSuperItemPos;
            }
        }
        else if (rightClosestToARP != -1){
            SuperItem rightSuperItem = getSuperItemFromOriginal(database, rightClosestToARP);
            pos[1] = rightSuperItem.getPos();
            pos[0] = rightSuperItem.getSplitAlignPos();
            if (pos[0] == -1 && (leftClosestToARP != -1 && leftClosestToARP < rightClosestToARP)){
                SuperItem leftSuperItem = getSuperItemFromOriginal(database, leftClosestToARP);
                pos[0] = leftSuperItem.getPos();
                int leftSuperItemSplitPos = leftSuperItem.getSplitAlignPos();
                
                if (leftSuperItemSplitPos != -1 && pos[1] == 0){
                    pos[1] = leftSuperItemSplitPos;
                }
                
            }else if (leftARPItem.getWeight() > 5 && pos[0] == -1){
                pos[0] = leftARPpatternSuperItemPos;
            }
        }
        
        if (pos[0] == -1 || pos[1] == -1){
            
            pos[0] = leftARPpatternSuperItemPos;
            pos[1] = rightARPpatternSuperItemPos;                        
        }      
        
        return pos;
        
    }        
    
}
